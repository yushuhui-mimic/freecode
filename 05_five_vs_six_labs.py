"""5 指标（去乳酸）敏感性分析：量化"乳酸不可得"对表型的影响。"""
import os
import numpy as np
import pandas as pd
from sklearn.experimental import enable_iterative_imputer  # noqa
from sklearn.impute import IterativeImputer
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, silhouette_score

BASE = r"./data/pooled_core.parquet"
OUT = r"./output"

SIX = ["creatinine_max","platelets_min","wbc_max","hemoglobin_min","lactate_max","bilirubin_max"]
FIVE = ["creatinine_max","platelets_min","wbc_max","hemoglobin_min","bilirubin_max"]
SKEWED = ["creatinine_max","wbc_max","lactate_max","bilirubin_max"]


def pipeline(labs, k=2, seed=11):
    d = pd.read_parquet(BASE)
    d["hosp_mort"] = pd.to_numeric(d["hosp_mort"], errors="coerce")
    d = d[d.def_level.isin(["A","B","C","C_flag"]) & d.hosp_mort.notna()].copy()
    for v in labs:
        d[v] = pd.to_numeric(d[v], errors="coerce")
        d[v] = d.groupby("source")[v].transform(lambda s: s.clip(s.quantile(.01), s.quantile(.99)))
        if v in SKEWED:
            d[v] = np.log1p(d[v])
    parts = []
    for src, g in d.groupby("source"):
        g = g.copy()
        cols = [v for v in labs if g[v].notna().sum() >= 50]
        if cols:
            g[cols] = IterativeImputer(max_iter=10, random_state=seed).fit_transform(g[cols])
        g[labs] = g[labs].fillna(0)
        parts.append(g)
    d = pd.concat(parts, ignore_index=True)
    for v in labs:
        d[v] = d.groupby("source")[v].transform(lambda s: (s - s.mean()) / (s.std(ddof=0) or 1))
    X = d[labs].values
    km = KMeans(n_clusters=k, n_init=50, random_state=seed).fit(X)
    d["cl"] = km.labels_
    order = d.groupby("cl").hosp_mort.mean().sort_values().index.tolist()
    d["cl"] = d.cl.map({c: i+1 for i, c in enumerate(order)})
    return d, X


def per_source_or(d, k=2):
    rows = []
    for src, g in d.groupby("source"):
        base = g[g.cl == 1]
        for c in range(2, k+1):
            sub = g[g.cl == c]
            if len(sub) < 50: continue
            a, n1 = sub.hosp_mort.sum(), len(sub)
            b, n0 = base.hosp_mort.sum(), len(base)
            a, n1, b, n0 = a+.5, n1+1, b+.5, n0+1
            or_ = (a/(n1-a))/(b/(n0-b))
            se = np.sqrt(1/a + 1/(n1-a) + 1/b + 1/(n0-b))
            rows.append(dict(source=src, cluster=c, n=len(sub), mort=round(sub.hosp_mort.mean()*100,1),
                             OR=round(or_,2), lci=round(np.exp(np.log(or_)-1.96*se),2),
                             uci=round(np.exp(np.log(or_)+1.96*se),2), logOR=np.log(or_), se=se))
    return pd.DataFrame(rows)


def dl_meta(R):
    y, s = R.logOR.values, R.se.values
    w = 1/s**2; fe = (w*y).sum()/w.sum(); Q = (w*(y-fe)**2).sum(); dfree = len(y)-1
    C = w.sum() - (w**2).sum()/w.sum(); tau2 = max(0, (Q-dfree)/C) if C > 0 else 0
    wr = 1/(s**2+tau2); est = (wr*y).sum()/wr.sum(); se = np.sqrt(1/wr.sum())
    I2 = max(0, (Q-dfree)/Q*100) if Q > 0 else 0
    return np.exp(est), np.exp(est-1.96*se), np.exp(est+1.96*se), I2, tau2


print("=" * 70)
print("模型 A：6 指标（含乳酸）")
d6, X6 = pipeline(SIX)
R6 = per_source_or(d6)
print(R6[["source","n","mort","OR","lci","uci"]].to_string(index=False))
print("合并 OR = %.2f (%.2f-%.2f), I2=%.1f%%" % dl_meta(R6)[:4])

print("\n" + "=" * 70)
print("模型 B：5 指标（去乳酸）")
d5, X5 = pipeline(FIVE)
R5 = per_source_or(d5)
print(R5[["source","n","mort","OR","lci","uci"]].to_string(index=False))
print("合并 OR = %.2f (%.2f-%.2f), I2=%.1f%%" % dl_meta(R5)[:4])

# 一致性
key = ["source","patient_id","stay_id"]
for c in key:
    d6[c] = d6[c].astype(str); d5[c] = d5[c].astype(str)
m = d6[key + ["cl","hosp_mort"]].merge(d5[key + ["cl"]], on=key, suffixes=("_6","_5"))
print("\n" + "=" * 70)
print("两模型聚类一致性")
print("总体 ARI = %.3f (n=%d)" % (adjusted_rand_score(m.cl_6, m.cl_5), len(m)))
for src, g in m.groupby("source"):
    print("  %-16s ARI=%.3f  n=%d" % (src, adjusted_rand_score(g.cl_6, g.cl_5), len(g)))

# 只在 6 指标里"乳酸真实存在"的病人
has_lac = d6.lactate_max.notna() & (d6.source != "ziGongCC")
print("\n乳酸真实可得人群（MIMIC/SICdb/eICU/Amsterdam）:")
for src, g in m.groupby("source"):
    if src == "ziGongCC": continue
    print("  %-16s ARI=%.3f" % (src, adjusted_rand_score(g.cl_6, g.cl_5)))
print("  ziGongCC（乳酸全缺）ARI=%.3f" % adjusted_rand_score(
    m[m.source=="ziGongCC"].cl_6, m[m.source=="ziGongCC"].cl_5))

out = pd.concat([R6.assign(model="6 labs"), R5.assign(model="5 labs")])
out.to_csv(os.path.join(OUT, "sensitivity_5vs6_labs.csv"), index=False, encoding="utf-8-sig")
d5[key + ["cl","hosp_mort"]].rename(columns={"cl":"cluster_5lab"}).to_csv(
    os.path.join(OUT, "assignment_5lab.csv"), index=False, encoding="utf-8-sig")
print("\nWROTE", OUT)
