import os
import numpy as np
import pandas as pd
from sklearn.experimental import enable_iterative_imputer  # noqa
from sklearn.impute import IterativeImputer
from sklearn.cluster import KMeans

BASE = r"./data/pooled_core.parquet"
OUT = r"./output"

LABS = ["creatinine_max", "platelets_min", "wbc_max", "hemoglobin_min", "lactate_max", "bilirubin_max"]
SKEWED = ["creatinine_max", "wbc_max", "lactate_max", "bilirubin_max"]

raw = pd.read_parquet(BASE)
raw["hosp_mort"] = pd.to_numeric(raw["hosp_mort"], errors="coerce")


def run(q=0.01, exclude_zi=False, primary_only=False, k=2, seed=11, verbose=True):
    d = raw[raw.def_level.isin(["A", "B", "C", "C_flag"]) & raw.hosp_mort.notna()].copy()
    if exclude_zi:
        d = d[d.source != "ziGongCC"]
    if primary_only:
        d = d[pd.to_numeric(d.in_primary, errors="coerce") == 1]
    for v in LABS:
        d[v] = pd.to_numeric(d[v], errors="coerce")
        d[v] = d.groupby("source")[v].transform(lambda s: s.clip(s.quantile(q), s.quantile(1 - q)))
        if v in SKEWED:
            d[v] = np.log1p(d[v])
    parts = []
    for src, g in d.groupby("source"):
        g = g.copy()
        cols = [v for v in LABS if g[v].notna().sum() >= 50]
        if cols:
            g[cols] = IterativeImputer(max_iter=10, random_state=seed).fit_transform(g[cols])
        g[LABS] = g[LABS].fillna(0)
        parts.append(g)
    d = pd.concat(parts, ignore_index=True)
    for v in LABS:
        d[v] = d.groupby("source")[v].transform(lambda s: (s - s.mean()) / (s.std(ddof=0) or 1))

    X = d[LABS].values
    km = KMeans(n_clusters=k, n_init=50, random_state=seed).fit(X)
    d["cl"] = km.labels_
    order = d.groupby("cl").hosp_mort.mean().sort_values().index.tolist()
    d["cl"] = d.cl.map({c: i + 1 for i, c in enumerate(order)})

    rows = []
    for src, g in d.groupby("source"):
        base = g[g.cl == 1]
        for c in range(2, k + 1):
            sub = g[g.cl == c]
            if len(sub) < 50: continue
            a, n1 = sub.hosp_mort.sum(), len(sub)
            b, n0 = base.hosp_mort.sum(), len(base)
            a, n1, b, n0 = a + .5, n1 + 1, b + .5, n0 + 1
            or_ = (a / (n1 - a)) / (b / (n0 - b))
            se = np.sqrt(1 / a + 1 / (n1 - a) + 1 / b + 1 / (n0 - b))
            rows.append(dict(source=src, cluster=c, n=len(sub), mort=sub.hosp_mort.mean(),
                             OR=or_, lci=np.exp(np.log(or_) - 1.96 * se),
                             uci=np.exp(np.log(or_) + 1.96 * se), logOR=np.log(or_), se=se))
    R = pd.DataFrame(rows)

    def dl(sub):
        y, s = sub.logOR.values, sub.se.values
        w = 1 / s ** 2
        fe = (w * y).sum() / w.sum()
        Q = (w * (y - fe) ** 2).sum()
        dfree = len(y) - 1
        C = w.sum() - (w ** 2).sum() / w.sum()
        tau2 = max(0, (Q - dfree) / C) if C > 0 else 0
        wr = 1 / (s ** 2 + tau2)
        est = (wr * y).sum() / wr.sum()
        se = np.sqrt(1 / wr.sum())
        I2 = max(0, (Q - dfree) / Q * 100) if Q > 0 else 0
        return est, se, I2

    res = []
    for c in range(2, k + 1):
        sub = R[R.cluster == c]
        if len(sub) < 2: continue
        est, se, I2 = dl(sub)
        res.append(dict(cluster=c, k_sources=len(sub), pooled_OR=np.exp(est),
                        lci=np.exp(est - 1.96 * se), uci=np.exp(est + 1.96 * se), I2=I2,
                        n_total=int(sub.n.sum())))
    out = dict(n=len(d), k=k, q=q, exclude_zi=exclude_zi, primary_only=primary_only,
               results=pd.DataFrame(res), per_source=R)
    if verbose:
        print(f"\n--- n={len(d)}, k={k}, q={q}, excl_zi={exclude_zi}, primary={primary_only}")
        print(out["results"].round(2).to_string(index=False))
    return out


scenarios = [
    ("主分析 k=2", dict(k=2)),
    ("S1 排除 ziGongCC", dict(k=2, exclude_zi=True)),
    ("S2 仅主分析集", dict(k=2, primary_only=True)),
    ("S3 缩尾 0.5%/99.5%", dict(k=2, q=0.005)),
    ("S4 k=3", dict(k=3)),
    ("S5 k=4", dict(k=4)),
]

summary = []
store = {}
for name, kw in scenarios:
    o = run(**kw)
    store[name] = o
    r = o["results"]
    top = r.iloc[-1]  # highest-mortality cluster
    summary.append(dict(scenario=name, n=o["n"], k=o["k"], cluster=int(top.cluster),
                        pooled_OR=round(top.pooled_OR, 2), lci=round(top.lci, 2),
                        uci=round(top.uci, 2), I2=round(top.I2, 1)))

S = pd.DataFrame(summary)
print("\n================ 敏感性分析汇总 ================")
print(S.to_string(index=False))
S.to_csv(os.path.join(OUT, "sensitivity_summary.csv"), index=False, encoding="utf-8-sig")
store["主分析 k=2"]["per_source"].to_csv(os.path.join(OUT, "per_source_OR_k2.csv"), index=False, encoding="utf-8-sig")
print("\nWROTE", OUT)
