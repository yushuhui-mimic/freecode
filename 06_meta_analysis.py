import os, numpy as np, pandas as pd, statsmodels.api as sm

OUT = r"./output"
BASE = r"./data/pooled_core.parquet"

cl = pd.read_csv(os.path.join(OUT, "lab_cluster_assignment.csv"))
core = pd.read_parquet(BASE)[["source", "patient_id", "stay_id", "age", "male"]]
for c in ["patient_id","stay_id"]:
    cl[c] = cl[c].astype(str); core[c] = core[c].astype(str)
d = cl.merge(core, on=["source", "patient_id", "stay_id"], how="left")
d["hosp_mort"] = pd.to_numeric(d["hosp_mort"], errors="coerce")
d = d[d.hosp_mort.notna()].copy()
d["age"] = pd.to_numeric(d["age"], errors="coerce")
print("merged", d.shape, d.source.value_counts().to_dict())

rows = []
for src, g in d.groupby("source"):
    g = g.dropna(subset=["age"])
    # 未校正
    for c in sorted(g.cluster.unique()):
        if c == 1: continue
        sub = g[g.cluster == c]; base = g[g.cluster == 1]
        if len(sub) < 20: continue
        y = pd.concat([sub.hosp_mort, base.hosp_mort])
        X = pd.DataFrame({"const": 1.0,
                          "cl": [1]*len(sub) + [0]*len(base)})
        X["age"] = pd.concat([sub.age, base.age]).values
        X["male"] = pd.concat([pd.to_numeric(sub.male, errors="coerce"),
                               pd.to_numeric(base.male, errors="coerce")]).values
        X = X.fillna(X.mean())
        m = sm.Logit(y.values, X.values).fit(disp=False)
        b, se = m.params[1], m.bse[1]
        rows.append(dict(source=src, cluster=c, n=len(sub), ref_n=len(base),
                         mort=round(sub.hosp_mort.mean()*100, 1),
                         ref_mort=round(base.hosp_mort.mean()*100, 1),
                         OR=round(np.exp(b), 2), lci=round(np.exp(b-1.96*se), 2),
                         uci=round(np.exp(b+1.96*se), 2), logOR=b, se=se))
R = pd.DataFrame(rows)
print("\n=== 各库校正 OR（vs 簇1，校正年龄性别） ===")
print(R[["source","cluster","n","mort","ref_mort","OR","lci","uci"]].to_string(index=False))

# DerSimonian-Laird 随机效应
def dl(sub):
    y, s = sub.logOR.values, sub.se.values
    w = 1/s**2
    fe = (w*y).sum()/w.sum()
    Q = (w*(y-fe)**2).sum()
    dfree = len(y)-1
    C = w.sum() - (w**2).sum()/w.sum()
    tau2 = max(0, (Q-dfree)/C) if C > 0 else 0
    wr = 1/(s**2+tau2)
    est = (wr*y).sum()/wr.sum()
    se = np.sqrt(1/wr.sum())
    I2 = max(0, (Q-dfree)/Q*100) if Q > 0 else 0
    return est, se, tau2, I2, Q, dfree

print("\n=== 随机效应合并（vs 簇1） ===")
for c in sorted(R.cluster.unique()):
    sub = R[R.cluster == c]
    if len(sub) < 2: continue
    est, se, tau2, I2, Q, dfree = dl(sub)
    print(f"簇{c}: k={len(sub)}  OR={np.exp(est):.2f} ({np.exp(est-1.96*se):.2f}-{np.exp(est+1.96*se):.2f})  "
          f"I2={I2:.1f}%  Q={Q:.1f}(df={dfree})  tau2={tau2:.3f}")

R.to_csv(os.path.join(OUT, "adjusted_OR_meta_input.csv"), index=False, encoding="utf-8-sig")

# 每库每簇死亡率表（宽表）
tab = d.groupby(["source","cluster"]).hosp_mort.agg(["size","mean"]).reset_index()
tab["mort%"] = (tab["mean"]*100).round(1)
print("\n=== 各库每簇死亡率(%) ===")
print(tab.pivot(index="cluster", columns="source", values="mort%").to_string())
print("\n=== 各库每簇人数 ===")
print(tab.pivot(index="cluster", columns="source", values="size").fillna(0).astype(int).to_string())
print("\nWROTE", OUT)
