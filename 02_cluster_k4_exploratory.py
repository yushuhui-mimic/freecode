import os
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler

BASE = r"./data/pooled_core.parquet"
OUT = r"./output"
os.makedirs(OUT, exist_ok=True)

LABS = ["creatinine_max", "platelets_min", "wbc_max",
        "hemoglobin_min", "lactate_max", "bilirubin_max"]
SKEWED = ["creatinine_max", "wbc_max", "lactate_max", "bilirubin_max"]

df = pd.read_parquet(BASE)
print("loaded", df.shape)

def prep(g):
    g = g.copy()
    for v in LABS:
        s = pd.to_numeric(g[v], errors="coerce")
        lo, hi = s.quantile(.01), s.quantile(.99)
        s = s.clip(lo, hi)
        if v in SKEWED:
            s = np.log1p(s)
        g[v] = s
    return g

df = df.groupby("source", group_keys=False).apply(prep)

print("\n缺失率（按库）：")
miss = df.groupby("source")[LABS].apply(lambda g: (g.isna().mean()*100).round(1))
print(miss.to_string())

# 主分析集：def_level A/B/C（排除 D/E 弱定义）
main = df[df.def_level.isin(["A", "B", "C", "C_flag"])].copy()
main = main[main.hosp_mort.notna()]
print("\n主分析集（A/B/C）:", main.shape, main.source.value_counts().to_dict())

# 完整病例
cc = main.dropna(subset=LABS)
print("完整病例:", cc.shape, cc.source.value_counts().to_dict())

# 库内标准化
z = cc.copy()
for v in LABS:
    z[v] = z.groupby("source")[v].transform(lambda x: (x - x.mean()) / x.std(ddof=0))

X = z[LABS].values

# 选 k
sample_idx = np.random.default_rng(42).choice(len(X), size=min(30000, len(X)), replace=False)
sil = {}
for k in range(2, 7):
    km = KMeans(n_clusters=k, n_init=10, random_state=42).fit(X[sample_idx])
    sil[k] = silhouette_score(X[sample_idx], km.labels_, sample_size=10000, random_state=42)
print("\n轮廓系数:", {k: round(v, 4) for k, v in sil.items()})
best_k = 4  # 临床可解释性优先，见 silhouette 表
print("选择 k =", best_k)

km = KMeans(n_clusters=best_k, n_init=50, random_state=42).fit(X)
z["cluster"] = km.labels_
# 按死亡率重排：1 = 最轻
ordr = z.groupby("cluster").hosp_mort.mean().sort_values().index.tolist()
z["cluster"] = z.cluster.map({c:i+1 for i,c in enumerate(ordr)})

prof = z.groupby("cluster")[LABS].mean().round(2)
sizes = z.cluster.value_counts().sort_index()
mort = z.groupby(["source", "cluster"]).hosp_mort.agg(["size", "mean"]).round(3)
overall = z.groupby("cluster").hosp_mort.agg(["size", "mean"]).round(3)
overall["pct"] = (overall["mean"] * 100).round(1)

print("\n=== 簇画像（库内 z 分数均值） ===")
print(prof.to_string())
print("\n=== 各库每簇死亡率 ===")
print((mort.reset_index().pivot(index="cluster", columns="source", values="mean") * 100).round(1).to_string())
print("\n=== 每簇样本量（各库） ===")
print(mort.reset_index().pivot(index="cluster", columns="source", values="size").fillna(0).astype(int).to_string())
print("\n=== 总体 ===")
print(overall.to_string())

out = z[["source", "patient_id", "stay_id", "cluster", "hosp_mort"] + LABS].copy()
out.to_csv(os.path.join(OUT, "lab_cluster_assignment.csv"), index=False, encoding="utf-8-sig")
prof.to_csv(os.path.join(OUT, "cluster_profiles.csv"), encoding="utf-8-sig")
mort.reset_index().to_csv(os.path.join(OUT, "mortality_by_cluster_source.csv"), index=False, encoding="utf-8-sig")
pd.Series(sil).to_csv(os.path.join(OUT, "silhouette_by_k.csv"))
print("\nWROTE", OUT)
