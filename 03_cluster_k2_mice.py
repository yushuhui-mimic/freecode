import os
import numpy as np
import pandas as pd
from sklearn.experimental import enable_iterative_imputer  # noqa
from sklearn.impute import IterativeImputer
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, adjusted_rand_score

BASE = r"./data/pooled_core.parquet"
OUT = r"./output"

LABS = ["creatinine_max", "platelets_min", "wbc_max",
        "hemoglobin_min", "lactate_max", "bilirubin_max"]
SKEWED = ["creatinine_max", "wbc_max", "lactate_max", "bilirubin_max"]

df = pd.read_parquet(BASE)
df["hosp_mort"] = pd.to_numeric(df["hosp_mort"], errors="coerce")
d = df[df.def_level.isin(["A", "B", "C", "C_flag"]) & df.hosp_mort.notna()].copy()
print("分析集:", d.shape, d.source.value_counts().to_dict())

# 库内预处理
for v in LABS:
    d[v] = pd.to_numeric(d[v], errors="coerce")
    d[v] = d.groupby("source")[v].transform(lambda s: s.clip(s.quantile(.01), s.quantile(.99)))
    if v in SKEWED:
        d[v] = np.log1p(d[v])

# 每库单独 MICE（避免跨库分布污染），至少要有该列的库
parts = []
for src, g in d.groupby("source"):
    g = g.copy()
    cols = [v for v in LABS if g[v].notna().sum() >= 50]
    if "lactate_max" not in cols:
        g["lactate_max"] = np.nan
    if cols:
        imp = IterativeImputer(max_iter=10, random_state=42, sample_posterior=False)
        g[cols] = imp.fit_transform(g[cols])
    # 仍缺的（整列缺失）填 0（=该库均值，z 后）
    g[LABS] = g[LABS].fillna(0)
    parts.append(g)
d = pd.concat(parts, ignore_index=True)

# z 分数
z = d.copy()
for v in LABS:
    z[v] = z.groupby("source")[v].transform(lambda s: (s - s.mean()) / (s.std(ddof=0) or 1))

X = z[LABS].values
print("缺失情况（填补后，按库）：")
print(z.groupby("source")[LABS].apply(lambda g: g.isna().mean()).round(3).to_string())

sid = np.random.default_rng(11).choice(len(X), size=min(40000, len(X)), replace=False)
sil = {}
for k in range(2, 7):
    km = KMeans(n_clusters=k, n_init=10, random_state=11).fit(X[sid])
    sil[k] = round(silhouette_score(X[sid], km.labels_, sample_size=15000, random_state=11), 4)
print("轮廓系数:", sil)

km = KMeans(n_clusters=4, n_init=50, random_state=11).fit(X)
z["cluster"] = km.labels_
ordr = z.groupby("cluster").hosp_mort.mean().sort_values().index.tolist()
z["cluster"] = z.cluster.map({c: i+1 for i, c in enumerate(ordr)})

prof = z.groupby("cluster")[LABS].mean().round(2)
raw = d.copy(); raw["cluster"] = z["cluster"].values
prof_raw = raw.groupby("cluster")[LABS].median().round(2)
overall = z.groupby("cluster").hosp_mort.agg(["size", "mean"])
overall["pct"] = (overall["mean"]*100).round(1)

print("\n=== 簇画像（z） ==="); print(prof.to_string())
print("\n=== 簇画像（原始单位中位数，注意已对数变换的列） ==="); print(prof_raw.to_string())
print("\n=== 总体 ==="); print(overall.to_string())

tab = z.groupby(["source", "cluster"]).hosp_mort.agg(["size", "mean"]).reset_index()
tab["mort%"] = (tab["mean"]*100).round(1)
print("\n=== 各库每簇死亡率(%) ===")
print(tab.pivot(index="cluster", columns="source", values="mort%").to_string())
print("\n=== 各库每簇人数 ===")
print(tab.pivot(index="cluster", columns="source", values="size").fillna(0).astype(int).to_string())

# 与完整病例版的一致性
try:
    old = pd.read_csv(os.path.join(OUT, "lab_cluster_assignment.csv"))
    for c in ["patient_id", "stay_id"]:
        old[c] = old[c].astype(str); z[c] = z[c].astype(str)
    m = z.merge(old[["source","patient_id","stay_id","cluster"]], on=["source","patient_id","stay_id"],
                suffixes=("", "_cc"))
    if len(m) > 100:
        print("\n与完整病例版聚类一致性 ARI =", round(adjusted_rand_score(m.cluster_cc, m.cluster), 3),
              f"(n={len(m)})")
except Exception as e:
    print("ARI 跳过:", e)

z[["source","patient_id","stay_id","cluster","hosp_mort"]+LABS].to_csv(
    os.path.join(OUT, "lab_cluster_assignment_mice.csv"), index=False, encoding="utf-8-sig")
prof.to_csv(os.path.join(OUT, "cluster_profiles_mice.csv"), encoding="utf-8-sig")
overall.to_csv(os.path.join(OUT, "cluster_overall_mice.csv"), encoding="utf-8-sig")
tab.to_csv(os.path.join(OUT, "mortality_by_cluster_source_mice.csv"), index=False, encoding="utf-8-sig")
print("\nWROTE", OUT)
