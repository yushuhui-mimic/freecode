import os, sys, json
import numpy as np
import pandas as pd

BASE = r"./data"
OUT = r"./output/qc"
os.makedirs(OUT, exist_ok=True)

p_parq = os.path.join(BASE, "pooled_core.parquet")
p_csv = os.path.join(BASE, "pooled_core.csv")
try:
    df = pd.read_parquet(p_parq)
    print("read parquet", df.shape)
except Exception as e:
    print("parquet failed:", e, "-> csv")
    df = pd.read_csv(p_csv, low_memory=False)
    print("read csv", df.shape)

num_vars = ["age","male","bmi","sofa","apsiii","sapsii","oasis","charlson",
            "lactate_max","creatinine_max","bun_max","platelets_min","wbc_max",
            "hemoglobin_min","hematocrit_min","bilirubin_max","inr_max","ptt_max","pt_max",
            "fibrinogen_min","albumin_min","pao2_min","fio2","pao2fio2_min","ph_min",
            "glucose_mean","sodium_min","potassium_max","chloride_min","calcium_min",
            "bicarbonate_min","anion_gap_max","d_dimer_max",
            "hr_mean","sbp_mean","mbp_mean","rr_mean","temp_mean","spo2_mean","gcs_min",
            "urine_24h","vasoactive_any","vent_any","rrt_any",
            "hosp_mort","icu_mort","mort_28d","los_icu","los_hosp"]

num_vars = [v for v in num_vars if v in df.columns]
rows = []
for src, g in df.groupby("source"):
    for v in num_vars:
        s = pd.to_numeric(g[v], errors="coerce")
        n = len(s); k = s.notna().sum()
        if k == 0:
            rows.append(dict(source=src, variable=v, n=n, nonmiss=0, pct=0.0,
                             median=np.nan, p25=np.nan, p75=np.nan,
                             p01=np.nan, p99=np.nan, mean=np.nan))
            continue
        rows.append(dict(source=src, variable=v, n=n, nonmiss=int(k), pct=round(100*k/n,1),
                         median=float(s.median()), p25=float(s.quantile(.25)),
                         p75=float(s.quantile(.75)), p01=float(s.quantile(.01)),
                         p99=float(s.quantile(.99)), mean=float(s.mean())))

qc = pd.DataFrame(rows)
qc.to_csv(os.path.join(OUT, "pooled_qc_long.csv"), index=False, encoding="utf-8-sig")

# 可比性矩阵：每个变量有多少库的非缺失率 >= 50%
piv = qc.pivot_table(index="variable", columns="source", values="pct")
piv = piv.round(1)
piv.to_csv(os.path.join(OUT, "pooled_missing_matrix.csv"), encoding="utf-8-sig")

med = qc.pivot_table(index="variable", columns="source", values="median").round(2)
med.to_csv(os.path.join(OUT, "pooled_median_matrix.csv"), encoding="utf-8-sig")

# 判定：中位数跨库比值（最大/最小），用于发现单位不一致
ratio = med.replace(0, np.nan)
r = (ratio.max(axis=1) / ratio.min(axis=1)).round(2)

summary = pd.DataFrame({
    "n_sources": (qc[qc.pct >= 50].groupby("variable").source.nunique()),
    "median_ratio_maxmin": r,
})
summary["comparable"] = np.where(
    (summary.n_sources >= 4) & (summary.median_ratio_maxmin.between(0.5, 2.0)), "YES",
    np.where((summary.n_sources >= 4) & (summary.median_ratio_maxmin.between(0.3, 3.0)), "CHECK",
             np.where(summary.n_sources >= 4, "UNIT_MISMATCH", "TOO_SPARSE")))
summary = summary.sort_values(["comparable","variable"])
summary.to_csv(os.path.join(OUT, "variable_comparability.csv"), encoding="utf-8-sig")

print("\n=== 每个库的样本量与变量非缺失率前 12 ===")
print(qc.groupby("source").n.first().to_string())
print("\n=== 可比性判定 ===")
print(summary.to_string())
print("\nWROTE", OUT)
