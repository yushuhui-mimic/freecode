import os
import numpy as np
import pandas as pd

BASE = r"./data/pooled_core.parquet"
OUT = r"./output"

d = pd.read_parquet(BASE)
d["hosp_mort"] = pd.to_numeric(d["hosp_mort"], errors="coerce")
d = d[d.def_level.isin(["A","B","C","C_flag"]) & d.hosp_mort.notna()].copy()
for v in ["age","male","sofa","apsiii","vent_any","vasoactive_any","rrt_any","los_icu",
          "lactate_max","creatinine_max","platelets_min","bilirubin_max","wbc_max","hemoglobin_min"]:
    d[v] = pd.to_numeric(d[v], errors="coerce")

cl = pd.read_csv(os.path.join(OUT, "k2_assignment.csv"))
for c in ["patient_id","stay_id"]:
    cl[c] = cl[c].astype(str); d[c] = d[c].astype(str)
d = d.merge(cl[["source","patient_id","stay_id","k2"]], on=["source","patient_id","stay_id"], how="left")

def fmt_cont(s):
    k = s.notna().sum()
    if k == 0: return "NA"
    return f"{s.mean():.1f} ± {s.sd if False else s.std():.1f}"

def fmt_bin(s):
    k = s.notna().sum()
    if k == 0: return "NA"
    return f"{int(s.sum())} ({100*s.mean():.1f})"

rows = []
order = ["MIMIC-IV","eICU","SICdb","ziGongCC","AmsterdamUMCdb"]
for src in order:
    g = d[d.source == src]
    rows.append(dict(
        Database=src, n=len(g),
        Age=fmt_cont(g.age), Male=fmt_bin(g.male),
        SOFA=fmt_cont(g.sofa), APACHE_III=fmt_cont(g.apsiii),
        Lactate=fmt_cont(g.lactate_max), Creatinine=fmt_cont(g.creatinine_max),
        Bilirubin=fmt_cont(g.bilirubin_max), Platelets=fmt_cont(g.platelets_min),
        Vasoactive=fmt_bin(g.vasoactive_any), Ventilation=fmt_bin(g.vent_any),
        RRT=fmt_bin(g.rrt_any), ICU_LOS=fmt_cont(g.los_icu),
        Mortality=fmt_bin(g.hosp_mort)))
T1 = pd.DataFrame(rows)
print("=== Table 1（按数据库） ===")
print(T1.to_string(index=False))

# 按表型
rows2 = []
for k in [1,2]:
    g = d[d.k2 == k]
    rows2.append(dict(Phenotype=f"Phenotype {k}", n=len(g),
        Age=fmt_cont(g.age), Male=fmt_bin(g.male), SOFA=fmt_cont(g.sofa),
        Lactate=fmt_cont(g.lactate_max), Creatinine=fmt_cont(g.creatinine_max),
        Bilirubin=fmt_cont(g.bilirubin_max), Platelets=fmt_cont(g.platelets_min),
        WBC=fmt_cont(g.wbc_max), Hemoglobin=fmt_cont(g.hemoglobin_min),
        Vasoactive=fmt_bin(g.vasoactive_any), Ventilation=fmt_bin(g.vent_any),
        RRT=fmt_bin(g.rrt_any), ICU_LOS=fmt_cont(g.los_icu), Mortality=fmt_bin(g.hosp_mort)))
T2 = pd.DataFrame(rows2)
print("\n=== 按表型 ===")
print(T2.to_string(index=False))

# 各库表型分布
print("\n=== 各库表型分布 ===")
tab = pd.crosstab(d.source, d.k2, margins=True)
print(tab.to_string())

T1.to_csv(os.path.join(OUT, "Table1_by_database.csv"), index=False, encoding="utf-8-sig")
T2.to_csv(os.path.join(OUT, "Table1_by_phenotype.csv"), index=False, encoding="utf-8-sig")
tab.to_csv(os.path.join(OUT, "phenotype_distribution_by_db.csv"), encoding="utf-8-sig")
print("\nWROTE", OUT)
