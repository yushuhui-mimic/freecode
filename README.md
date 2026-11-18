# Transportability of a six-laboratory sepsis phenotype across five critical care databases

Analysis code accompanying the manuscript:

> **Transportability of a simple six-laboratory sepsis phenotype across five critical care
> databases: a multicohort retrospective study**

## What this repository contains

Reproducible analysis code only. **No patient-level data are included**, because the data use
agreements of all five source databases prohibit redistribution of row-level records.

| File | Purpose |
|---|---|
| `01_pooled_data_qc.py` | Missingness, distributions and unit checks across the five databases |
| `02_cluster_k4_exploratory.py` | Exploratory k-means solution (k = 4, complete cases) |
| `03_cluster_k2_mice.py` | Main analysis: within-cohort MICE, z-scoring, k-means (k = 2) |
| `04_sensitivity.py` | Six sensitivity analyses (exclusion, primary set, winsorising, k = 3/4) |
| `05_five_vs_six_labs.py` | Five-variable model (lactate excluded) and adjusted Rand index |
| `06_meta_analysis.py` | Per-database logistic regression and DerSimonian–Laird pooling |
| `07_table1.py` | Baseline characteristics (Table 1 and Table 2) |
| `08_latent_class_validation.R` | Latent class analysis and agreement with k-means |
| `09_figures.R` | All manuscript figures (Arial, Okabe-Ito, 300 dpi) |
| `10_sae_cohort_sensitivity.R` | Sepsis-associated encephalopathy cohort sensitivity analyses |
| `11_hrv_mortality.R` | Heart-rate variability and mortality analyses |

## Data sources and access

The analysis uses five de-identified critical care databases, each available to credentialed
users under its own data use agreement:

| Database | Access | Version |
|---|---|---|
| MIMIC-IV | https://physionet.org/content/mimiciv/ | 3.1 |
| eICU-CRD | https://physionet.org/content/eicu-crd/ | 2.0 |
| SICdb | https://physionet.org/content/sicdb/ | 1.0.8 |
| AmsterdamUMCdb | https://doi.org/10.17026/dans-22u-f8vd | 1.0.2 |
| Zigong infection database | https://physionet.org/content/icu-infection-zigong-fourth/1.1/ | 1.1 |

Access requires completion of human-subjects research training (CITI Program, "Data or
Specimens Only Research") and signature of the corresponding data use agreement.

## Input files expected

The scripts expect a harmonised pooled table in Parquet format with one row per sepsis
episode and the following columns (among others):

```
source, patient_id, stay_id, hadm_id, age, male, sofa, apsiii,
creatinine_max, platelets_min, wbc_max, hemoglobin_min, lactate_max, bilirubin_max,
vasoactive_any, vent_any, rrt_any, hosp_mort, los_icu, in_primary, def_level
```

Set the input path at the top of each script (`BASE` / `POOL`) to point at your local copy.
The harmonisation itself was performed from the five source databases in PostgreSQL; the
pooling script is database-specific and is not distributed here because it embeds local
connection details.

## Environment

```
Python 3.14
  pandas 2.3.3
  numpy
  scikit-learn 1.9.0
  statsmodels 0.14.6
  matplotlib 3.11.1
R 4.6.1
  ggplot2, data.table, scales, patchwork, poLCA
```

Install Python dependencies with:

```bash
pip install pandas numpy scikit-learn statsmodels matplotlib
```

## Reproducibility

- All random processes use fixed seeds (`random_state = 11` in Python, `set.seed(42)` in R).
- k-means uses 50 random starts.
- Figures are exported as vector PDF and as 300 dpi PNG.

## Citation

If you use this code, please cite the manuscript (DOI to be added upon publication) and the
Zenodo archive of this repository.

## Licence

MIT Licence for the code. The underlying data remain subject to their original data use
agreements and are **not** covered by this licence.

## Contact

Shuhui Yu, Department of Critical Care Medicine, Renmin Hospital of Wuhan University
(yushuhui@whu.edu.cn, ORCID 0009-0009-8683-0630)
