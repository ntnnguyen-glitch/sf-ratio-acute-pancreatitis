# SpO2/FiO2-defined persistent respiratory failure in critically ill patients with acute pancreatitis

Code for the retrospective cohort study "Persistent respiratory failure defined by the SpO2/FiO2 ratio and mortality in critically ill patients with acute pancreatitis: a retrospective cohort study in two intensive care databases".

The repository contains code only. No patient-level data are included.

## Data access

Both databases are available from PhysioNet after credentialing, completion of a human-subjects research course and signing of the data use agreement:

- MIMIC-IV v3.1: https://doi.org/10.13026/kpb9-mt58
- eICU Collaborative Research Database v2.0: https://doi.org/10.13026/C2WM1R

## Requirements

- Python 3.10 or later with the packages in `requirements.txt`
- R 4.6.1 with the packages listed in `R/packages.R`

## How to run

```bash
export MIMIC_DIR=/path/to/mimic-iv-3.1      # folder with hosp/ and icu/
export EICU_DIR=/path/to/eicu-crd-2.0       # folder with the .csv.gz files
bash run_all.sh
```

`run_all.sh` runs the steps below from the repository root.

| Step | Script | Output |
|---|---|---|
| 1 | `python/01_extract_mimic.py` | `data/mimic.duckdb` |
| 2 | `python/02_extract_eicu.py` | `data/eicu.duckdb` |
| 3 | `python/03_build_variables.py` | `data/analytic.csv`, `data/pairs.csv` |
| 4 | `R/04_analysis.R` | `results/` (all estimates, Table 1, agreement, balance) |
| 5 | `R/05_figures.R` | `figures/Figure1-4` (TIFF and PNG) |

The `data/` folder contains derived patient-level data and must not be shared or uploaded.

## Definitions

- Cohort: adults at their first ICU stay with acute pancreatitis (MIMIC-IV: ICD-10 K85 or ICD-9 577.0 among the first three diagnoses; eICU-CRD: APACHE admission diagnosis or problem-list diagnosis of acute pancreatitis, first ICU stay per patient).
- Respiratory failure: worst SpO2/FiO2 <=315 within 24 h of ICU admission.
- Persistent or transient respiratory failure: worst SpO2/FiO2 at 48-72 h <=315 or >315; patients without a value at 48-72 h are not classifiable.
- Landmark: patients alive 72 h after ICU admission.
- Primary outcome: 28-day mortality (MIMIC-IV) and in-hospital mortality (eICU-CRD), analysed separately.

## License

MIT, see `LICENSE`.
