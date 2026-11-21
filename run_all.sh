#!/usr/bin/env bash
set -euo pipefail
: "${MIMIC_DIR:?Set MIMIC_DIR}"
: "${EICU_DIR:?Set EICU_DIR}"
python3 python/01_extract_mimic.py
python3 python/02_extract_eicu.py
python3 python/03_build_variables.py
Rscript R/04_analysis.R
Rscript R/05_figures.R
