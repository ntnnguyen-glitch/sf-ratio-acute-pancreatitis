## Install the R packages used by R/04_analysis.R and R/05_figures.R
pkgs <- c("survival", "pROC", "WeightIt", "MatchIt", "cobalt", "mice", "EValue", "sandwich", "lmtest", "jsonlite",
          "ggplot2", "survminer", "forestploter", "gridExtra", "ragg")
install.packages(setdiff(pkgs, rownames(installed.packages())))
