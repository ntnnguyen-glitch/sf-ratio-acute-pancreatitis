suppressPackageStartupMessages({library(ggplot2);library(survival);library(survminer);library(forestploter);library(grid);library(jsonlite);library(gridExtra)})
## Figures 1-4. Run after R/04_analysis.R.
dir.create("figures", showWarnings = FALSE)
save_fig <- function(p, name, w, h, draw = FALSE) {
  ragg::agg_tiff(file.path("figures", paste0(name, ".tiff")), width = w, height = h, units = "in", res = 300, compression = "lzw"); if (draw) grid.draw(p) else print(p); dev.off()
  png(file.path("figures", paste0(name, ".png")), width = w, height = h, units = "in", res = 120); if (draw) grid.draw(p) else print(p); dev.off() }
cm <- function(v) format(v, big.mark = ",", trim = TRUE)
fl <- fromJSON("results/flow.json")$flow; res <- read.csv("results/out_results.csv")
box <- function(db, x, n0) { f <- fl[[db]]
  data.frame(x = x, y = c(11, 9, 7, 5, 2.4), h = c(1.2, 1.0, 1.0, 1.4, 1.6),
    label = c(paste0(db, ": first ICU stay, age \u226518 years,\nacute pancreatitis (n = ", cm(n0), ")"),
      paste0("Excluded: no SpO2/FiO2 within 24 h (n = ", cm(n0 - f$sf_day1), ");\ndied within 72 h (n = ", cm(f$early), ")"),
      paste0("Alive at 72-h landmark (n = ", cm(f$landmark), ")\nRespiratory failure within 24 h: ", cm(f$rf)),
      paste0("Excluded, not classifiable (n = ", cm(f$unc_n), "): RF within 24 h but\nno SpO2/FiO2 at 48-72 h (ICU discharge\nbefore 72 h, n = ", cm(f$unc_icuout), ")"),
      paste0("Analytic cohort (n = ", cm(f$g0 + f$g1 + f$g2), ")\nNo RF: ", cm(f$g0), "\nTransient RF: ", cm(f$g1), "\nPersistent RF: ", cm(f$g2)))) }
b <- rbind(box("MIMIC-IV", 3, fl[["MIMIC-IV"]]$n_eligible), box("eICU-CRD", 9.2, fl[["eICU-CRD"]]$n_eligible)); b$w <- 5.4
f1 <- ggplot() + geom_rect(data = b, aes(xmin = x - w/2, xmax = x + w/2, ymin = y - h/2, ymax = y + h/2), fill = "white", colour = "black", linewidth = 0.4) +
  geom_text(data = b, aes(x, y, label = label), size = 2.9, lineheight = 0.95) +
  annotate("segment", x = rep(c(3, 9.2), each = 4), xend = rep(c(3, 9.2), each = 4), y = rep(c(10.4, 8.5, 6.5, 4.3), 2), yend = rep(c(9.5, 7.5, 5.7, 3.2), 2), arrow = arrow(length = unit(0.13, "cm")), linewidth = 0.4) +
  annotate("text", x = 6.1, y = 0.7, size = 2.5, label = "RF, respiratory failure (worst SpO2/FiO2 \u2264315 within 24 h); persistent/transient RF, SpO2/FiO2 \u2264315 or >315 observed at 48-72 h") +
  xlim(0, 12.2) + ylim(0.3, 11.8) + theme_void() + theme(plot.background = element_rect(fill = "white", colour = NA))
save_fig(f1, "Figure1", 8.5, 7)
lab <- c("No RF", "Transient RF", "Persistent RF"); pal <- c("#00468B", "#42B540", "#AD002A")
km <- read.csv("results/out_km_MIMIC-IV.csv"); km$grp <- factor(km$g, 0:2, lab)
fit <- survfit(Surv(tlm, y) ~ grp, data = km)
k <- ggsurvplot(fit, data = km, palette = pal, legend.labs = lab, legend.title = "", title = "A  MIMIC-IV: 28-day mortality", xlim = c(0, 25), break.time.by = 5, ylim = c(0.6, 1),
  censor = FALSE, size = 0.7, pval = TRUE, pval.coord = c(1, 0.65), risk.table = TRUE, risk.table.height = 0.28, fontsize = 3.2,
  xlab = "Days after the 72-h landmark", ylab = "Survival probability", ggtheme = theme_classic(base_size = 11), tables.theme = theme_cleantable())
ke <- read.csv("results/out_km_eICU-CRD.csv"); ag <- aggregate(y ~ g, ke, function(z) c(sum(z), length(z)))
wil <- function(x, n) { p <- x / n; z <- 1.96; c((p + z^2/(2*n) - z*sqrt(p*(1-p)/n + z^2/(4*n^2))) / (1 + z^2/n), (p + z^2/(2*n) + z*sqrt(p*(1-p)/n + z^2/(4*n^2))) / (1 + z^2/n)) }
eb <- data.frame(grp = factor(lab, lab), x = ag$y[, 1], n = ag$y[, 2]); eb$p <- eb$x / eb$n; ci <- t(mapply(wil, eb$x, eb$n)); eb$lo <- ci[, 1]; eb$hi <- ci[, 2]
pb <- ggplot(eb, aes(grp, 100 * p, fill = grp)) + geom_col(width = 0.6, colour = "black", linewidth = 0.3) + geom_errorbar(aes(ymin = 100 * lo, ymax = 100 * hi), width = 0.15) +
  geom_text(aes(label = sprintf("%d/%d", x, n), y = 100 * hi + 2), size = 3.3) + scale_fill_manual(values = pal, guide = "none") +
  labs(x = NULL, y = "In-hospital mortality, % (95% CI)", title = "B  eICU-CRD: in-hospital mortality") + theme_classic(base_size = 11)
f2 <- arrangeGrob(arrangeGrob(k$plot, k$table, ncol = 1, heights = c(3, 1.1)), pb, ncol = 2, widths = c(1.25, 1))
save_fig(f2, "Figure2", 11.5, 5.6, draw = TRUE)
mk <- function(db, ttl) { pp <- read.csv(paste0("results/out_pairs_", db, ".csv")); a <- read.csv(paste0("results/out_agree_", db, ".csv")); rho <- a[a$stat == "rho", ]
  ggplot(pp, aes(pf, sf)) + geom_point(alpha = 0.18, size = 0.8, colour = "#00468B") + geom_vline(xintercept = 300, linetype = "dashed", colour = "grey40") + geom_hline(yintercept = 315, linetype = "dotted", colour = "grey40") +
    annotate("text", x = 690, y = 60, hjust = 1, size = 3, label = sprintf("Spearman rho %.2f (95%% CI %.2f-%.2f)\n%s pairs, %s patients", rho$est, rho$lo, rho$hi, cm(rho$npairs), cm(rho$npat))) +
    labs(x = "PaO2/FiO2", y = "SpO2/FiO2", title = ttl) + theme_classic(base_size = 11) }
rocp <- function(db, ttl) { r <- read.csv(paste0("results/out_roc_", db, ".csv")); a <- read.csv(paste0("results/out_agree_", db, ".csv")); u <- a[a$stat == "auc", ]
  ggplot(r[order(r$fpr, r$tpr), ], aes(fpr, tpr)) + geom_step(colour = "#AD002A", linewidth = 0.8) + geom_abline(linetype = "dashed", colour = "grey50") +
    annotate("text", x = 0.97, y = 0.08, hjust = 1, size = 3, label = sprintf("AUC %.3f (95%% CI %.3f-%.3f)", u$est, u$lo, u$hi)) + labs(x = "1 - Specificity", y = "Sensitivity", title = ttl) + coord_equal() + theme_classic(base_size = 11) }
save_fig(arrangeGrob(mk("MIMIC-IV", "A  MIMIC-IV"), rocp("MIMIC-IV", "B  MIMIC-IV"), mk("eICU-CRD", "C  eICU-CRD"), rocp("eICU-CRD", "D  eICU-CRD"), ncol = 2), "Figure3", 10, 9, draw = TRUE)
nmz <- c(old = "Age", male = "Sex", vaso = "Vasopressor within 24 h", cr2 = "Creatinine \u22652 mg/dL", copd = "COPD")
lv <- list(old = c("<65 years", "\u226565 years"), male = c("Female", "Male"), vaso = c("No", "Yes"), cr2 = c("No", "Yes"), copd = c("No", "Yes"))
fp_db <- function(db) { sg <- res[res$sec == "subgroup" & res$db == db, ]; rows <- list()
  for (v in names(nmz)) { s <- sg[sub(" .*", "", sg$lab) == v, ]
    rows[[length(rows) + 1]] <- data.frame(Subgroup = nmz[[v]], N = "", Events = "", est = NA, lo = NA, hi = NA, OR = "", Pint = sprintf("%.3f", as.numeric(s$extra[1])))
    for (k in 1:2) { ok <- is.finite(s$hi[k]) && s$hi[k] < 1000
      rows[[length(rows) + 1]] <- data.frame(Subgroup = paste0("   ", lv[[v]][k]), N = cm(s$n[k]), Events = as.character(s$ev[k]), est = ifelse(ok, s$est[k], NA), lo = ifelse(ok, s$lo[k], NA), hi = ifelse(ok, s$hi[k], NA),
        OR = ifelse(ok, sprintf("%.2f (%.2f-%.2f)", s$est[k], s$lo[k], s$hi[k]), "Not estimable"), Pint = "") } }
  fd <- do.call(rbind, rows); fd$` ` <- paste(rep(" ", 20), collapse = " "); names(fd)[names(fd) == "OR"] <- "OR (95% CI)"; names(fd)[names(fd) == "Pint"] <- "P for interaction"
  fp <- forest(fd[, c("Subgroup", "N", "Events", " ", "OR (95% CI)", "P for interaction")], est = fd$est, lower = fd$lo, upper = fd$hi, ci_column = 4, ref_line = 1, x_trans = "log",
    xlim = c(0.1, 200), ticks_at = c(0.1, 1, 10, 100), xlab = "Adjusted OR (95% CI)", title = db, theme = forest_theme(base_size = 9, ci_pch = 15, ci_col = "black", refline_gp = gpar(lty = "dashed")))
  edit_plot(fp, row = which(fd$N == ""), gp = gpar(fontface = "bold"), part = "body", col = 1) }
save_fig(arrangeGrob(fp_db("MIMIC-IV"), fp_db("eICU-CRD"), ncol = 1), "Figure4", 8.5, 9.5, draw = TRUE)
cat("done\n")
