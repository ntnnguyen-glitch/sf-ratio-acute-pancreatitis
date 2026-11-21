## Statistical analysis (R 4.6.1).
## Input: data/analytic.csv and data/pairs.csv from python/03_build_variables.py.
## Output: results/*.csv and results/flow.json.
## Each database is analysed separately; persistence of respiratory failure requires an observed SpO2/FiO2 at 48-72 h.
suppressPackageStartupMessages({library(survival);library(pROC);library(WeightIt);library(MatchIt);library(cobalt)
  library(mice);library(EValue);library(sandwich);library(lmtest);library(jsonlite)})
set.seed(2026); options(warn = -1)
dir.create("results", showWarnings = FALSE)
d0 <- read.csv("data/analytic.csv"); pr <- read.csv("data/pairs.csv")
n_elig <- table(d0$db)
d0 <- d0[!is.na(d0$sf1), ]
d0$early <- as.integer(!is.na(d0$tdeath) & d0$tdeath <= 3)
cont <- c("map_min","gcs","lactate","creat","bun","wbc","plt","bili","hr","rr")
## exposure groups: 0 no RF; 1 transient (observed S/F >315 at 48-72 h); 2 persistent (observed S/F <=315 at 48-72 h); NA unclassifiable
grp <- function(s1, s3, cut = 315) ifelse(s1 > cut, 0, ifelse(is.na(s3), NA, ifelse(s3 <= cut, 2, 1)))
d0$g <- grp(d0$sf1, d0$sf3)
d0$g_orig <- ifelse(d0$sf1 > 315, 0, ifelse(!is.na(d0$sf3) & d0$sf3 <= 315, 2, 1))
d0$g235 <- grp(d0$sf1, d0$sf3, 235)
s1_97 <- ifelse(is.na(d0$sf1_97), 999, d0$sf1_97); s3_97 <- ifelse(is.na(d0$sf3_97) & !is.na(d0$sf3), 999, d0$sf3_97); d0$g97 <- grp(s1_97, s3_97)
d0$g_pf <- ifelse(is.na(d0$pf1), NA, ifelse(d0$pf1 > 300, 0, ifelse(is.na(d0$pf3), NA, ifelse(d0$pf3 <= 300, 2, 1))))
d0$mm_sf <- as.numeric(as.character(cut(d0$sf1, c(-1, 148, 232, 316, 400, 9999), labels = c(4, 3, 2, 1, 0))))
d0$y <- ifelse(d0$db == "MIMIC-IV", as.integer(!is.na(d0$tdeath) & d0$tdeath <= 28), d0$dhosp)
d0$d90 <- as.integer(!is.na(d0$tdeath) & d0$tdeath <= 90)
d0$tlm <- pmin(ifelse(is.na(d0$tdeath), 28, d0$tdeath), 28) - 3
d0$icu_after <- pmax(0, d0$icu_h - 72) / 24
d0$ifd <- ifelse(d0$y == 1 & (d0$db == "eICU-CRD" | (!is.na(d0$tdeath) & d0$tdeath <= 28)), 0, pmax(0, 25 - pmin(25, d0$icu_after)))
d0$rrt_new_rs <- ifelse(d0$rrt_pre72 == 1, NA, d0$rrt_new); d0$imv_new_rs <- ifelse(d0$imv_0_72 == 1, NA, d0$imv_new72)
d0$old <- as.integer(d0$age >= 65); d0$cr2 <- as.integer(!is.na(d0$creat) & d0$creat >= 2)
imp <- function(df) { for (v in cont) { df[[paste0(v, "_m")]] <- as.integer(is.na(df[[v]])); df[[v]][is.na(df[[v]])] <- median(df[[v]], na.rm = TRUE) }; df }
M2 <- "age + male + copd + chf + ckd + liver + dm + cancer"
M3 <- paste(M2, "+ alc + bil + vaso + imv24 + map_min + gcs + log(lactate) + log(creat) + log(bun) + log(wbc) + log(plt) + log(bili)")
bins <- c("copd","chf","ckd","liver","dm","cancer","alc","bil","vaso","imv24","male","old","cr2")
fx <- function(df, f) { for (v in bins) if (v %in% names(df) && length(unique(df[[v]])) < 2) f <- gsub(paste0("\\+ ", v, "\\b"), "", f); f }
res <- data.frame(); add <- function(sec, lab, db, n, ev, est, lo, hi, p, kind = "OR", extra = NA) res <<- rbind(res, data.frame(sec, lab, db, n, ev, est, lo, hi, p, kind, extra))
rob <- function(m, df, term) { v <- if (df$db[1] == "eICU-CRD") vcovCL(m, cluster = df$hosp) else vcovHC(m, type = "HC0"); s <- coeftest(m, vcov. = v)[term, ]; s }
lg <- function(df, f, term = "e", cluster = TRUE, family = binomial) { m <- glm(as.formula(fx(df, f)), data = df, family = family)
  s <- if (cluster && df$db[1] == "eICU-CRD") rob(m, df, term) else summary(m)$coef[term, ]
  c(exp(s[1]), exp(s[1] - 1.96 * s[2]), exp(s[1] + 1.96 * s[2]), s[4]) }
pe <- function(df, g = "g", a = 0, b = 2) { x <- df[!is.na(df[[g]]) & df[[g]] %in% c(a, b), ]; x$e <- as.integer(x[[g]] == b); x }
flow <- list(); DBn <- c("MIMIC-IV", "eICU-CRD"); A <- list()
for (nm in DBn) {
  raw <- d0[d0$db == nm, ]; lm_ <- raw[raw$early == 0, ]; rf <- lm_[lm_$sf1 <= 315, ]
  flow[[nm]] <- list(n_eligible = as.integer(n_elig[[nm]]), sf_day1 = nrow(raw), early = sum(raw$early), landmark = nrow(lm_), rf = nrow(rf), rf_obs = sum(!is.na(rf$sf3)),
    unc_icuout = sum(is.na(rf$sf3) & rf$icu_h < 72), unc_inicu = sum(is.na(rf$sf3) & rf$icu_h >= 72),
    unc_died = sum(is.na(rf$sf3) & rf$y == 1), g0 = sum(lm_$g == 0, na.rm = TRUE), g1 = sum(lm_$g == 1, na.rm = TRUE), g2 = sum(lm_$g == 2, na.rm = TRUE),
    pf_day1 = sum(!is.na(lm_$pf1)), rf_pf_obs = sum(!is.na(rf$pf3)), unc_y = sum(rf$y[is.na(rf$sf3)]), unc_n = sum(is.na(rf$sf3)))
  L <- imp(lm_[!is.na(lm_$g), ]); A[[nm]] <- L
  for (m in list(c("Crude", "y ~ e"), c("Model 2", paste("y ~ e +", M2)), c("Model 3", paste("y ~ e +", M3)))) {
    x <- pe(L); o <- lg(x, m[2]); add("primary", paste("Persistent vs no RF |", m[1]), nm, nrow(x), sum(x$y), o[1], o[2], o[3], o[4]) }
  x <- pe(L, a = 1, b = 2); o <- lg(x, paste("y ~ e +", M3)); add("primary", "Persistent vs transient RF | Model 3", nm, nrow(x), sum(x$y), o[1], o[2], o[3], o[4])
  L$tr <- L$g; o <- lg(L, paste("y ~ tr +", M3), "tr"); add("primary", "Per category (trend) | Model 3", nm, nrow(L), sum(L$y), o[1], o[2], o[3], o[4])
  L$sfd <- -L$sf1 / 50; o <- lg(L, paste("y ~ sfd +", M3), "sfd"); add("primary", "Worst S/F per 50-unit decrease | Model 3", nm, nrow(L), sum(L$y), o[1], o[2], o[3], o[4])
  o <- lg(L, paste("y ~ mm_sf +", M3), "mm_sf"); add("primary", "Marshall respiratory score (S/F) per point | Model 3", nm, nrow(L), sum(L$y), o[1], o[2], o[3], o[4])
  if (nm == "MIMIC-IV") { x <- pe(L); cx <- coxph(as.formula(fx(x, paste("Surv(tlm, y) ~ e +", M3))), data = x); s <- summary(cx)
    add("primary", "Persistent vs no RF | Cox Model 3", nm, nrow(x), sum(x$y), s$conf.int["e", 1], s$conf.int["e", 3], s$conf.int["e", 4], s$coef["e", 5], "HR") }
  Lp <- imp(lm_[!is.na(lm_$g_pf), ]); x <- pe(Lp, "g_pf"); o <- tryCatch(lg(x, paste("y ~ e +", M3)), error = function(e) rep(NA, 4))
  add("pf", "Persistent vs no RF (P/F definition) | Model 3", nm, nrow(x), sum(x$y), o[1], o[2], o[3], o[4])
  ## secondary
  sec <- list(c("rrt_new_rs", "New RRT 72 h-7 d"), c("imv_new_rs", "New IMV 72 h-28 d"))
  if (nm == "MIMIC-IV") sec <- c(list(c("dhosp", "In-hospital mortality"), c("d90", "90-day mortality")), sec)
  for (s_ in sec) { x <- pe(L); x <- x[!is.na(x[[s_[1]]]), ]; o <- tryCatch(lg(x, paste(s_[1], "~ e +", M3)), error = function(e) rep(NA, 4))
    add("secondary", s_[2], nm, nrow(x), sum(x[[s_[1]]]), o[1], o[2], o[3], o[4]) }
  x <- pe(L); m <- lm(as.formula(fx(x, paste("ifd ~ e +", M3))), data = x); s <- rob(m, x, "e")
  add("secondary", "ICU-free days to day 28 (mean difference)", nm, nrow(x), NA, s[1], s[1] - 1.96 * s[2], s[1] + 1.96 * s[2], s[4], "MD")
  ## confounding
  x <- pe(L); for (v in c("lactate","creat","bun","wbc","plt","bili")) x[[paste0("l", v)]] <- log(x[[v]])
  psf <- as.formula(fx(x, paste("e ~", M3)))
  W1 <- weightit(psf, data = x, method = "glm", estimand = "ATO"); m <- glm(y ~ e, data = x, weights = W1$weights, family = quasibinomial); s <- rob(m, x, "e")
  ess <- function(w, t) c(sum(w[t == 1])^2 / sum(w[t == 1]^2), sum(w[t == 0])^2 / sum(w[t == 0]^2))
  es <- ess(W1$weights, x$e); add("confounding", "Overlap weighting", nm, nrow(x), sum(x$y), exp(s[1]), exp(s[1] - 1.96 * s[2]), exp(s[1] + 1.96 * s[2]), s[4], extra = paste(round(es), collapse = "/"))
  W2 <- trim(weightit(psf, data = x, method = "glm", estimand = "ATE", stabilize = TRUE), at = 0.99)
  m <- glm(as.formula(fx(x, paste("y ~ e +", M3))), data = x, weights = W2$weights, family = quasibinomial); s <- rob(m, x, "e")
  add("confounding", "IPTW with outcome-model adjustment (doubly robust)", nm, nrow(x), sum(x$y), exp(s[1]), exp(s[1] - 1.96 * s[2]), exp(s[1] + 1.96 * s[2]), s[4])
  mt <- matchit(psf, data = x, method = "nearest", distance = "glm", caliper = 0.2, std.caliper = TRUE); md <- match.data(mt)
  m <- glm(y ~ e, data = md, weights = weights, family = quasibinomial); s <- coeftest(m, vcov. = vcovCL(m, cluster = ~subclass))["e", ]
  add("confounding", "PSM 1:1 (caliper 0.2 SD)", nm, nrow(md), sum(md$y), exp(s[1]), exp(s[1] - 1.96 * s[2]), exp(s[1] + 1.96 * s[2]), s[4], extra = sum(md$e))
  cv <- c("age","male","copd","chf","ckd","liver","dm","cancer","alc","bil","vaso","imv24","map_min","gcs","llactate","lcreat","lbun","lwbc","lplt","lbili")
  bt <- function(dd, w = NULL) { b <- bal.tab(dd[, cv], treat = dd$e, weights = w, s.d.denom = "pooled", binary = "std")$Balance; if (is.null(w)) b$Diff.Un else b$Diff.Adj }
  smd <- data.frame(var = cv, before = bt(x), ow = bt(x, W1$weights), iptw = bt(x, W2$weights), psm = bt(md)); smd$db <- nm
  write.csv(smd, paste0("results/out_smd_", nm, ".csv"), row.names = FALSE)
  p3 <- tail(res[res$sec == "primary" & res$db == nm & res$lab == "Persistent vs no RF | Model 3", ], 1)
  ev <- evalues.OR(p3$est, p3$lo, p3$hi, rare = mean(x$y) < 0.15); add("evalue", "E-value point/lower", nm, NA, NA, ev[2, 1], ev[2, 2], NA, NA, "E")
  ## subgroups
  for (v in c("old","male","vaso","cr2","copd")) { fa <- fx(x, paste("y ~ e +", M3)); if (!grepl(paste0("\\b", v, "\\b"), fa)) fa <- paste(fa, "+", v)
    fi <- sub("y ~ e", paste0("y ~ e * ", v), fa); pint <- tryCatch(anova(glm(as.formula(fa), data = x, family = binomial), glm(as.formula(fi), data = x, family = binomial), test = "LRT")[2, "Pr(>Chi)"], error = function(e) NA)
    for (k in 0:1) { s_ <- x[x[[v]] == k, ]; o <- tryCatch(lg(s_, paste("y ~ e +", M3)), error = function(e) rep(NA, 4)); add("subgroup", paste(v, k), nm, nrow(s_), sum(s_$y), o[1], o[2], o[3], o[4], extra = pint) } }
  ## sensitivity
  sens <- function(lab, df, g = "g", f = paste("y ~ e +", M3), cl = TRUE) { x <- pe(df, g); o <- tryCatch(lg(x, f, cluster = cl), error = function(e) rep(NA, 4)); add("sensitivity", lab, nm, nrow(x), sum(x$y), o[1], o[2], o[3], o[4]) }
  Lall <- imp(lm_)
  sens("Original rule (ICU discharge before 72 h counted as resolved)", Lall, "g_orig")
  sens("All groups restricted to patients with S/F observed at 48-72 h", L[!is.na(L$sf3), ])
  sens("S/F <=235 for respiratory failure", imp(lm_[!is.na(lm_$g235), ]), "g235")
  sens("S/F computed only when SpO2 <=97% (both windows)", imp(lm_[!is.na(lm_$g97), ]), "g97")
  sens("Acute pancreatitis as principal diagnosis", L[L$sq1 == 1, ]); sens("Excluding COPD", L[L$copd == 0, ]); sens("Excluding IMV at ICU admission", L[L$imv0 == 0, ])
  Lx <- L; av <- c("paco2","alb","ca","glu")[sapply(c("paco2","alb","ca","glu"), function(v) any(!is.na(Lx[[v]])))]
  for (v in av) Lx[[v]][is.na(Lx[[v]])] <- median(Lx[[v]], na.rm = TRUE)
  lab_ <- c(paco2 = "PaCO2", alb = "albumin", ca = "calcium", glu = "glucose")[av]
  sens(paste("Further adjusted for", paste(lab_, collapse = ", ")), Lx, f = paste("y ~ e +", M3, "+", paste(av, collapse = " + ")))
  sens("Missing-indicator adjustment", L, f = paste("y ~ e +", M3, "+ lactate_m + bili_m + gcs_m"))
  if (nm == "eICU-CRD") sens("Model-based standard errors without hospital clustering", L, cl = FALSE)
  xr <- pe(lm_[!is.na(lm_$g), ]); im <- mice(xr[, c("y","e","age","male","copd","chf","ckd","liver","dm","cancer","alc","bil","vaso","imv24", cont)], m = 5, method = "pmm", maxit = 10, printFlag = FALSE, seed = 2026)
  fit <- with(im, glm(as.formula(fx(xr, paste("y ~ e +", M3))), family = binomial)); pl <- summary(pool(fit), conf.int = TRUE); s <- pl[pl$term == "e", ]
  add("sensitivity", "Multiple imputation (MICE, m = 5)", nm, nrow(xr), sum(xr$y), exp(s$estimate), exp(s$`2.5 %`), exp(s$`97.5 %`), s$p.value)
  ## agreement on simultaneous pairs, patient-level cluster bootstrap
  q <- pr[pr$db == nm & pr$pf >= 30 & pr$pf <= 700 & pr$sf >= 50 & pr$sf <= 500 & pr$id %in% lm_$id, ]
  st <- function(qq) { k <- function(a, b) { t <- table(factor(a, c(FALSE, TRUE)), factor(b, c(FALSE, TRUE))); n <- sum(t); po <- sum(diag(t)) / n; pe_ <- sum(rowSums(t) * colSums(t)) / n^2; (po - pe_) / (1 - pe_) }
    c(rho = cor(qq$sf, qq$pf, method = "spearman"), auc = as.numeric(auc(roc(as.integer(qq$pf <= 300), qq$sf, direction = ">", quiet = TRUE))),
      kappa = k(qq$sf <= 315, qq$pf <= 300), sens = mean(qq$sf[qq$pf <= 300] <= 315), spec = mean(qq$sf[qq$pf > 300] > 315)) }
  est <- st(q); ids <- unique(q$id); sp <- split(q, q$id)
  B <- t(replicate(1000, { bb <- sample(ids, replace = TRUE); st(do.call(rbind, sp[as.character(bb)])) }))
  ag <- data.frame(stat = names(est), est = est, lo = apply(B, 2, quantile, .025), hi = apply(B, 2, quantile, .975), db = nm, npairs = nrow(q), npat = length(ids))
  write.csv(ag, paste0("results/out_agree_", nm, ".csv"), row.names = FALSE); write.csv(q[, c("pf","sf")], paste0("results/out_pairs_", nm, ".csv"), row.names = FALSE)
  rr <- roc(as.integer(q$pf <= 300), q$sf, direction = ">", quiet = TRUE); write.csv(data.frame(fpr = 1 - rr$specificities, tpr = rr$sensitivities), paste0("results/out_roc_", nm, ".csv"), row.names = FALSE)
  ## KM / mortality by group
  write.csv(L[, c("tlm","y","g")], paste0("results/out_km_", nm, ".csv"), row.names = FALSE)
  if (nm == "MIMIC-IV") { sd <- survdiff(Surv(tlm, y) ~ g, data = L); add("km", "Log-rank", nm, nrow(L), sum(L$y), NA, NA, NA, 1 - pchisq(sd$chisq, 2), "test") }
  ## Table 1
  vars_c <- c("age","map_min","hr","rr","gcs","lactate","creat","bun","wbc","plt","bili","sf1","pf1","icu_h","hosp_los")
  vars_b <- c("male","copd","chf","ckd","liver","dm","cancer","alc","bil","imv24","vaso","rrt_pre72","y")
  g <- lm_[!is.na(lm_$g), ]; t1 <- data.frame()
  q1 <- function(z) sprintf("%.1f (%.1f, %.1f)", median(z, na.rm = TRUE), quantile(z, .25, na.rm = TRUE), quantile(z, .75, na.rm = TRUE))
  q2 <- function(z) sprintf("%d (%.1f)", sum(z, na.rm = TRUE), 100 * mean(z, na.rm = TRUE))
  for (v in vars_c) t1 <- rbind(t1, data.frame(var = v, g0 = q1(g[[v]][g$g == 0]), g1 = q1(g[[v]][g$g == 1]), g2 = q1(g[[v]][g$g == 2]), miss = sum(is.na(g[[v]])), p = tryCatch(kruskal.test(g[[v]] ~ g$g)$p.value, error = function(e) NA)))
  for (v in vars_b) t1 <- rbind(t1, data.frame(var = v, g0 = q2(g[[v]][g$g == 0]), g1 = q2(g[[v]][g$g == 1]), g2 = q2(g[[v]][g$g == 2]), miss = 0, p = suppressWarnings(chisq.test(table(g[[v]], g$g))$p.value)))
  t1 <- rbind(data.frame(var = "N", g0 = sum(g$g == 0), g1 = sum(g$g == 1), g2 = sum(g$g == 2), miss = 0, p = NA), t1); t1$db <- nm
  write.csv(t1, paste0("results/out_t1_", nm, ".csv"), row.names = FALSE)
}
write.csv(res, "results/out_results.csv", row.names = FALSE)
write_json(list(flow = flow, R = R.version.string), "results/flow.json", auto_unbox = TRUE, digits = 8)
print(res[, c("sec","lab","db","n","ev","est","lo","hi","p","extra")], digits = 3, row.names = FALSE)
str(flow)
