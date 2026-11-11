# =============================================================================
# 03_metric_power.R - can the usual group-comparison metrics detect a real difference?
# (a) Coverage curves: share of each group's cold-attributable deaths captured when the coldest
#     x% of days are flagged (regional Tmin percentile thresholds), x = 0.5 ... 60%.
# (b) Power of the metrics.  Parametric simulation at the second stage: group curves are drawn
#     from N(true curve, the group's REAL pooled covariance).  Truth: group B = all-cause curve;
#     group A = all-cause curve scaled so that its cold RR at P1 (vs the common MMT) is 1.00, 1.05,
#     1.10 or 1.20 times B's.  Metrics and what "correct" means (A is the more cold-vulnerable):
#       RRR at P1 (common MMT)        detection: two-sided p < 0.05 with RRR > 1
#       cold attributable fraction    A > B
#       own MMT (Cao et al. 2026)     not ordered by construction; spread reported
#       days colder than own MMT      A > B   (Cao exposure logic)
#       days with RR >= 1.05          A > B
#       "significant" cold days       A > B   (lower 95% bound > 1)
#       risk-weighted cold days       A > B   (sum of 1 - 1/RR)
#     Validated against the full null simulation (01, MMT spread) when its output exists.
# (c) Shrinkage of group curves toward each other (multivariate random-effects meta-analysis
#     across the 7 margins, BLUPs): how much of the spread in group MMTs survives?
# Input : models/pooled_Tmin.rds; tables/07_methods_supplement/null_by_size.csv (optional)
# Output: tables/07_methods_supplement/coverage_curves.csv, metric_power.csv,
#         metric_power_validation.csv, shrinkage_group_mmt.csv; figures/07_methods_supplement/*
# =============================================================================
source("2. Scripts/00_setup/config.R"); source("2. Scripts/00_setup/dlnm_tools.R")
STAGE <- "07_methods_supplement"; set.seed(SEED); NS <- 2000
panel <- load_panel(); setorder(panel, Region, Date); PL <- read_model("pooled_Tmin")
P <- setNames(lapply(GROUPS, function(g) PL[[paste(g, "cum", sep = "|")]]), GROUPS); pa <- P$All; mmt_all <- pa$mmt; r0 <- pidx(mmt_all)
w <- c(diff(PCT), 0)

# ---- (a) coverage curves ------------------------------------------------------------------------------
rb <- function(x, at) onebasis(at, fun = "ns", knots = quantile(x, c(.1, .75, .9), na.rm = TRUE), Boundary.knots = range(x, na.rm = TRUE))
cc <- rbindlist(lapply(GROUPS, function(g) { p <- P[[g]]
  d <- rbindlist(lapply(names(p$blup), function(r) { dd <- panel[Region == r]; x <- dd$Tmin_f; cen <- quantile(x, mmt_all)
    l <- drop(sweep(rb(x, x), 2, rb(x, cen), "-") %*% p$blup[[r]]$blup); a <- (1 - exp(-l)) * frollmean(dd[[g]], MAX_LAG + 1, align = "left", na.rm = TRUE)
    data.table(pct = frank(x) / length(x), an = fifelse(x < cen, a, 0)) }))
  tot <- sum(d$an, na.rm = TRUE)
  data.table(group = g, flagged = c(.005, .01, .025, .05, .10, .15, .20, .30, .40, .50, .60), coverage = sapply(c(.005, .01, .025, .05, .10, .15, .20, .30, .40, .50, .60), function(q) sum(d$an[d$pct <= q], na.rm = TRUE) / tot)) }))
write_tab(cc, STAGE, "coverage_curves")

# ---- (b) power of the metrics --------------------------------------------------------------------------
metrics <- function(cf, p) { lin <- unname(drop(p$basis %*% cf)); l0 <- lin - lin[r0]                  # vs the common MMT
  own <- INRANGE[which.min(lin[INRANGE])]; lo <- seq_along(PCT) < r0
  se <- sqrt(rowSums((sweep(p$basis, 2, p$basis[r0, ], "-") %*% p$vcov) * sweep(p$basis, 2, p$basis[r0, ], "-")))
  c(lrr_p1 = l0[pidx(.01)], AF = sum(w[lo] * (1 - exp(-l0[lo]))) , own_mmt = PCT[own], days_below_own_mmt = sum(w[seq_along(PCT) < own]),
    days_rr105 = sum(w[lo & l0 >= log(1.05)]), sig_days = sum(w[lo & l0 - 1.96 * se > 0]), riskweighted = sum(w[lo] * pmax(0, 1 - exp(-l0[lo])))) }
pairs <- list(sex = c("Female", "Male"), race = c("Black", "White"), age_old_mid = c("Age75plus", "Age50_74"), age_old_young = c("Age75plus", "Under50"))
lp1 <- drop((pa$basis[pidx(.01), ] - pa$basis[r0, ]) %*% pa$coef)
pw <- rbindlist(lapply(names(pairs), function(nm) { pr <- pairs[[nm]]; Va <- P[[pr[1]]]$vcov; Vb <- P[[pr[2]]]$vcov
  rbindlist(lapply(c(1, 1.05, 1.10, 1.20), function(eff) { k <- log(eff) / lp1; ta <- pa$coef * (1 + k); tb <- pa$coef
    pA <- modifyList(P[[pr[1]]], list(vcov = Va)); pB <- modifyList(P[[pr[2]]], list(vcov = Vb)); pA$basis <- pa$basis; pB$basis <- pa$basis
    MA <- t(apply(mvrnorm(NS, ta, Va), 1, metrics, p = pA)); MB <- t(apply(mvrnorm(NS, tb, Vb), 1, metrics, p = pB))
    xa <- pa$basis[pidx(.01), ] - pa$basis[r0, ]; se_d <- sqrt(drop(t(xa) %*% (Va + Vb) %*% xa)); z <- (MA[, "lrr_p1"] - MB[, "lrr_p1"]) / se_d
    data.table(pair = nm, true_RRR_P1 = eff,
      RRR_test_detects = mean(abs(z) > 1.96 & z > 0), RRR_test_wrong_sign = mean(abs(z) > 1.96 & z < 0),
      AF_correct = mean(MA[, "AF"] > MB[, "AF"]), cao_days_correct = mean(MA[, "days_below_own_mmt"] > MB[, "days_below_own_mmt"]),
      rr105_days_correct = mean(MA[, "days_rr105"] > MB[, "days_rr105"]), sig_days_correct = mean(MA[, "sig_days"] > MB[, "sig_days"]),
      riskweighted_correct = mean(MA[, "riskweighted"] > MB[, "riskweighted"]),
      own_mmt_gap_q05 = 100 * quantile(MA[, "own_mmt"] - MB[, "own_mmt"], .05), own_mmt_gap_q95 = 100 * quantile(MA[, "own_mmt"] - MB[, "own_mmt"], .95)) })) }))
write_tab(pw, STAGE, "metric_power")
# validation of the parametric shortcut against the full null simulation (MMT spread, identical curves)
nf <- out_path("tables", STAGE, "null_by_size.csv")
if (file.exists(nf)) { nul <- fread(nf)
map <- c(sB = "Black", sW = "White", sF = "Female", sM = "Male", sU50 = "Under50", s5074 = "Age50_74", s75 = "Age75plus")
val <- rbindlist(lapply(names(map), function(s) { g <- map[[s]]; pg <- P[[g]]; pg$basis <- pa$basis
  m <- apply(mvrnorm(NS, pa$coef, pg$vcov), 1, function(cf) PCT[INRANGE[which.min(drop(pa$basis %*% cf)[INRANGE])]])
  data.table(group = g, parametric_q05 = 100 * quantile(m, .05), parametric_q95 = 100 * quantile(m, .95),
             full_sim_q05 = nul[group == s, mmt_q05], full_sim_q95 = nul[group == s, mmt_q95]) }))
write_tab(val, STAGE, "metric_power_validation") }

# ---- (c) shrinkage of group curves ----------------------------------------------------------------------
gs <- setdiff(GROUPS, "All"); co <- t(sapply(gs, function(g) P[[g]]$coef)); S <- lapply(gs, function(g) P[[g]]$vcov)
mm <- mixmeta(co ~ 1, S = S, method = "reml"); bl <- blup(mm, vcov = FALSE)
cao <- function(cf) { lin <- unname(drop(pa$basis %*% cf)); own <- INRANGE[which.min(lin[INRANGE])]; c(mmt = 100 * PCT[own], hot_days = 100 * sum(w[seq_along(PCT) > own])) }
sh <- rbindlist(lapply(seq_along(gs), function(i) { a <- cao(co[i, ]); b <- cao(bl[i, ])
  data.table(group = gs[i], MMT_raw = a[["mmt"]], MMT_shrunk = b[["mmt"]], hot_days_raw = a[["hot_days"]], hot_days_shrunk = b[["hot_days"]]) }))
sh <- rbind(sh, data.table(group = "All", MMT_raw = 100 * mmt_all, MMT_shrunk = 100 * mmt_all, hot_days_raw = 100 * (1 - mmt_all), hot_days_shrunk = 100 * (1 - mmt_all)))
write_tab(sh, STAGE, "shrinkage_group_mmt")
cat("between-group heterogeneity I2 =", round(summary(mm)$i2stat[1], 1), "%; Cochran Q p =", signif(summary(mm)$qstat$pvalue[1], 2), "\n")

# ---- figures ---------------------------------------------------------------------------------------------------
fig <- function(g, n, w = 10, h = 4.5) save_fig(g, STAGE, n, w, h)
fig(ggplot(cc, aes(100 * flagged, 100 * coverage, colour = group)) + geom_abline(slope = 1, intercept = 0, colour = COL[["grey"]], linetype = 2) + geom_line(linewidth = .8) +
      geom_point(size = 1.4) + scale_colour_manual(values = GROUP_COL, labels = GROUP_LABEL, name = NULL) + scale_x_continuous(breaks = c(0, 2.5, 10, 20, 30, 40, 50, 60)) +
      labs(title = "How many cold days must be flagged to capture cold-attributable deaths?",
           subtitle = "Coldest x% of region-days (regional Tmin percentile) vs share of the group's cold-attributable deaths occurring on them. Dashed = proportional.",
           x = "% of days flagged (coldest first)", y = "% of cold-attributable deaths covered"), "coverage_curves", 9, 4.6)
pl <- melt(pw, id.vars = c("pair", "true_RRR_P1"), measure.vars = c("RRR_test_detects", "AF_correct", "riskweighted_correct", "rr105_days_correct", "cao_days_correct", "sig_days_correct"))
ML <- c(RRR_test_detects = "RRR test detects (p < 0.05)", AF_correct = "attributable fraction ranks correctly", riskweighted_correct = "risk-weighted days rank correctly",
        rr105_days_correct = "days RR ≥ 1.05 rank correctly", cao_days_correct = "days below own MMT rank correctly (Cao)", sig_days_correct = "'significant' days rank correctly")
PLB <- c(sex = "Female vs Male", race = "Black vs White", age_old_mid = "75+ vs 50–74", age_old_young = "75+ vs Under 50")
fig(ggplot(pl, aes(true_RRR_P1, 100 * value, colour = variable)) + geom_hline(yintercept = c(5, 50), colour = COL[["grey"]], linetype = 3) + geom_line(linewidth = .8) + geom_point(size = 1.6) +
      facet_wrap(~factor(PLB[pair], PLB), nrow = 1) + scale_x_continuous(breaks = c(1, 1.05, 1.1, 1.2)) +
      scale_colour_manual(values = c(COL[["blue"]], COL[["orange"]], COL[["aqua"]], "#eda100", "#e34948", "#6250d6"), labels = ML, name = NULL) +
      guides(colour = guide_legend(nrow = 2)) +
      labs(title = "Can each metric tell which group is more cold-vulnerable?", subtitle = "2,000 draws per scenario with the real group covariances; group A's true cold RR at P1 is 1.00–1.20 times group B's. 50% = coin toss.",
           x = "True ratio of cold RR at P1 (A / B)", y = "% of draws"), "metric_power", 13, 4.8)
sl <- melt(sh[group != "All"], id.vars = "group", measure.vars = c("MMT_raw", "MMT_shrunk"))
fig(ggplot(sl, aes(value, factor(GROUP_LABEL[group], rev(GROUP_LABEL)), colour = variable)) + geom_vline(xintercept = 100 * mmt_all, colour = COL[["grey"]], linetype = 2) +
      geom_point(size = 3) + scale_colour_manual(values = c(MMT_raw = COL[["orange"]], MMT_shrunk = COL[["blue"]]), labels = c("group's own curve", "after empirical-Bayes shrinkage"), name = NULL) +
      labs(title = "Group MMTs before and after shrinkage toward the common curve", subtitle = "Multivariate random-effects model across the 7 margins (BLUPs). Dashed: all-cause MMT.", x = "MMT (percentile of Tmin)", y = NULL), "shrinkage_mmt", 8.5, 3.8)
print(pw[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)]); print(sh)
