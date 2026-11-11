# =============================================================================
# 03_threshold_model.R - estimated temperature thresholds (hockey-stick model) vs the MMT.
#   model (00_setup/threshold_tools.R; Armstrong 2006; Baccini et al. 2008): conditional Poisson
#   with strata region x year x month x day of week; cold = (tau_c - Tmin)_+ over lags 0-21 (ns lag
#   basis), heat = (Tmin - tau_h)_+ over lags 0-2; tau = a regional Tmin percentile common to the
#   12 regions.  Thresholds by profile likelihood on a 1-percentile grid (cold P30-P99, heat
#   P50-P99); 95% CI = { q : (D(q) - D_min) / phi <= 3.84 }.  Slopes = cumulative % change in deaths
#   per 1 C beyond the threshold, SEs clustered by region x year-month.  Safe band = heat - cold threshold.
#   group differences: likelihood-ratio test of a common threshold (each group keeps its own strata
#   and slopes).
#   population at risk (Cao): person-days below / above each group's threshold vs the all-deaths
#   threshold, intervals from likelihood-weighted threshold draws.
#   known null: all deaths split at random (binomial) into two pseudo-groups with the SAME true
#   curve, shares 50/50 and 20/80, NREP replicates each; per replicate the A - B difference in the
#   MMT and significance threshold (two-stage DLNM), in the estimated thresholds and the LR test, in
#   the 10% risk-level temperature and the slope test at the all-deaths thresholds (used by 04).
# Output: tables/06_thresholds/tm_thresholds.csv, tm_profiles.csv, tm_pair_tests.csv,
#         tm_population_at_risk.csv, tm_partition_totals.csv, null_split_replicates.csv,
#         null_split_summary.csv; models/threshold_models.rds
# =============================================================================
source("2. Scripts/00_setup/config.R"); source("2. Scripts/00_setup/dlnm_tools.R"); source("2. Scripts/00_setup/threshold_tools.R")
STAGE <- "06_thresholds"; set.seed(SEED); NREP <- 50; PROBS <- c(0.5, 0.2); ND <- 2000
panel <- load_panel(); setorder(panel, Region, Date); regs <- unique(panel$Region)
rp <- fread(file.path(DERIVED, "regional_populations.csv"))[Year %in% YEARS]; EX <- exposure_tools(panel, rp)
S <- tm_setup(panel); QC <- S$QC; QH <- S$QH

# ---- (1) all deaths first (its thresholds are the common rule), then every margin ---------------------------------------
fit_group <- function(y, qc0 = NULL, qh0 = NULL) { p <- tm_grid_slopes(S, y, tm_profile(S, y))
  p$common <- if (is.null(qc0)) NULL else tm_fit(S, y, qc0, qh0, se = TRUE)[c("bc", "sc", "bh", "sh")]; p }
P_ALL <- fit_group(panel$All); QC0 <- P_ALL$qc; QH0 <- P_ALL$qh; P_ALL$common <- P_ALL[c("bc", "sc", "bh", "sh")]
PR <- c(list(All = P_ALL), mclapply(setNames(setdiff(GROUPS, "All"), setdiff(GROUPS, "All")), function(g) fit_group(panel[[g]], QC0, QH0), mc.cores = NCORES))
save_model(list(models = PR, common = c(qc = QC0, qh = QH0), grids = list(QC = QC, QH = QH)), "threshold_models"); cat("groups fitted\n")
pc_ <- function(b) 100 * (exp(b) - 1)
TG <- rbindlist(lapply(GROUPS, function(g) { p <- PR[[g]]; cc <- tm_ci(p$pc, QC, p$phi); ch <- tm_ci(p$ph, QH, p$phi)
  data.table(group = g, deaths = sum(panel[[g]], na.rm = TRUE),
             cold_thr_pct = 100 * p$qc, cold_lo = 100 * cc[1], cold_hi = 100 * cc[2], cold_thr_C = EX$pct_C(g, p$qc), cold_lo_C = EX$pct_C(g, cc[1]), cold_hi_C = EX$pct_C(g, cc[2]),
             cold_slope = pc_(p$bc), cold_slope_lo = pc_(p$bc - 1.96 * p$sc), cold_slope_hi = pc_(p$bc + 1.96 * p$sc), p_cold_slope = 2 * pnorm(-abs(p$bc / p$sc)),
             heat_thr_pct = 100 * p$qh, heat_lo = 100 * ch[1], heat_hi = 100 * ch[2], heat_thr_C = EX$pct_C(g, p$qh), heat_lo_C = EX$pct_C(g, ch[1]), heat_hi_C = EX$pct_C(g, ch[2]),
             heat_slope = pc_(p$bh), heat_slope_lo = pc_(p$bh - 1.96 * p$sh), heat_slope_hi = pc_(p$bh + 1.96 * p$sh), p_heat_slope = 2 * pnorm(-abs(p$bh / p$sh)),
             band_C = EX$pct_C(g, p$qh) - EX$pct_C(g, p$qc), phi = p$phi) }))
write_tab(TG, STAGE, "tm_thresholds")
write_tab(rbindlist(lapply(GROUPS, function(g) rbind(data.table(group = g, side = "cold", q = QC, dD = (PR[[g]]$pc - min(PR[[g]]$pc)) / PR[[g]]$phi),
                                                     data.table(group = g, side = "heat", q = QH, dD = (PR[[g]]$ph - min(PR[[g]]$ph)) / PR[[g]]$phi)))), STAGE, "tm_profiles")

# ---- (2) pair tests: common threshold (LR) -------------------------------------------------------------------------------------
TP <- rbindlist(lapply(PAIRS, function(pr) { a <- PR[[pr[1]]]; b <- PR[[pr[2]]]; phi <- mean(c(a$phi, b$phi))
  dc <- tm_draws(a$pc, QC, a$phi, ND) - tm_draws(b$pc, QC, b$phi, ND); dh <- tm_draws(a$ph, QH, a$phi, ND) - tm_draws(b$ph, QH, b$phi, ND)
  data.table(pair = paste(pr, collapse = "/"), cold_diff = 100 * (a$qc - b$qc), cold_diff_lo = 100 * quantile(dc, .025), cold_diff_hi = 100 * quantile(dc, .975), cold_LR_p = tm_lr_common(a$pc, b$pc, phi)[["p"]],
             heat_diff = 100 * (a$qh - b$qh), heat_diff_lo = 100 * quantile(dh, .025), heat_diff_hi = 100 * quantile(dh, .975), heat_LR_p = tm_lr_common(a$ph, b$ph, phi)[["p"]]) }))

# ---- (3) Cao population at risk with the estimated thresholds ---------------------------------------------------------------------------------
DRW <- lapply(PR, function(p) list(cold = tm_draws(p$pc, QC, p$phi, 1000), heat = tm_draws(p$ph, QH, p$phi, 1000)))
PAR <- rbindlist(lapply(c("cold", "heat"), function(side) { qk <- if (side == "cold") "qc" else "qh"
  rbindlist(lapply(setdiff(GROUPS, "All"), function(g) { own <- PR[[g]][[qk]]; com <- PR$All[[qk]]; po <- EX$pd_pct(g, own, side); pc <- EX$pd_pct(g, com, side)
    rel <- sapply(1:1000, function(s) EX$pd_pct(g, DRW[[g]][[side]][s], side) / EX$pd_pct(g, DRW$All[[side]][s], side) - 1)
    data.table(side = side, group = g, own_pct = 100 * own, common_pct = 100 * com, share_own = 100 * po / EX$PY[[g]], share_common = 100 * pc / EX$PY[[g]],
               delta_pct = 100 * (po / pc - 1), delta_lo = 100 * quantile(rel, .025), delta_hi = 100 * quantile(rel, .975)) })) }))
write_tab(PAR, STAGE, "tm_population_at_risk")
PTOT <- rbindlist(lapply(c("cold", "heat"), function(side) { qk <- if (side == "cold") "qc" else "qh"
  rbindlist(lapply(names(PARTITIONS), function(nm) { gs <- PARTITIONS[[nm]]; com <- sum(sapply(gs, EX$pd_pct, p = PR$All[[qk]], side = side)); own <- sum(sapply(gs, function(g) EX$pd_pct(g, PR[[g]][[qk]], side)))
    rel <- sapply(1:1000, function(s) sum(sapply(gs, function(g) EX$pd_pct(g, DRW[[g]][[side]][s], side))) / sum(sapply(gs, EX$pd_pct, p = DRW$All[[side]][s], side = side)) - 1)
    data.table(side = side, partition = nm, E_common_millions = com / 1e6, E_own_millions = own / 1e6, delta_pct = 100 * (own / com - 1), delta_lo = 100 * quantile(rel, .025), delta_hi = 100 * quantile(rel, .975)) })) }))
write_tab(PTOT, STAGE, "tm_partition_totals"); cat("population at risk done\n")

# ---- (4) known null: random split of all deaths -------------------------------------------------------------------------------------------------
dlnm_metrics <- function(d, col) { s1 <- Filter(Negate(is.null), lapply(regs, function(r) { z <- tryCatch(fit_region(d[Region == r], "Tmin_f", col), error = function(e) NULL)
  if (is.null(z)) NULL else c(list(region = r, group = col), z) }))
  thr_pair(pool(s1, col, "cum"), pool(s1, col, "acute")) }
t10 <- function(p, side) EX$temp_C("All", risk_level_temp(S, if (side == "cold") p$qc else p$qh, if (side == "cold") p$bc else p$bh, .10, side))
jobs <- CJ(prob = PROBS, rep = seq_len(NREP))
NUL <- rbindlist(mclapply(seq_len(nrow(jobs)), function(j) { set.seed(SEED + 4300 + j); sp <- split_null(panel$All, jobs$prob[j])
  d <- copy(panel)[, `:=`(gA = sp$a, gB = sp$b)]
  pa <- tm_profile(S, sp$a); pb <- tm_profile(S, sp$b); phi <- mean(c(pa$phi, pb$phi))
  ca <- tm_fit(S, sp$a, QC0, QH0, se = TRUE); cb <- tm_fit(S, sp$b, QC0, QH0, se = TRUE)
  ma <- dlnm_metrics(d, "gA"); mb <- dlnm_metrics(d, "gB")
  data.table(prob = jobs$prob[j], rep = jobs$rep[j],
             thr_cold_diff = 100 * (pa$qc - pb$qc), thr_heat_diff = 100 * (pa$qh - pb$qh), LR_p_cold = tm_lr_common(pa$pc, pb$pc, phi)[["p"]], LR_p_heat = tm_lr_common(pa$ph, pb$ph, phi)[["p"]],
             mmt_cold_diff = 100 * (ma[["cold_mmt"]] - mb[["cold_mmt"]]), mmt_heat_diff = 100 * (ma[["heat_mmt"]] - mb[["heat_mmt"]]),
             sig_cold_diff = 100 * (ma[["cold_sig"]] - mb[["cold_sig"]]), sig_heat_diff = 100 * (ma[["heat_sig"]] - mb[["heat_sig"]]),
             T10_cold_diff = t10(pa, "cold") - t10(pb, "cold"), T10_heat_diff = t10(pa, "heat") - t10(pb, "heat"),
             p_slope_cold = 2 * pnorm(-abs((ca$bc - cb$bc) / sqrt(ca$sc^2 + cb$sc^2))), p_slope_heat = 2 * pnorm(-abs((ca$bh - cb$bh) / sqrt(ca$sh^2 + cb$sh^2)))) }, mc.cores = NCORES))
write_tab(NUL, STAGE, "null_split_replicates")
MET <- c(mmt_cold_diff = "MMT, cold", sig_cold_diff = "significance threshold, cold", thr_cold_diff = "estimated threshold, cold",
         mmt_heat_diff = "MMT, heat", sig_heat_diff = "significance threshold, heat", thr_heat_diff = "estimated threshold, heat")
NS <- melt(NUL, id.vars = "prob", measure.vars = names(MET))[, .(metric = MET[[as.character(variable[1])]], q95_abs_diff = quantile(abs(value), .95, na.rm = TRUE),
  lo = quantile(value, .025, na.rm = TRUE), hi = quantile(value, .975, na.rm = TRUE), share_missing = mean(is.na(value))), by = .(prob, variable)]
FP <- NUL[, .(variable = c("thr_cold_diff", "thr_heat_diff"), LR_false_positive = c(mean(LR_p_cold < .05), mean(LR_p_heat < .05))), by = prob]
NS <- merge(NS, FP, by = c("prob", "variable"), all.x = TRUE); write_tab(NS, STAGE, "null_split_summary")
# null-calibrated p of the real threshold differences (50/50 null for sex, 20/80 for the other pairs)
TP[, prob_null := fifelse(pair == "Female/Male", .5, .2)]
TP[, `:=`(cold_p_null = mapply(function(pp, x) mean(abs(NUL[prob == pp]$thr_cold_diff) >= abs(x)), prob_null, cold_diff),
          heat_p_null = mapply(function(pp, x) mean(abs(NUL[prob == pp]$thr_heat_diff) >= abs(x)), prob_null, heat_diff))]
write_tab(TP, STAGE, "tm_pair_tests")
options(width = 220); print(TG[, .(group, cold = sprintf("P%g (%g-%g) %.1fC", cold_thr_pct, cold_lo, cold_hi, cold_thr_C), cold_slope = round(cold_slope, 2), heat = sprintf("P%g %.1fC", heat_thr_pct, heat_thr_C), heat_slope = round(heat_slope, 2), p_heat = round(p_heat_slope, 3))])
print(TP[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 3) else v)]); print(NS[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 2) else v)])
