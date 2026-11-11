# =============================================================================
# 01_mmt_and_significance.R - group thresholds read off the pooled two-stage curves, and the
# population at risk they imply (the Cao et al. 2026 comparison).
#   thresholds, all deaths and each margin (cold: cumulative lag 0-21 curve; heat: lag 0-2 curve):
#     MMT                    minimum of the curve on P1-P99 (Tobias et al. 2017)
#     significance threshold coldest / warmest uninterrupted run where RR vs the MMT has its lower
#                            95% bound > 1 (cold: from P1 upwards; heat: from P99 downwards)
#   95% intervals: 1,000 draws of the pooled coefficients (MMT re-located in each draw).
#   noise benchmark: the margin's curve redrawn around the ALL-DEATHS curve with the margin's own
#     covariance (identical true curves) -> spread expected from sampling noise alone.
#   population at risk: person-days beyond the margin's own threshold vs the all-deaths threshold,
#     Delta = (PD_own - PD_common) / PD_common, per margin and summed over each partition, with the
#     same quantity under the noise benchmark.
# Input : models/pooled_Tmin.rds; derived/regional_populations.csv
# Output: tables/06_thresholds/curve_thresholds.csv, curve_threshold_population_at_risk.csv,
#         curve_threshold_partition_totals.csv
# =============================================================================
source("2. Scripts/00_setup/config.R"); source("2. Scripts/00_setup/dlnm_tools.R"); source("2. Scripts/00_setup/threshold_tools.R")
STAGE <- "06_thresholds"; set.seed(SEED); NS <- 1000
PL <- read_model("pooled_Tmin"); panel <- load_panel()
rp <- fread(file.path(DERIVED, "regional_populations.csv"))[Year %in% YEARS]; EX <- exposure_tools(panel, rp)
SIDES <- list(cold = list(side = "cold", window = "cum"), heat = list(side = "heat", window = "acute"))
q95 <- function(v) quantile(v, c(.025, .975), na.rm = TRUE)

# ---- thresholds, draws and noise benchmark ---------------------------------------------------------------------
TH <- list(); DR <- list()
for (sn in names(SIDES)) { sp <- SIDES[[sn]]; pa <- PL[[paste("All", sp$window, sep = "|")]]
  for (g in GROUPS) { p <- PL[[paste(g, sp$window, sep = "|")]]; B <- p$basis
    pt <- thr_set(B, p$coef, p$vcov, sp$side)
    dr <- t(sapply(seq_len(NS), function(s) thr_set(B, p$sims[s, ], p$vcov, sp$side)))
    nz <- if (g == "All") dr else { cfn <- mvrnorm(NS, pa$coef, p$vcov); t(sapply(seq_len(NS), function(s) thr_set(B, cfn[s, ], p$vcov, sp$side))) }
    DR[[paste(sn, g)]] <- list(own = dr, noise = nz)
    for (k in c("mmt", "sig")) {
      pn <- if (g == "All" || is.na(pt[[k]])) NA_real_ else { v <- nz[, k]; v[is.na(v)] <- if (sp$side == "cold") 0 else 1; m <- median(v); mean(abs(v - m) >= abs(pt[[k]] - m)) }
      TH[[paste(sn, g, k)]] <- data.table(side = sn, group = g, threshold = k, pct = pt[[k]], lo = q95(dr[, k])[[1]], hi = q95(dr[, k])[[2]],
        share_draws_with_threshold = mean(!is.na(dr[, k])), noise_lo = q95(nz[, k])[[1]], noise_hi = q95(nz[, k])[[2]], p_vs_noise = pn) } }
  cat(sn, "thresholds done\n") }
TH <- rbindlist(TH)[, `:=`(degC = mapply(EX$pct_C, group, pct), degC_lo = mapply(EX$pct_C, group, lo), degC_hi = mapply(EX$pct_C, group, hi))]
write_tab(TH, STAGE, "curve_thresholds")

# ---- population at risk -----------------------------------------------------------------------------------------------
rel <- function(a, b) if (b > 0) a / b - 1 else NA_real_
PR <- rbindlist(lapply(names(SIDES), function(sn) { side <- SIDES[[sn]]$side
  rbindlist(lapply(c("mmt", "sig"), function(k) { cA <- TH[side == sn & group == "All" & threshold == k]$pct; dA <- DR[[paste(sn, "All")]]$own[, k]
    rbindlist(lapply(setdiff(GROUPS, "All"), function(g) { own <- TH[side == sn & group == g & threshold == k]$pct; po <- EX$pd_pct(g, own, side); pc <- EX$pd_pct(g, cA, side)
      dg <- DR[[paste(sn, g)]]$own[, k]; dn <- DR[[paste(sn, g)]]$noise[, k]
      r_own <- sapply(seq_len(NS), function(s) rel(EX$pd_pct(g, dg[s], side), EX$pd_pct(g, dA[s], side)))
      r_noise <- sapply(seq_len(NS), function(s) rel(EX$pd_pct(g, dn[s], side), pc))
      data.table(side = sn, rule = k, group = g, own_pct = own, common_pct = cA, share_own = po / EX$PY[[g]], share_common = pc / EX$PY[[g]],
                 delta_pct = 100 * rel(po, pc), delta_lo = 100 * quantile(r_own, .025, na.rm = TRUE), delta_hi = 100 * quantile(r_own, .975, na.rm = TRUE),
                 noise_lo = 100 * quantile(r_noise, .025, na.rm = TRUE), noise_hi = 100 * quantile(r_noise, .975, na.rm = TRUE)) })) })) }))
write_tab(PR, STAGE, "curve_threshold_population_at_risk")
TOT <- rbindlist(lapply(names(SIDES), function(sn) { side <- SIDES[[sn]]$side; rbindlist(lapply(c("mmt", "sig"), function(k) { cA <- TH[side == sn & group == "All" & threshold == k]$pct
  rbindlist(lapply(names(PARTITIONS), function(nm) { gs <- PARTITIONS[[nm]]; ec <- sum(sapply(gs, EX$pd_pct, p = cA, side = side)); if (ec == 0) return(NULL)
    eo <- sum(sapply(gs, function(g) EX$pd_pct(g, TH[side == sn & group == g & threshold == k]$pct, side)))
    v <- 100 * (sapply(seq_len(NS), function(s) sum(sapply(gs, function(g) EX$pd_pct(g, DR[[paste(sn, g)]]$noise[s, k], side)))) / ec - 1)
    data.table(side = sn, rule = k, partition = nm, E_common_millions = ec / 1e6, E_own_millions = eo / 1e6, delta_pct = 100 * (eo / ec - 1),
               noise_median = median(v), noise_lo = quantile(v, .025), noise_hi = quantile(v, .975)) })) })) }))
write_tab(TOT, STAGE, "curve_threshold_partition_totals")
print(TH[, .(side, group, threshold, pct = 100 * pct, degC = round(degC, 1), lo = 100 * lo, hi = 100 * hi, p_vs_noise = round(p_vs_noise, 2))], nrows = 50)
print(TOT[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 1) else v)])
