# =============================================================================
# 04_risk_level.R - the temperature at which mortality is raised by at least k% (k = 5, 10), by
# group, and what it changes for the definition of days at risk.
#   From the threshold models of 03: beyond the threshold the cumulative log RR is b x degrees, so
#   excess risk k is reached at  T_k,r = tau_r - ln(1 + k) / b  (cold)  or  tau_r + ln(1 + k) / b  (heat).
#   (a) T5, T10 by group (population-weighted deg C), 95% intervals from joint draws of (threshold,
#       slope): q ~ profile-likelihood weights, b | q ~ N(b_hat(q), se(q)^2); b <= 0 = never reached.
#   (b) pair tests: T10 difference (draws) and slope difference at the all-deaths thresholds
#       (z-test), each calibrated on the random-split null of 03 (50/50 for sex, 20/80 otherwise).
#   (c) days at risk at each level: person-days below the group's OWN T_k vs below the all-deaths
#       T_k; "missed" days and deaths = days on which the group's own excess risk is >= k but the
#       population rule does not flag them, with the group's temperature-attributable deaths on those
#       days (forward attributable deaths from its own threshold model, Gasparrini & Leone 2014);
#       "extra" days = flagged by the population rule although the group is below k.
# Input : models/threshold_models.rds; tables/06_thresholds/null_split_replicates.csv
# Output: tables/06_thresholds/risk_level_temperatures.csv, risk_level_pair_tests.csv,
#         risk_level_days_at_risk.csv, risk_level_partition_totals.csv
# =============================================================================
source("2. Scripts/00_setup/config.R"); source("2. Scripts/00_setup/dlnm_tools.R"); source("2. Scripts/00_setup/threshold_tools.R")
STAGE <- "06_thresholds"; set.seed(SEED); ND <- 2000; KS <- c(.05, .10)
panel <- load_panel(); setorder(panel, Region, Date)
rp <- fread(file.path(DERIVED, "regional_populations.csv"))[Year %in% YEARS]; EX <- exposure_tools(panel, rp)
S <- tm_setup(panel); TM <- read_model("threshold_models"); PR <- TM$models
NUL <- fread(out_path("tables", STAGE, "null_split_replicates.csv"))
q_of <- function(p, side) if (side == "cold") p$qc else p$qh
b_of <- function(p, side) if (side == "cold") p$bc else p$bh
DR <- lapply(PR, function(p) list(cold = tm_draws_qb(p, "cold", ND), heat = tm_draws_qb(p, "heat", ND)))
tk_draws <- function(g, side, k) sapply(seq_len(ND), function(s) EX$temp_C(g, risk_level_temp(S, DR[[g]][[side]]$q[s], DR[[g]][[side]]$b[s], k, side)))

# ---- (a) risk-level temperatures -------------------------------------------------------------------------------------
RL <- rbindlist(lapply(GROUPS, function(g) { p <- PR[[g]]; rbindlist(lapply(c("cold", "heat"), function(side) rbindlist(lapply(KS, function(k) { v <- tk_draws(g, side, k)
  data.table(group = g, side = side, k = 100 * k, onset_C = EX$temp_C(g, TAU(S, q_of(p, side))), slope_pct = 100 * (exp(b_of(p, side)) - 1),
             T_k = EX$temp_C(g, risk_level_temp(S, q_of(p, side), b_of(p, side), k, side)), lo = quantile(v, .025), hi = quantile(v, .975), share_reached = mean(is.finite(v))) })))) }))
write_tab(RL, STAGE, "risk_level_temperatures")

# ---- (b) pair tests -----------------------------------------------------------------------------------------------------
PT <- rbindlist(lapply(PAIRS, function(pr) rbindlist(lapply(c("cold", "heat"), function(side) { a <- PR[[pr[1]]]; b <- PR[[pr[2]]]; S_ <- side
  d <- tk_draws(pr[1], side, .10) - tk_draws(pr[2], side, .10)
  pt <- RL[group == pr[1] & side == S_ & k == 10]$T_k - RL[group == pr[2] & side == S_ & k == 10]$T_k
  ba <- if (side == "cold") a$common$bc else a$common$bh; bb <- if (side == "cold") b$common$bc else b$common$bh
  sa <- if (side == "cold") a$common$sc else a$common$sh; sb <- if (side == "cold") b$common$sc else b$common$sh
  pn <- if (pr[1] == "Female") .5 else .2; nv <- NUL[prob == pn][[paste0("T10_", side, "_diff")]]; nv <- nv[is.finite(nv)]
  data.table(pair = paste(pr, collapse = "/"), side = side, T10_diff = pt, T10_diff_lo = quantile(d, .025, na.rm = TRUE), T10_diff_hi = quantile(d, .975, na.rm = TRUE),
             p_T10_draws = 2 * min(mean(d > 0, na.rm = TRUE), mean(d < 0, na.rm = TRUE)), p_T10_null = if (is.finite(pt)) mean(abs(nv) >= abs(pt)) else NA_real_,
             slope_A_common = 100 * (exp(ba) - 1), slope_B_common = 100 * (exp(bb) - 1), slope_ratio = (exp(ba) - 1) / (exp(bb) - 1),
             p_slope_common = 2 * pnorm(-abs((ba - bb) / sqrt(sa^2 + sb^2))), slope_test_false_positive_null = mean(NUL[prob == pn][[paste0("p_slope_", side)]] < .05)) }))))
write_tab(PT, STAGE, "risk_level_pair_tests")

# ---- (c) days at risk, missed days and deaths ---------------------------------------------------------------------------------
an_fwd <- function(g, p, side) { t <- TAU(S, q_of(p, side))[match(panel$Region, S$regs)]
  x <- if (side == "cold") pmax(t - panel$Tmin_f, 0) else pmax(panel$Tmin_f - t, 0)
  fm <- panel[, .(f = frollmean(get(g), if (side == "cold") MAX_LAG + 1 else 3, align = "left", na.rm = TRUE)), by = Region]$f
  ifelse(is.na(x) | is.na(fm), 0, (1 - exp(-b_of(p, side) * x)) * fm) }
DAR <- rbindlist(lapply(KS, function(kk) rbindlist(lapply(c("cold", "heat"), function(side) { pa <- PR$All; common <- risk_level_temp(S, q_of(pa, side), b_of(pa, side), kk, side)
  rbindlist(lapply(setdiff(GROUPS, "All"), function(g) { p <- PR[[g]]; own <- risk_level_temp(S, q_of(p, side), b_of(p, side), kk, side)
    an <- an_fwd(g, p, side); fo <- EX$pd_temp(g, own, side, flag = TRUE); fc <- EX$pd_temp(g, common, side, flag = TRUE); miss <- fo & !fc; extra <- fc & !fo; pop <- EX$pop[, g]
    rel <- sapply(seq_len(1000), function(s) { o <- EX$pd_temp(g, risk_level_temp(S, DR[[g]][[side]]$q[s], DR[[g]][[side]]$b[s], kk, side), side)
      cc <- EX$pd_temp(g, risk_level_temp(S, DR$All[[side]]$q[s], DR$All[[side]]$b[s], kk, side), side); if (cc > 0) o / cc - 1 else NA })
    data.table(k = 100 * kk, side = side, group = g, own_Tk = EX$temp_C(g, own), common_Tk = EX$temp_C(g, common),
               share_days_own = 100 * sum(pop[fo]) / EX$PY[[g]], share_days_common = 100 * sum(pop[fc]) / EX$PY[[g]],
               delta_pct = 100 * (sum(pop[fo]) / max(sum(pop[fc]), 1) - 1), delta_lo = 100 * quantile(rel, .025, na.rm = TRUE), delta_hi = 100 * quantile(rel, .975, na.rm = TRUE),
               share_attr_own = 100 * sum(an[fo]) / sum(an), share_attr_common = 100 * sum(an[fc]) / sum(an),
               missed_days_pct = 100 * sum(pop[miss]) / EX$PY[[g]], missed_attr_pct = 100 * sum(an[miss]) / sum(an), extra_days_pct = 100 * sum(pop[extra]) / EX$PY[[g]]) })) })) ))
write_tab(DAR, STAGE, "risk_level_days_at_risk")
DTOT <- rbindlist(lapply(KS, function(kk) { K <- 100 * kk; rbindlist(lapply(c("cold", "heat"), function(side) { S_ <- side
  rbindlist(lapply(names(PARTITIONS), function(nm) { x <- DAR[k == K & side == S_ & group %in% PARTITIONS[[nm]]]; w <- EX$PY[x$group]
    data.table(k = K, side = side, partition = nm, delta_pct = 100 * (sum(x$share_days_own * w) / sum(x$share_days_common * w) - 1),
               missed_days_pct = sum(x$missed_days_pct * w) / sum(w), extra_days_pct = sum(x$extra_days_pct * w) / sum(w)) })) })) }))
write_tab(DTOT, STAGE, "risk_level_partition_totals")
options(width = 220); print(RL[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 2) else v)], nrows = 40); print(PT[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 3) else v)])
print(DAR[, .(k, side, group, own = round(own_Tk, 1), common = round(common_Tk, 1), days_own = round(share_days_own, 1), days_common = round(share_days_common, 1), missed_days = round(missed_days_pct, 1), missed_attr = round(missed_attr_pct, 1), extra = round(extra_days_pct, 1))], nrows = 40)
