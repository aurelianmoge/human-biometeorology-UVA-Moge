# =============================================================================
# 02_real_vs_null.R - summaries of the null simulation and the real estimates inside it.
#   (1) calibration: false-positive rates of the whole-curve and pointwise tests;
#   (2) by group size: spread of group MMTs, Cao "hot days" and significance-based cold days;
#   (3) real Virginia statistics vs the null distributions (empirical two-sided p);
#   (4) SE calibration: SD of log RRR across replicates / median analytic SE, used to
#       correct the minimum detectable RRR and the years needed to detect 1.10.
# Input : tables/07_methods_supplement/null_*.csv (script 01); stage 03 outputs.
# Output: tables/07_methods_supplement/null_calibration.csv, null_by_size.csv,
#         real_vs_null_groups.csv, real_vs_null_pairs.csv, power_calibrated.csv
# =============================================================================
source("2. Scripts/00_setup/config.R"); source("2. Scripts/00_setup/dlnm_tools.R")
STAGE <- "07_methods_supplement"
f <- out_path("tables", STAGE, "null_pairs.csv")
if (!file.exists(f)) { message("No null simulation found (run 01_null_simulation.R); skipping."); quit(save = "no") }
G <- fread(out_path("tables", STAGE, "null_groups.csv")); Pn <- fread(f)
PL <- read_model("pooled_Tmin"); P <- function(g) PL[[paste(g, "cum", sep = "|")]]
map <- c(Black = "sB", White = "sW", Female = "sF", Male = "sM", Under50 = "sU50", Age50_74 = "s5074", Age75plus = "s75")
pairs <- list(race = c("Black", "White"), sex = c("Female", "Male"), age_old_young = c("Age75plus", "Under50"), age_old_mid = c("Age75plus", "Age50_74"))
q95 <- function(x) paste(round(quantile(x, c(.025, .975), na.rm = TRUE), 3), collapse = " to ")

calib <- Pn[, .(replicates = .N, false_pos_curve = mean(p_curve < .05, na.rm = TRUE), false_pos_P1 = mean(p_p1 < .05), false_pos_P99 = mean(p_p99 < .05),
                null95_RRR_P1 = q95(rrr_p1), null95_own_RRR_P1 = q95(own_rrr_p1), abs_dMMT_90pct = 100 * quantile(abs(dmmt), .9)), by = pair]
write_tab(calib, STAGE, "null_calibration")
size <- G[, .(mmt_q05 = 100 * quantile(mmt, .05), mmt_median = 100 * median(mmt), mmt_q95 = 100 * quantile(mmt, .95), mmt_at_boundary = mean(mmt <= .011 | mmt >= .989),
              hot_days_q05 = 100 * quantile(hot_days, .05), hot_days_q95 = 100 * quantile(hot_days, .95), sig_cold_mean = 100 * mean(sig_cold)), by = .(group, size)][order(size)]
write_tab(size, STAGE, "null_by_size")

# ---- real statistics ---------------------------------------------------------------------------------------
w_pct <- c(diff(PCT), 0)
real_g <- rbindlist(lapply(names(map), function(g) { p <- P(g); rc <- rr_curve(p, p$mmt); lo <- rc$logrr - 1.96 * rc$se
  data.table(group = g, sim = map[[g]], mmt = p$mmt, own_p1 = lrr(p, .01, p$mmt)[[1]], own_p99 = lrr(p, .99, p$mmt)[[1]],
             sig_cold = sum(w_pct[lo > 0 & PCT < p$mmt]), sig_heat = sum(w_pct[lo > 0 & PCT > p$mmt])) }))
emp <- function(null, obs) { m <- median(null, na.rm = TRUE); mean(abs(null - m) >= abs(obs - m), na.rm = TRUE) }
cg <- rbindlist(lapply(seq_len(nrow(real_g)), function(i) { r <- real_g[i]; n <- G[group == r$sim]
  rbindlist(lapply(c("mmt", "sig_cold", "sig_heat", "own_p1", "own_p99"), function(v) data.table(group = r$group, statistic = v, observed = r[[v]],
    null_median = median(n[[v]]), null_q025 = quantile(n[[v]], .025), null_q975 = quantile(n[[v]], .975), empirical_p = emp(n[[v]], r[[v]])))) }))
write_tab(cg, STAGE, "real_vs_null_groups")
S1 <- read_model("stage1_Tmin"); ref <- P("All")$mmt
real_p <- rbindlist(lapply(names(pairs), function(nm) { a <- pairs[[nm]][1]; b <- pairs[[nm]][2]; ra <- real_g[group == a]; rb <- real_g[group == b]
  data.table(pair = nm, p_curve = curve_test(S1, pairs[[nm]])[["p"]], rrr_p1 = rr_ratio(P(a), P(b), .01, ref)[["RRR"]], rrr_p99 = rr_ratio(P(a), P(b), .99, ref)[["RRR"]],
             own_rrr_p1 = exp(ra$own_p1 - rb$own_p1), own_rrr_p99 = exp(ra$own_p99 - rb$own_p99), dmmt = ra$mmt - rb$mmt, dsig_cold = ra$sig_cold - rb$sig_cold) }))
cp <- rbindlist(lapply(seq_len(nrow(real_p)), function(i) { r <- real_p[i]; n <- Pn[pair == r$pair]
  rbindlist(lapply(setdiff(names(real_p), "pair"), function(v) data.table(pair = r$pair, statistic = v, observed = r[[v]], null_median = median(n[[v]], na.rm = TRUE),
    null_q025 = quantile(n[[v]], .025, na.rm = TRUE), null_q975 = quantile(n[[v]], .975, na.rm = TRUE),
    empirical_p = if (v == "p_curve") mean(n[[v]] <= r[[v]], na.rm = TRUE) else emp(n[[v]], r[[v]])))) }))
write_tab(cp, STAGE, "real_vs_null_pairs")

# ---- SE calibration and corrected power -------------------------------------------------------------------------
cal <- Pn[is.finite(p_p1) & p_p1 > 0 & p_p1 < 1, .(se_ratio = sd(log(rrr_p1)) / median(abs(log(rrr_p1)) / qnorm(1 - p_p1 / 2))), by = pair]
pw <- fread(out_path("tables", "03_temperature", "power_minimum_detectable.csv"))[pct == 1]
pw[, sim := c("Female/Male" = "sex", "Black/White" = "race", "Age75plus/Under50" = "age_old_young", "Age75plus/Age50_74" = "age_old_mid")[pair]]
pw <- merge(pw[!is.na(sim)], cal, by.x = "sim", by.y = "pair")
pw[, se := log(MDR80) / 2.80][, `:=`(MDR80_calibrated = exp(2.80 * se * se_ratio), years_calibrated = 16 * (2.80 * se * se_ratio / log(1.10))^2)]
write_tab(pw[, .(pair, RRR, lo, hi, MDR80, MDR80_calibrated, years_to_detect_1.10, years_calibrated, se_ratio)], STAGE, "power_calibrated")
print(calib); print(size); print(pw[, .(pair, MDR80 = round(MDR80, 3), MDR80_calibrated = round(MDR80_calibrated, 3), years = round(years_to_detect_1.10), years_calibrated = round(years_calibrated))])
