# =============================================================================
# 02_group_contrasts.R - do the demographic margins differ in their temperature response?
#   (a) whole-curve test: multivariate meta-regression coef ~ group, random ~ 1 | region/group
#       (Gasparrini et al. 2012); within a region the two margins are disjoint deaths;
#   (b) ratio of relative risks (RRR) at P1, P2.5, P10, P90, P99 with a COMMON reference
#       (all-cause MMT; median as sensitivity);
#   (c) Black/White: expectation from age composition alone, using the age distribution of
#       Black and White deaths (DeGuzman et al. 2026, Table 1, same VDH files):
#                      <50    50-74   75+
#          Black      14.8%  43.1%  42.1%
#          White       8.1%  33.5%  58.4%
#   (d) dependence of the RRR at P1 on the reference temperature (P30-P99, own MMTs, 5 C);
#   (e) minimum detectable RRR at P1 (80% power, two-sided 5%) = exp(2.80 SE) and years of
#       data needed to detect 1.10 (SE proportional to 1/sqrt(years)).
# Input : models/stage1_*.rds, pooled_*.rds.   Output: tables/03_temperature/
# =============================================================================
source("2. Scripts/00_setup/config.R"); source("2. Scripts/00_setup/dlnm_tools.R")
STAGE <- "03_temperature"
S1 <- list(Tmin = read_model("stage1_Tmin"), Tmean = read_model("stage1_Tmean"))
PL <- list(Tmin = read_model("pooled_Tmin"), Tmean = read_model("pooled_Tmean"))
P <- function(e, g, w = "cum") PL[[e]][[paste(g, w, sep = "|")]]

# ---- (a) whole-curve tests ------------------------------------------------------------------------
tests <- rbindlist(lapply(names(S1), function(e) rbindlist(lapply(c("cum", "acute", "delayed"), function(w) rbindlist(lapply(PAIRS, function(pr) {
  ct <- curve_test(S1[[e]], pr, w); data.table(exposure = e, window = w, pair = paste(pr, collapse = "/"), wald = ct[["wald"]], df = ct[["df"]], p = ct[["p"]]) }))))))
write_tab(tests, STAGE, "wholecurve_tests")

# ---- (b) pointwise RRR ---------------------------------------------------------------------------------
rrr <- list()
for (e in names(PL)) for (w in c("cum", "acute", "delayed")) for (pr in PAIRS) for (cen in c("MMT_all", "median")) {
  ref <- if (cen == "median") .5 else P(e, "All")$mmt
  for (q in c(.01, .025, .10, .90, .99)) {
    a <- lrr(P(e, pr[1], w), q, ref); b <- lrr(P(e, pr[2], w), q, ref); r <- rr_ratio(P(e, pr[1], w), P(e, pr[2], w), q, ref)
    rrr[[length(rrr) + 1]] <- data.table(exposure = e, window = w, pair = paste(pr, collapse = "/"), centre = cen, pct = 100 * q,
      RR_a = exp(a[1]), RR_b = exp(b[1]), RRR = r[["RRR"]], lo = r[["lo"]], hi = r[["hi"]], p = r[["p"]], se = r[["se"]])
  }
}
rrr <- rbindlist(rrr); write_tab(rrr, STAGE, "pointwise_RRR")

# ---- (c) age-composition expectation for Black / White ------------------------------------------------
w_age <- rbind(Black = c(.148, .431, .421), White = c(.081, .335, .584)); ages <- c("Under50", "Age50_74", "Age75plus")
comp <- rbindlist(lapply(names(PL), function(e) rbindlist(lapply(c("cum", "acute", "delayed"), function(w) {
  ref <- P(e, "All")$mmt; qs <- c(.01, .025, .10, .90, .99)
  sims <- function(g) { p <- P(e, g, w); X <- sweep(p$basis[sapply(qs, pidx), , drop = FALSE], 2, p$basis[pidx(ref), ], "-"); exp(p$sims %*% t(X)) }
  A <- lapply(ages, sims); ex <- Reduce(`+`, Map(`*`, A, w_age["Black", ])) / Reduce(`+`, Map(`*`, A, w_age["White", ]))
  obs <- sims("Black") / sims("White"); exc <- obs / ex
  rbindlist(lapply(seq_along(qs), function(j) data.table(exposure = e, window = w, pct = 100 * qs[j], RRR_observed = median(obs[, j]),
    RRR_age_expected = median(ex[, j]), age_lo = quantile(ex[, j], .025), age_hi = quantile(ex[, j], .975),
    within_age_RRR = median(exc[, j]), within_lo = quantile(exc[, j], .025), within_hi = quantile(exc[, j], .975)))) }))))
write_tab(comp, STAGE, "race_age_composition")

# ---- (d) reference dependence of the RRR at P1 (Tmin, cumulative) -----------------------------------------
refdep <- rbindlist(lapply(PAIRS, function(pr) { a <- P("Tmin", pr[1]); b <- P("Tmin", pr[2]); p5 <- PCT[which.min(abs(a$tavg - 5))]
  one <- function(ra, rb, kind, refp) { la <- lrr(a, .01, ra); lb <- lrr(b, .01, rb); d <- la[1] - lb[1]; s <- sqrt(la[2]^2 + lb[2]^2)
    data.table(pair = paste(pr, collapse = "/"), kind = kind, ref_pct = refp, RRR = exp(d), lo = exp(d - 1.96 * s), hi = exp(d + 1.96 * s)) }
  rbind(rbindlist(lapply(c(seq(.30, .95, .05), .99), function(q) one(q, q, "common percentile", 100 * q))),
        one(a$mmt, b$mmt, "own MMTs (Cao)", NA), one(p5, p5, "5 C (DeGuzman)", 100 * p5)) }))
write_tab(refdep, STAGE, "reference_dependence")

# ---- (e) power ----------------------------------------------------------------------------------------------
pw <- rrr[exposure == "Tmin" & window == "cum" & centre == "MMT_all" & pct %in% c(1, 99)]
pw[, `:=`(MDR80 = exp(2.80 * se), years_to_detect_1.10 = 16 * (2.80 * se / log(1.10))^2)]
write_tab(pw[, .(pair, pct, RRR, lo, hi, MDR80, years_to_detect_1.10)], STAGE, "power_minimum_detectable")

options(width = 200)
print(tests[exposure == "Tmin", .(window, pair, wald = round(wald, 1), p = signif(p, 2))])
print(rrr[exposure == "Tmin" & window == "cum" & centre == "MMT_all" & pct %in% c(1, 99), .(pair, pct, RRR = round(RRR, 3), lo = round(lo, 3), hi = round(hi, 3), p = signif(p, 2))])
