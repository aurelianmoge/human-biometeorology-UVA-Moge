# =============================================================================
# 01_published_design.R - the design of DeGuzman et al. (2026), applied to every margin.
#   * one statewide daily series; exposure = mortality-weighted mean of the 12 regional
#     stations (DeGuzman 2026, s.2); quasi-Poisson DLNM, lag 0-21 (ns knots P10/P50/P90),
#     ns(time, 3 df/yr; 7 df/yr as sensitivity), day of week;
#   * Tmin: RR at -10, -5, 0 C (cold) and 22, 24 C (heat) vs 5 C (their reference);
#   * DTR: RR at P1, P5, P95, P99 vs the median, cumulative 0-21 and acute 0-2, adjusted for a
#     Tmin cross-basis;
#   * reported as in the literature (each group's RR with CI; excess-risk ratio "x times")
#     plus the direct test of the RR ratio (groups = disjoint series, independent SEs).
#   * geography: Black / White ratio implied by WHERE deaths occur, using region-specific
#     (BLUP) all-cause curves weighted by each race's regional share of deaths.
#   * exact replication of DeGuzman et al.'s exposure basis (their Methods): quadratic spline with
#     two equally spaced knots, lag natural spline with two log-spaced knots, ns(time, 3 df/yr),
#     RRs vs 5 C at -12, -10, 0, 22, 24 C, compared with their published Black / White RRs.
# Output: tables/05_group_differences/published_design_RR.csv, published_design_contrasts.csv,
#         geography_expectation.csv, deguzman_exact_RR.csv, deguzman_exact_contrasts.csv
# =============================================================================
source("2. Scripts/00_setup/config.R"); source("2. Scripts/00_setup/dlnm_tools.R")
STAGE <- "05_group_differences"
panel <- load_panel()
st <- panel[, c(lapply(.SD, function(v) if (anyNA(v)) NA_real_ else as.numeric(sum(v))),
                list(Tmin_w = weighted.mean(Tmin_f, fifelse(is.na(All), 0, All)), DTR_w = weighted.mean(DTR_f, fifelse(is.na(All), 0, All)))),
            by = Date, .SDcols = GROUPS][order(Date)]
st[, `:=`(time = seq_len(.N), dow = factor(wday(Date)), year = year(Date))]

fit <- function(y, expo, dfy, control = NULL) {
  x <- st[[expo]]; cb <- crossbasis(x, lag = MAX_LAG, argvar = argvar_spec(x, "nsm"), arglag = list(knots = logknots(MAX_LAG, 3)))
  rhs <- c("cb", sprintf("ns(time, df = %d)", dfy * 16), "dow")
  if (!is.null(control)) { z <- st[[control]]; cc <- crossbasis(z, lag = MAX_LAG, argvar = argvar_spec(z, "ns3"), arglag = list(knots = logknots(MAX_LAG, 3))); rhs <- c(rhs, "cc") }
  m <- glm(as.formula(paste("yy ~", paste(rhs, collapse = " + "))), data = cbind(st, yy = st[[y]]), family = quasipoisson, na.action = "na.exclude")
  list(cb = cb, m = m)
}
pred <- function(f, at, cen, lag = NULL) { cb <- f$cb                         # dlnm matches coefficients by the basis name
  if (is.null(lag)) { cp <- crosspred(cb, f$m, at = at, cen = cen); return(data.table(est = unname(cp$allfit), se = unname(cp$allse))) }
  r <- crossreduce(cb, f$m, type = "overall", lag = lag, at = at, cen = cen); data.table(est = unname(r$fit), se = unname(r$se)) }

specs <- list(list(expo = "Tmin", var = "Tmin_w", ctl = NULL, at = c(-10, -5, 0, 22, 24), cen = 5, lags = list(cum = NULL)),
              list(expo = "DTR", var = "DTR_w", ctl = "Tmin_w", at = NULL, cen = NULL, lags = list(cum = NULL, acute = c(0, 2))))
rows <- list()
for (sp in specs) for (dfy in c(3, 7)) {
  fits <- lapply(setNames(GROUPS, GROUPS), fit, expo = sp$var, dfy = dfy, control = sp$ctl)
  x <- st[[sp$var]]; at <- if (is.null(sp$at)) unname(quantile(x, c(.01, .05, .95, .99))) else sp$at
  cen <- if (is.null(sp$cen)) median(x) else sp$cen
  pts <- if (is.null(sp$at)) c("P1", "P5", "P95", "P99") else paste0(at, "C")
  for (w in names(sp$lags)) for (g in GROUPS) rows[[length(rows) + 1]] <- data.table(exposure = sp$expo, df_year = dfy, window = w, group = g,
    point = pts, value = at, reference = cen, pred(fits[[g]], at, cen, sp$lags[[w]]))
  cat(sp$expo, dfy, "df/yr done\n")
}
R <- rbindlist(rows)[, `:=`(RR = exp(est), lo = exp(est - 1.96 * se), hi = exp(est + 1.96 * se))]
write_tab(R, STAGE, "published_design_RR")
C <- rbindlist(lapply(PAIRS, function(pr) {
  m <- merge(R[group == pr[1]], R[group == pr[2]], by = c("exposure", "df_year", "window", "point", "value", "reference"), suffixes = c("_a", "_b"))
  m[, `:=`(pair = paste(pr, collapse = "/"), d = est_a - est_b, sd = sqrt(se_a^2 + se_b^2))]
  m[, .(exposure, df_year, window, point, value, pair, RR_a, lo_a, hi_a, RR_b, lo_b, hi_b, excess_ratio = (RR_a - 1) / (RR_b - 1),
        RRR = exp(d), RRR_lo = exp(d - 1.96 * sd), RRR_hi = exp(d + 1.96 * sd), p_RRR = 2 * pnorm(-abs(d / sd)))] }))
write_tab(C, STAGE, "published_design_contrasts")

# ---- geography-only Black / White expectation (regional all-cause BLUP curves) -------------------------------
pa <- read_model("pooled_Tmin")[["All|cum"]]
v <- panel[, .(Black = sum(Black, na.rm = TRUE), White = sum(White, na.rm = TRUE)), by = Region][, `:=`(vB = Black / sum(Black), vW = White / sum(White))]
geo <- rbindlist(lapply(c(.01, .025, .10), function(q) { X <- pa$basis[pidx(q), ] - pa$basis[pidx(pa$mmt), ]
  rr <- sapply(v$Region, function(r) exp(sum(X * pa$blup[[r]]$blup)))
  data.table(pct = 100 * q, RR_Black_geo = sum(v$vB * rr), RR_White_geo = sum(v$vW * rr), RRR_geo = sum(v$vB * rr) / sum(v$vW * rr)) }))
write_tab(geo, STAGE, "geography_expectation")
print(C[df_year == 3 & p_RRR < .1, .(exposure, window, point, pair, RR_a = round(RR_a, 3), RR_b = round(RR_b, 3), RRR = round(RRR, 3), p = signif(p_RRR, 2))])

# ---- exact DeGuzman et al. (2026) exposure basis ---------------------------------------------------------------
x <- st$Tmin_w; kn <- seq(min(x), max(x), length.out = 4)[2:3]
cbE <- crossbasis(x, lag = MAX_LAG, argvar = list(fun = "bs", degree = 2, knots = kn), arglag = list(knots = logknots(MAX_LAG, 2)))
TX <- c(-12, -10, 0, 22, 24)
EX <- rbindlist(lapply(GROUPS, function(g) {
  m <- glm(st[[g]] ~ cbE + ns(time, df = 3 * 16) + dow, data = st, family = quasipoisson, na.action = "na.exclude")
  cp <- crosspred(cbE, m, at = TX, cen = 5)
  data.table(group = g, temp = TX, est = unname(cp$allfit), se = unname(cp$allse)) }))[, `:=`(RR = exp(est), lo = exp(est - 1.96 * se), hi = exp(est + 1.96 * se))]
pub <- data.table(group = c("Black", "White", "Black", "White"), temp = c(0, 0, -10, -10), published_RR = c(1.03, 1.04, 1.13, 1.02))
EX <- merge(EX, pub, by = c("group", "temp"), all.x = TRUE)
write_tab(EX, STAGE, "deguzman_exact_RR")
EC <- rbindlist(lapply(PAIRS, function(pr) { m <- merge(EX[group == pr[1]], EX[group == pr[2]], by = "temp", suffixes = c("_a", "_b"))
  m[, `:=`(d = est_a - est_b, sd = sqrt(se_a^2 + se_b^2))]
  m[, .(pair = paste(pr, collapse = "/"), temp, RR_a, lo_a, hi_a, RR_b, lo_b, hi_b, excess_ratio = (RR_a - 1) / (RR_b - 1),
        RRR = exp(d), RRR_lo = exp(d - 1.96 * sd), RRR_hi = exp(d + 1.96 * sd), p_RRR = 2 * pnorm(-abs(d / sd)))] }))
write_tab(EC, STAGE, "deguzman_exact_contrasts")
print(EX[group %in% c("Black", "White") & temp %in% c(-10, 0), .(group, temp, RR = round(RR, 3), lo = round(lo, 3), hi = round(hi, 3), published_RR)])
