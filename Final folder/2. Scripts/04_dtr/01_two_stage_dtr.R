# =============================================================================
# 01_two_stage_dtr.R - diurnal temperature range (DTR = Tmax - Tmin) and mortality.
# Same two-stage design as the temperature model; DTR gets its own cross-basis
# (ns knots at regional P10/P50/P90, lag 0-21) and the model ALWAYS includes a full
# cross-basis of Tmean, so DTR effects are net of temperature level (Davis et al. 2025).
# Reference = regional median DTR.
# Variants
#   base      cumulative 0-21 and acute 0-2, all months
#   cold / warm   outcome days restricted to Nov-Mar / May-Sep (exposure histories kept)
#   humidity  + ns(7 am RH, 3) + ns(7 am dew point, 3)      (low DTR = cloudy / humid days)
#   anomaly   DTR minus its region x calendar-month mean     (removes the DTR seasonal cycle)
#   lag7      lag 0-7
# Also: lag-response at the regional P10 of DTR by sex, and region-specific Female/Male ratios.
# Output: models/pooled_DTR.rds (base, cumulative; used by stage 06); tables/04_dtr/
# =============================================================================
source("2. Scripts/00_setup/config.R"); source("2. Scripts/00_setup/dlnm_tools.R")
STAGE <- "04_dtr"; set.seed(SEED)
panel <- load_panel()
panel[, DTRa_f := DTR_f - mean(DTR_f, na.rm = TRUE), by = .(Region, month)]
variants <- list(base = list(), cold = list(months = c(11, 12, 1, 2, 3)), warm = list(months = 5:9),
                 humidity = list(extra = c("ns(RH07_f, 3)", "ns(Td07_f, 3)")), anomaly = list(expo = "DTRa_f"), lag7 = list(lag = 7))
QS <- c(.01, .05, .10, .25, .75, .90, .95, .99)
rows <- list(); tests <- list(); s1_base <- NULL
for (nm in names(variants)) {
  v <- variants[[nm]]; d <- copy(panel)
  if (!is.null(v$months)) for (g in GROUPS) d[!(month %in% v$months), (g) := NA]
  lag <- if (is.null(v$lag)) MAX_LAG else v$lag
  s1 <- run_stage1(d, if (is.null(v$expo)) "DTR_f" else v$expo, lag = lag, vtype = "nsm", control = "Tmean_f", extra = v$extra,
                   windows = list(cum = c(0, NA), acute = c(0, 2)), lag_at = if (nm == "base") .10 else NULL)
  if (nm == "base") s1_base <- s1
  for (w in if (nm == "base") c("cum", "acute") else "cum") {
    pl <- setNames(lapply(GROUPS, function(g) pool(s1, g, w, "nsm", nsim = if (nm == "base" && w == "cum") 1000 else 0)), GROUPS)
    if (nm == "base" && w == "cum") save_model(pl, "pooled_DTR")
    for (g in GROUPS) for (q in QS) { l <- lrr(pl[[g]], q, .5)
      rows[[length(rows) + 1]] <- data.table(variant = nm, window = w, group = g, pct = 100 * q, DTR = pl[[g]]$tavg[pidx(q)],
        RR = exp(l[1]), lo = exp(l[1] - 1.96 * l[2]), hi = exp(l[1] + 1.96 * l[2]), I2 = pl[[g]]$i2) }
    for (pr in PAIRS) { ct <- curve_test(s1, pr, w)
      for (q in QS) { r <- rr_ratio(pl[[pr[1]]], pl[[pr[2]]], q, .5)
        tests[[length(tests) + 1]] <- data.table(variant = nm, window = w, pair = paste(pr, collapse = "/"), pct = 100 * q,
          RRR = r[["RRR"]], lo = r[["lo"]], hi = r[["hi"]], p_point = r[["p"]], p_curve = ct[["p"]]) } }
  }
  cat(nm, "done\n")
}
write_tab(rbindlist(rows), STAGE, "dtr_RR"); tests <- rbindlist(tests); write_tab(tests, STAGE, "dtr_contrasts")

# ---- lag-response at P10 DTR, by sex (pooled) -----------------------------------------------------
lb <- onebasis(0:MAX_LAG, fun = "ns", knots = logknots(MAX_LAG, 3), intercept = TRUE)
lagres <- rbindlist(lapply(c("All", "Female", "Male"), function(g) { s <- Filter(function(z) z$group == g, s1_base)
  m <- mixmeta(t(sapply(s, function(z) z$red$lag_response$coef)) ~ 1, S = lapply(s, function(z) z$red$lag_response$vcov), method = "reml")
  est <- drop(lb %*% coef(m)); se <- sqrt(rowSums((lb %*% vcov(m)) * lb))
  data.table(group = g, lag = 0:MAX_LAG, RR = exp(est), lo = exp(est - 1.96 * se), hi = exp(est + 1.96 * se)) }))
write_tab(lagres, STAGE, "dtr_lag_response_P10_by_sex")

# ---- region-specific Female / Male ratio at the regional P10 of DTR ------------------------------------
reg <- rbindlist(lapply(sort(unique(panel$Region)), function(r) {
  one <- function(g) { z <- Filter(function(z) z$region == r && z$group == g, s1_base)[[1]]
    b <- onebasis(z$pct, fun = "ns", knots = z$pct[sapply(c(.1, .5, .9), pidx)], Boundary.knots = range(z$pct))
    x <- b[pidx(.10), ] - b[pidx(.5), ]; c(sum(x * z$red$cum$coef), sqrt(drop(t(x) %*% z$red$cum$vcov %*% x))) }
  f <- one("Female"); m <- one("Male"); d <- f[1] - m[1]; se <- sqrt(f[2]^2 + m[2]^2)
  data.table(region = r, RRR = exp(d), lo = exp(d - 1.96 * se), hi = exp(d + 1.96 * se)) }))
write_tab(reg, STAGE, "dtr_female_male_by_region")
print(tests[pair == "Female/Male" & pct %in% c(10, 99), .(variant, window, pct, RRR = round(RRR, 3), lo = round(lo, 3), hi = round(hi, 3), p_curve = signif(p_curve, 2))])
cat("Regions with Female/Male < 1 at P10:", sum(reg$RRR < 1), "/", nrow(reg), "\n")
