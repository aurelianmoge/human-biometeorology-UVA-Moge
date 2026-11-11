# =============================================================================
# 04_burden_decomposition.R - why groups differ in temperature-attributable deaths.
#   For group g:  Rate_g = sum_r  w_g,r x m_g,r x AF_g,r   (per 100,000 person-years)
#     w  = share of g's person-years in region r        -> GEOGRAPHY (where the group lives)
#     m  = g's baseline death rate in region r            -> FRAILTY (baseline mortality, incl. age mix)
#     AF = share of g's deaths in r attributable to cold or heat, from g's own region-specific
#          (BLUP) curve; within a region all groups share the same weather -> SENSITIVITY
#   Cold: cumulative lag 0-21 curves; heat: lag 0-2 curves; both centred at the common all-cause
#   MMT; forward attributable deaths (Gasparrini & Leone 2014).  The gap between two groups is
#   split exactly into the three components by a Shapley decomposition (mean over the 6 orders
#   of substitution).  Uncertainty: 200 draws of the regional curves (w, m are complete counts).
#   Robustness: the same with the nearest-station crosswalk.
#   Regions: precision-weighted meta-regression (mixmeta, REML) of the cold-attributable rate,
#   the cold-attributable fraction and the baseline rate on regional poverty; and a two-factor
#   Shapley split of the high- vs low-poverty tercile gap into frailty and temperature risk.
# Output: tables/05_group_differences/decomposition.csv, decomposition_regions.csv,
#         regional_poverty_metaregression.csv, regional_tercile_decomposition.csv
# =============================================================================
source("2. Scripts/00_setup/config.R"); source("2. Scripts/00_setup/dlnm_tools.R")
STAGE <- "05_group_differences"; ND <- 200
panel <- load_panel(); setorder(panel, Region, Date)
PL <- read_model("pooled_Tmin"); mmt_all <- PL[["All|cum"]]$mmt; regs <- sort(unique(panel$Region))

# ---- person-years (two crosswalks) -----------------------------------------------------------------------------
py <- fread(file.path(DERIVED, "regional_populations.csv"))[Year %in% YEARS, .(PY = sum(Population)), by = .(Region, Group)]
pop <- fread(file.path(DERIVED, "county_joint_populations.csv"))[Year %in% YEARS]
pg <- rbind(pop[, .(Group = "All", P = sum(Population)), by = FIPS], pop[, .(P = sum(Population)), by = .(FIPS, Group = Sex)],
            pop[Race %in% c("Black", "White"), .(P = sum(Population)), by = .(FIPS, Group = Race)], pop[, .(P = sum(Population)), by = .(FIPS, Group = Age)])
cwn <- fread(file.path(DERIVED, "county_region_crosswalk.csv"))[, .(FIPS, Region = Nearest_any_station)]
py_near <- merge(pg, cwn, by = "FIPS")[, .(PY = sum(P)), by = .(Region, Group)]

# ---- attributable fractions by group x region (point and draws) ----------------------------------------------------------
rb <- function(x, at) onebasis(at, fun = "ns", knots = quantile(x, c(.1, .75, .9), na.rm = TRUE), Boundary.knots = range(x, na.rm = TRUE))
BX <- lapply(setNames(regs, regs), function(r) { x <- panel[Region == r]$Tmin_f; cen <- quantile(x, mmt_all)
  list(X = sweep(rb(x, x), 2, rb(x, cen), "-"), cold = x < cen, heat = x > cen) })
af_region <- function(g, r, cfc, cfh) { y <- panel[Region == r][[g]]; D <- sum(y, na.rm = TRUE)
  ac <- (1 - exp(-drop(BX[[r]]$X %*% cfc))) * frollmean(y, MAX_LAG + 1, align = "left", na.rm = TRUE)   # forward: deaths over the next 0-21 days
  ah <- (1 - exp(-drop(BX[[r]]$X %*% cfh))) * frollmean(y, 3, align = "left", na.rm = TRUE)
  c(cold = sum(ac[BX[[r]]$cold], na.rm = TRUE) / D, heat = sum(ah[BX[[r]]$heat], na.rm = TRUE) / D, D = D) }
AFp <- list(); AFd <- list()
for (g in GROUPS) { pc <- PL[[paste(g, "cum", sep = "|")]]; ph <- PL[[paste(g, "acute", sep = "|")]]
  AFp[[g]] <- rbindlist(lapply(regs, function(r) as.data.table(as.list(af_region(g, r, pc$blup[[r]]$blup, ph$blup[[r]]$blup)))[, Region := r]))[, Group := g]
  AFd[[g]] <- rbindlist(mclapply(seq_len(ND), function(s) { set.seed(SEED + s)
    rbindlist(lapply(regs, function(r) as.data.table(as.list(af_region(g, r, mvrnorm(1, pc$blup[[r]]$blup, pc$blup[[r]]$vcov),
                                                                         mvrnorm(1, ph$blup[[r]]$blup, ph$blup[[r]]$vcov))))[, `:=`(Region = r, draw = s)])) }, mc.cores = NCORES))[, Group := g]
  cat(g, "done\n") }
AFp <- rbindlist(AFp); AFd <- rbindlist(AFd)

# ---- Shapley decomposition of group gaps ------------------------------------------------------------------------------------
rate <- function(w, m, a) 1e5 * sum(w * m * a)
shapley3 <- function(A, B) {                     # A, B: one row per region (same order) with w, m, a
  f <- function(s) rate(if (s[1]) A$w else B$w, if (s[2]) A$m else B$m, if (s[3]) A$a else B$a)
  ct <- c(0, 0, 0)
  for (o in list(1:3, c(1, 3, 2), c(2, 1, 3), c(2, 3, 1), c(3, 1, 2), c(3, 2, 1))) { s <- c(FALSE, FALSE, FALSE); prev <- f(s)
    for (k in o) { s[k] <- TRUE; cur <- f(s); ct[k] <- ct[k] + (cur - prev) / 6; prev <- cur } }
  c(rate_A = f(c(TRUE, TRUE, TRUE)), rate_B = f(c(FALSE, FALSE, FALSE)), gap = f(c(TRUE, TRUE, TRUE)) - f(c(FALSE, FALSE, FALSE)),
    geography = ct[1], frailty = ct[2], sensitivity = ct[3]) }
decompose <- function(P_, lab) {
  base <- merge(AFp[, .(Region, Group, D)], P_, by = c("Region", "Group"))[, `:=`(w = PY / sum(PY), m = D / PY), by = Group]
  comp <- function(g, src, side) merge(base[Group == g, .(Region, w, m)], src[Group == g, .(Region, a = get(side))], by = "Region")[order(Region)]
  rbindlist(lapply(c("cold", "heat"), function(side) rbindlist(lapply(PAIRS, function(pr) {
    pt <- shapley3(comp(pr[1], AFp, side), comp(pr[2], AFp, side))
    dr <- t(sapply(seq_len(ND), function(s) shapley3(comp(pr[1], AFd[draw == s], side), comp(pr[2], AFd[draw == s], side))))
    data.table(crosswalk = lab, side = side, pair = paste(pr, collapse = "/"), component = names(pt), estimate = pt,
               lo = apply(dr, 2, quantile, .025), hi = apply(dr, 2, quantile, .975)) })))) }
DEC <- rbind(decompose(py, "rebuilt"), decompose(py_near, "nearest-station"))
write_tab(DEC, STAGE, "decomposition")

# ---- regions: poverty, frailty and temperature risk ---------------------------------------------------------------------------
cty <- fread(file.path(DATA, "reference", "county_baseline_analysis_data.csv"))[FIPS == 51515, FIPS := 51019]   # Bedford city -> Bedford County
cw <- fread(file.path(DERIVED, "county_region_crosswalk.csv"))
pov <- merge(cty, cw[, .(FIPS, Region)], by = "FIPS")[, .(poverty = weighted.mean(poverty_pct, Population)), by = Region]
RG <- merge(merge(AFp[Group == "All"], py[Group == "All", .(Region, PY)], by = "Region"), pov, by = "Region")[, m := D / PY]
rd <- merge(AFd[Group == "All", .(Region, draw, cold, heat)], RG[, .(Region, m)], by = "Region")
se <- rd[, .(se_cold_rate = sd(1e5 * m * cold), se_cold_af = sd(100 * cold), se_heat_rate = sd(1e5 * m * heat)), by = Region]
RG <- merge(RG, se, by = "Region")[, `:=`(baseline_rate = 1e5 * m, cold_rate = 1e5 * m * cold, cold_af = 100 * cold, heat_rate = 1e5 * m * heat,
                                         se_log_baseline = 1 / sqrt(D), tercile = cut(frank(poverty), c(0, 4, 8, 12), labels = c("low", "middle", "high")))]
write_tab(RG[, .(Region, poverty, tercile, PY, deaths = D, baseline_rate, cold_af, se_cold_af, cold_rate, se_cold_rate, heat_rate, se_heat_rate)], STAGE, "decomposition_regions")
mr <- function(y, s, lab, scale = 1) { f <- mixmeta(RG[[y]] ~ I(RG$poverty / 10), S = RG[[s]]^2, method = "reml"); cf <- coef(f)[2]; v <- vcov(f)[2, 2]
  data.table(outcome = lab, slope_per_10pp_poverty = scale * cf, lo = scale * (cf - 1.96 * sqrt(v)), hi = scale * (cf + 1.96 * sqrt(v)),
             p = 2 * pnorm(-abs(cf / sqrt(v))), I2 = summary(f)$i2stat[1], spearman = cor(RG$poverty, RG[[y]], method = "spearman")) }
RG[, log_baseline := log(baseline_rate)]
MR <- rbind(mr("cold_rate", "se_cold_rate", "cold-attributable deaths per 100,000"), mr("cold_af", "se_cold_af", "cold-attributable fraction (%)"),
            mr("heat_rate", "se_heat_rate", "heat-attributable deaths per 100,000 (lag 0-2)"),
            mr("log_baseline", "se_log_baseline", "baseline death rate (log; slope = ratio - 1)")[, `:=`(slope_per_10pp_poverty = exp(slope_per_10pp_poverty) - 1, lo = exp(lo) - 1, hi = exp(hi) - 1)])
write_tab(MR, STAGE, "regional_poverty_metaregression")
# high vs low poverty tercile: rate = M (baseline) x AF; 2-factor Shapley split, with draws
terc <- function(afcol) { d <- merge(RG[, .(Region, tercile, PY, m)], afcol, by = "Region")[tercile != "middle"]
  t <- d[, .(M = sum(PY * m) / sum(PY), AF = sum(PY * m * a) / sum(PY * m)), by = tercile]; H <- t[tercile == "high"]; L <- t[tercile == "low"]
  c(rate_high = 1e5 * H$M * H$AF, rate_low = 1e5 * L$M * L$AF, gap = 1e5 * (H$M * H$AF - L$M * L$AF),
    frailty = 1e5 * (H$M - L$M) * (H$AF + L$AF) / 2, temperature_risk = 1e5 * (H$AF - L$AF) * (H$M + L$M) / 2) }
TD <- rbindlist(lapply(c("cold", "heat"), function(side) {
  pt <- terc(RG[, .(Region, a = get(side))]); dr <- t(sapply(seq_len(ND), function(s) terc(AFd[Group == "All" & draw == s, .(Region, a = get(side))])))
  data.table(side = side, component = names(pt), estimate = pt, lo = apply(dr, 2, quantile, .025), hi = apply(dr, 2, quantile, .975)) }))
write_tab(TD, STAGE, "regional_tercile_decomposition")
options(width = 200)
print(DEC[crosswalk == "rebuilt" & component %in% c("gap", "geography", "frailty", "sensitivity"), .(side, pair, component, est = round(estimate, 2), lo = round(lo, 2), hi = round(hi, 2))], nrows = 100)
print(MR[, lapply(.SD, function(v) if (is.numeric(v)) signif(v, 3) else v)]); print(TD[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 2) else v)])
