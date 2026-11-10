# =============================================================================
# 02_exposure_choice.R - which temperature metric fits daily mortality best?
# For each region and margin, the stage-1 DLNM (ns P10/75/90 x lag 0-21, 7 df/yr, dow)
# is fitted with Tmin, Tmean, Tmax and the 1 pm temperature on the SAME days, and
# compared by quasi-AIC (Peng et al. 2006; Gasparrini's fqaic): QAIC = -2 logL / phi + 2k.
# Lower is better; differences are summed over the 12 regions.
# DeGuzman et al. (2026) selected Tmin after a similar search; Davis et al. (2023) used
# a Tmax consensus model.  Input: analysis_panel.rds.  Output: tables/02_descriptives/
# =============================================================================
source("2. Scripts/00_setup/config.R"); source("2. Scripts/00_setup/dlnm_tools.R")
STAGE <- "02_descriptives"
panel <- load_panel()
EXPO <- c(Tmin = "Tmin_f", Tmean = "Tmean_f", Tmax = "Tmax_f", Afternoon = "Tafternoon_f")
qaic <- function(d, expo, y) {
  x <- d[[expo]]; cb <- crossbasis(x, lag = MAX_LAG, argvar = argvar_spec(x, "ns3"), arglag = list(knots = logknots(MAX_LAG, 3)))
  m <- glm(d[[y]] ~ cb + ns(time, df = DF_YEAR * 16) + dow, data = d, family = quasipoisson)
  ll <- sum(dpois(m$y, m$fitted.values, log = TRUE)); phi <- summary(m)$dispersion
  c(qaic = -2 * ll / phi + 2 * length(coef(m)), dispersion = phi, n = length(m$y))
}
jobs <- CJ(Region = sort(unique(panel$Region)), Group = GROUPS)
res <- rbindlist(mclapply(seq_len(nrow(jobs)), function(i) {
  d <- panel[Region == jobs$Region[i]][order(Date)]
  ok <- complete.cases(d[, EXPO, with = FALSE]); d[!ok, (EXPO) := NA]          # identical days for every metric
  rbindlist(lapply(names(EXPO), function(e) as.data.table(as.list(qaic(d, EXPO[[e]], jobs$Group[i])))[, `:=`(Region = jobs$Region[i], Group = jobs$Group[i], exposure = e)]))
}, mc.cores = NCORES))
res[, dQAIC := qaic - min(qaic), by = .(Region, Group)]
write_tab(res, STAGE, "T10_exposure_qaic_by_region")
summ <- res[, .(dQAIC_sum = sum(dQAIC), regions_best = sum(dQAIC == 0)), by = .(Group, exposure)]
write_tab(summ[order(Group, dQAIC_sum)], STAGE, "T10_exposure_qaic_summary")
print(dcast(summ, Group ~ exposure, value.var = "dQAIC_sum")[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 1) else x)])
print(dcast(summ, Group ~ exposure, value.var = "regions_best"))
