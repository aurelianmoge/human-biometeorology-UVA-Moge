# =============================================================================
# 01_two_stage_dlnm.R - main temperature-mortality model, by demographic margin.
# Stage 1: region x margin quasi-Poisson DLNM on the region's own station Tmin
#          (ns knots at regional P10/P75/P90; lag 0-21, 3 log-spaced knots;
#          ns(time, 7 df/yr); day of week), reduced to cumulative (0-21),
#          acute (0-2) and delayed (3-21) curves.  Tmean as secondary exposure.
# Stage 2: multivariate REML meta-analysis across the 12 regions (percentile scale);
#          MMT on P1-P99 with a 1,000-draw parametric bootstrap (Tobias et al. 2017).
# Output : 3. Outputs/models/stage1_<expo>.rds, pooled_<expo>.rds;
#          3. Outputs/tables/03_temperature/pooled_RR.csv
# =============================================================================
source("2. Scripts/00_setup/config.R"); source("2. Scripts/00_setup/dlnm_tools.R")
STAGE <- "03_temperature"; set.seed(SEED)
panel <- load_panel()
rows <- list()
for (e in c(Tmin = "Tmin_f", Tmean = "Tmean_f")) {
  en <- sub("_f", "", e)
  s1 <- run_stage1(panel, e)
  save_model(s1, paste0("stage1_", en))
  pooled <- list()
  for (g in GROUPS) for (w in c("cum", "acute", "delayed")) pooled[[paste(g, w, sep = "|")]] <- pool(s1, g, w, boot_mmt = TRUE)
  save_model(pooled, paste0("pooled_", en))
  for (k in names(pooled)) { p <- pooled[[k]]; gw <- strsplit(k, "\\|")[[1]]
    for (cen in c("MMT", "median")) { ref <- if (cen == "MMT") p$mmt else .5
      for (q in c(.01, .025, .10, .90, .975, .99)) { l <- lrr(p, q, ref)
        rows[[length(rows) + 1]] <- data.table(exposure = en, group = gw[1], window = gw[2], centre = cen, pct = 100 * q, temp = p$tavg[pidx(q)],
          RR = exp(l[1]), lo = exp(l[1] - 1.96 * l[2]), hi = exp(l[1] + 1.96 * l[2]),
          MMT_pct = 100 * p$mmt, MMT_temp = p$tavg[pidx(p$mmt)], MMT_lo = 100 * quantile(p$mmt_boot, .025), MMT_hi = 100 * quantile(p$mmt_boot, .975), I2 = p$i2) } } }
  cat(en, "done\n")
}
res <- rbindlist(rows); write_tab(res, STAGE, "pooled_RR")
print(res[exposure == "Tmin" & window == "cum" & centre == "MMT" & pct %in% c(1, 99),
          .(group, pct, RR = round(RR, 3), lo = round(lo, 3), hi = round(hi, 3), MMT_pct, MMT_lo, MMT_hi, I2 = round(I2))])
