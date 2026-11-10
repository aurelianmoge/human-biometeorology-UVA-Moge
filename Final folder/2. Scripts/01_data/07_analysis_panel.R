# =============================================================================
# 07_analysis_panel.R - the region-day panel used by every model.
#   * deaths by margin as extracted (never imputed; ambiguous region-days stay NA);
#   * suffix _f = gap-filled weather: Tmin_f, Tmax_f from script 05; Tmean_f = (Tmin_f + Tmax_f) / 2;
#     DTR_f = Tmax_f - Tmin_f; afternoon temperature, 7 am RH and dew point filled from
#     the other stations (month-specific calibration) for sensitivity models only;
#   * calendar terms: year, month, day of week, running day index.
# Input : 1. Data/derived/regional_panel.csv, weather_filled.csv
# Output: 1. Data/derived/analysis_panel.rds; 3. Outputs/tables/01_data/panel_audit.csv
# =============================================================================
source("2. Scripts/00_setup/config.R")
STAGE <- "01_data"
panel <- fread(file.path(DERIVED, "regional_panel.csv"))[, Date := as.IDate(Date)]
wf <- fread(file.path(DERIVED, "weather_filled.csv"))[, Date := as.IDate(Date)]
panel <- merge(panel, wf[, .(Region, Date, Tmin_f = Tmin, Tmax_f = Tmax, Tmin_source)], by = c("Region", "Date"))
panel[, `:=`(Tmean_f = (Tmin_f + Tmax_f) / 2, DTR_f = Tmax_f - Tmin_f)]

fill_from_stations <- function(dt, v) {                  # month-specific calibration on the other stations' mean
  w <- dcast(dt, Date ~ Region, value.var = v); m <- as.matrix(w[, -1]); mon <- factor(month(w$Date)); out <- m
  for (r in colnames(m)) { if (all(is.na(m[, r]))) next
    oth <- rowMeans(m[, setdiff(colnames(m), r), drop = FALSE], na.rm = TRUE); f <- lm(m[, r] ~ mon * oth)
    i <- is.na(m[, r]) & is.finite(oth); out[i, r] <- predict(f, data.frame(mon = mon, oth = oth))[i] }
  l <- melt(data.table(Date = w$Date, out), id.vars = "Date", variable.name = "Region", value.name = paste0(v, "_f"))
  l[, Region := as.character(Region)]; l
}
for (v in c("Tafternoon", "RH07", "Td07")) panel <- merge(panel, fill_from_stations(panel, v), by = c("Region", "Date"))
setorder(panel, Region, Date)
panel[, `:=`(year = year(Date), month = month(Date), dow = factor(wday(Date)), time = as.numeric(Date - min(Date)) + 1)]
saveRDS(panel, file.path(DERIVED, "analysis_panel.rds"))

audit <- panel[, .(days = .N, death_days_observed = sum(!is.na(All)), deaths = sum(All, na.rm = TRUE),
                   Tmin_missing = sum(is.na(Tmin)), Tmin_filled = sum(Tmin_source != "observed")), by = Region]
write_tab(audit, STAGE, "panel_audit")
cat("Deaths:", sum(audit$deaths), "| Tmin filled:", sum(audit$Tmin_filled), "of", sum(audit$days), "region-days\n")
