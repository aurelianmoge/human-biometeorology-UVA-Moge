# =============================================================================
# 03_exposure_by_group.R - the weather each group experiences, given where it lives.
#   * regional daily station series, weighted by the group's share of person-years in each
#     region (Census county populations, rebuilt county-region crosswalk);
#   * per person-year: mean Tmin, mean DTR, frost days (Tmin < 0 C), very cold days
#     (Tmin <= -10 C), warm nights (Tmin >= 22 C), days with DTR >= 15 C, days with heat
#     index >= 100 F (NWS Rothfusz heat index from Tmax and the 7 am dew point);
#   * sensitivity: weights from the nearest-station crosswalk, and from where each group's
#     deaths occur (needs no crosswalk).
# Output: tables/02_descriptives/exposure_by_group.csv, exposure_by_region.csv,
#         residential_distribution.csv
# =============================================================================
source("2. Scripts/00_setup/config.R"); source("2. Scripts/00_setup/dlnm_tools.R")
STAGE <- "02_descriptives"
panel <- load_panel()

# ---- person-years by region and group, under the two crosswalks --------------------------------------------
py <- fread(file.path(DERIVED, "regional_populations.csv"))[Year %in% YEARS, .(W = sum(Population)), by = .(Region, Group)]
pop <- fread(file.path(DERIVED, "county_joint_populations.csv"))[Year %in% YEARS]
pg <- rbind(pop[, .(Group = "All", P = sum(Population)), by = FIPS], pop[, .(P = sum(Population)), by = .(FIPS, Group = Sex)],
            pop[Race %in% c("Black", "White"), .(P = sum(Population)), by = .(FIPS, Group = Race)], pop[, .(P = sum(Population)), by = .(FIPS, Group = Age)])
cwn <- fread(file.path(DERIVED, "county_region_crosswalk.csv"))[, .(FIPS, Region = Nearest_any_station)]
py_near <- merge(pg, cwn, by = "FIPS")[, .(W = sum(P)), by = .(Region, Group)]
dw <- melt(panel[, lapply(.SD, sum, na.rm = TRUE), by = Region, .SDcols = GROUPS], id.vars = "Region", variable.name = "Group", value.name = "W")

# ---- regional exposure indicators ------------------------------------------------------------------------------
rh <- function(t, td) pmin(100, 100 * exp(17.625 * td / (243.04 + td) - 17.625 * t / (243.04 + t)))
heat_index_f <- function(tf, r) { s <- 0.5 * (tf + 61 + (tf - 68) * 1.2 + r * 0.094)
  ro <- -42.379 + 2.04901523 * tf + 10.14333127 * r - .22475541 * tf * r - .00683783 * tf^2 - .05481717 * r^2 +
    .00122874 * tf^2 * r + .00085282 * tf * r^2 - .00000199 * tf^2 * r^2
  ifelse((s + tf) / 2 >= 80, ro, s) }
w <- copy(panel[, .(Region, Tmin_f, Tmax_f, DTR_f, T07, Td07_f)])
w[Region == "VJI", Td07_f := (Td07_f - 32) * 5 / 9]                    # VJI dew point is stored in deg F (data issue, see report)
w[, Td := pmin(fifelse(is.na(Td07_f), Tmin_f, Td07_f), fifelse(is.na(T07), Tmax_f, T07), Tmax_f)]   # dew point <= temperature
w[, HI := heat_index_f(Tmax_f * 9 / 5 + 32, rh(Tmax_f, Td))]
EXV <- c("mean_Tmin", "mean_DTR", "frost_days", "very_cold_days", "warm_nights", "high_DTR_days", "heat_index_100F_days")
ry <- w[, .(mean_Tmin = mean(Tmin_f), mean_DTR = mean(DTR_f), frost_days = 365.25 * mean(Tmin_f < 0), very_cold_days = 365.25 * mean(Tmin_f <= -10),
            warm_nights = 365.25 * mean(Tmin_f >= 22), high_DTR_days = 365.25 * mean(DTR_f >= 15), heat_index_100F_days = 365.25 * mean(HI >= 100, na.rm = TRUE)), by = Region]
write_tab(ry, STAGE, "exposure_by_region")

# ---- group exposure = regional indicators weighted by where the group lives --------------------------------------
expo <- function(wt, lab) { m <- merge(wt, ry, by = "Region")[, s := W / sum(W), by = Group]
  m[, c(list(weights = lab), lapply(.SD, function(v) sum(s * v))), by = Group, .SDcols = EXV] }
E <- rbind(expo(py, "residents, rebuilt crosswalk"), expo(py_near, "residents, nearest-station crosswalk"), expo(dw, "deaths"))
write_tab(E, STAGE, "exposure_by_group")
geo <- dcast(py[Group != "All"][, share := W / sum(W), by = Group], Region ~ Group, value.var = "share")
write_tab(geo, STAGE, "residential_distribution")
print(E[weights == "residents, rebuilt crosswalk", lapply(.SD, function(v) if (is.numeric(v)) round(v, 2) else v)])
