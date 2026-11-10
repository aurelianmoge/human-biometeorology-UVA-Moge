# =============================================================================
# 01_descriptives.R - baseline description of deaths, populations and weather.
#   T1 deaths, person-years and crude rates by margin (statewide; Census denominators)
#   T2 indirectly age-standardised mortality ratios (SMR) by sex and race
#   T3 age structure of the Black and White populations (Census)
#   T4 regional deaths, populations, crude rates and SMR (crosswalk of script 01/06)
#   T5 regional climate: Tmin / Tmax / Tmean / DTR distribution
#   T6 annual trends (deaths, rates, mean Tmin); T7 seasonality (deaths per day by month)
#   T8 raw association: deaths relative to the stratum mean (region x year x month), by quintile of the
#      21-day mean Tmin within the same stratum (model-free analogue of the time-stratified design)
#   T9 correlation between candidate exposures (within region)
# Input : 1. Data/derived/analysis_panel.rds, statewide_population_margins.csv,
#         county_joint_populations.csv, regional_populations.csv
# Output: 3. Outputs/tables/02_descriptives/, 3. Outputs/figures/02_descriptives/
# =============================================================================
source("2. Scripts/00_setup/config.R"); source("2. Scripts/00_setup/dlnm_tools.R")
STAGE <- "02_descriptives"
panel <- load_panel()
spop <- fread(file.path(DERIVED, "statewide_population_margins.csv"))[Year %in% YEARS]
cpop <- fread(file.path(DERIVED, "county_joint_populations.csv"))[Year %in% YEARS]
rpop <- fread(file.path(DERIVED, "regional_populations.csv"))[Year %in% YEARS]
ages <- c("Under50", "Age50_74", "Age75plus")
long <- melt(panel[, c("Region", "Date", "year", "month", GROUPS), with = FALSE], id.vars = c("Region", "Date", "year", "month"),
             variable.name = "Group", value.name = "deaths")
ann  <- long[, .(deaths = mean(deaths, na.rm = TRUE) * .N, obs_days = sum(!is.na(deaths))), by = .(Group, year)]   # region-days with missing deaths imputed at the year mean, for rates only

# ---- T1 statewide deaths and crude rates ----------------------------------------------
t1 <- merge(long[, .(deaths = sum(deaths, na.rm = TRUE), deaths_per_day = sum(deaths, na.rm = TRUE) / uniqueN(Date)), by = Group],
            spop[, .(person_years = sum(Population)), by = Group], by = "Group")
t1[, `:=`(share_of_deaths = deaths / t1[Group == "All", deaths], crude_rate_per100k = 1e5 * deaths / person_years)]
write_tab(t1[match(GROUPS, Group)], STAGE, "T1_deaths_rates_statewide")

# ---- T2 indirect age standardisation by sex and race ---------------------------------------
arate <- merge(ann[Group %in% ages, .(D = sum(deaths)), by = .(Group, year)], spop[Group %in% ages, .(Group, year = Year, P = Population)])[, rate := D / P]
py <- rbind(cpop[, .(P = sum(Population)), by = .(year = Year, Age, margin = Sex)], cpop[Race %in% c("Black", "White"), .(P = sum(Population)), by = .(year = Year, Age, margin = Race)])
ex <- merge(py, arate[, .(year, Age = Group, rate)], by = c("year", "Age"))[, .(expected = sum(P * rate)), by = margin]
t2 <- merge(ann[Group %in% c("Female", "Male", "Black", "White"), .(observed = sum(deaths)), by = .(margin = Group)], ex, by = "margin")[, SMR := observed / expected]
t2 <- rbind(t2, data.table(margin = c("Female/Male", "Black/White"), observed = NA, expected = NA,
                           SMR = c(t2[margin == "Female", SMR] / t2[margin == "Male", SMR], t2[margin == "Black", SMR] / t2[margin == "White", SMR])))
write_tab(t2, STAGE, "T2_indirect_SMR_sex_race")

# ---- T3 age structure of populations and of deaths --------------------------------------------
t3 <- cpop[Race %in% c("Black", "White"), .(P = sum(Population)), by = .(Race, Age)][, share := P / sum(P), by = Race]
t3 <- merge(t3, rbind(cpop[, .(P = sum(Population)), by = Age][, `:=`(Race = "All", share = P / sum(P))]), all = TRUE)
write_tab(dcast(t3, Race ~ Age, value.var = "share")[, c("Race", ages), with = FALSE], STAGE, "T3_population_age_structure")

# ---- T4 regions ---------------------------------------------------------------------------------
cwv <- fread(file.path(OUT, "tables", "01_data", "crosswalk_validation.csv"))
t4 <- merge(long[Group == "All", .(deaths = sum(deaths, na.rm = TRUE), deaths_per_day = mean(deaths, na.rm = TRUE)), by = Region],
            rpop[Group == "All", .(person_years = sum(Population)), by = Region], by = "Region")
t4 <- merge(t4, long[Group %in% c("Black", "Age75plus"), .(d = sum(deaths, na.rm = TRUE)), by = .(Region, Group)][, dcast(.SD, Region ~ Group, value.var = "d")], by = "Region")
t4[, `:=`(crude_rate_per100k = 1e5 * deaths / person_years, share_Black = Black / deaths, share_75plus = Age75plus / deaths)]
t4 <- merge(t4[, .(Region, deaths, deaths_per_day, person_years, crude_rate_per100k, share_Black, share_75plus)], cwv[, .(Region, SMR, localities, inferred)], by = "Region")
write_tab(t4[order(-deaths)], STAGE, "T4_regions")

# ---- T5 climate by region ------------------------------------------------------------------------
q <- function(x) as.list(setNames(quantile(x, c(.01, .1, .5, .9, .99), na.rm = TRUE), c("P1", "P10", "P50", "P90", "P99")))
t5 <- panel[, c(list(Tmin_mean = mean(Tmin_f), Tmin_sd = sd(Tmin_f)), q(Tmin_f), list(Tmax_mean = mean(Tmax_f), Tmean_mean = mean(Tmean_f),
        DTR_mean = mean(DTR_f), DTR_P10 = quantile(DTR_f, .1), DTR_P90 = quantile(DTR_f, .9), frost_days_share = mean(Tmin_f < 0), warm_nights_share = mean(Tmin_f >= 20))), by = Region]
write_tab(t5, STAGE, "T5_climate_by_region")

# ---- T6 / T7 trends and seasonality -----------------------------------------------------------------
t6 <- merge(ann, spop[, .(Group, year = Year, Population)], by = c("Group", "year"))[, crude_rate_per100k := 1e5 * deaths / Population]
t6 <- merge(t6, panel[, .(mean_Tmin = mean(Tmin_f), mean_DTR = mean(DTR_f)), by = year], by = "year")
write_tab(t6, STAGE, "T6_annual")
t7 <- long[, .(deaths_per_day = mean(deaths, na.rm = TRUE)), by = .(Group, Region, month)][, .(deaths_per_day = sum(deaths_per_day)), by = .(Group, month)]
t7[, relative_to_annual := deaths_per_day / mean(deaths_per_day), by = Group]
write_tab(t7, STAGE, "T7_seasonality")

# ---- T8 raw association (no model), time-stratified: within each region x year x month, days are ranked by the
#      mean Tmin of lags 0-21 (quintiles) and deaths are divided by the stratum mean -------------------------------
panel[, Tmin21 := frollmean(Tmin_f, 22, align = "right"), by = Region]
panel[!is.na(Tmin21), Tmin_q := as.integer(ceiling(5 * frank(Tmin21) / .N)), by = .(Region, year, month)]
t8 <- melt(panel[!is.na(Tmin_q), c("Region", "year", "month", "Tmin_q", GROUPS), with = FALSE], id.vars = c("Region", "year", "month", "Tmin_q"),
           variable.name = "Group", value.name = "d")[!is.na(d)]
t8[, rel := d / mean(d), by = .(Group, Region, year, month)]
t8 <- t8[is.finite(rel), .(relative_deaths = mean(rel), se = sd(rel) / sqrt(.N), n = .N), by = .(Group, Tmin_q)]
write_tab(t8, STAGE, "T8_relative_deaths_by_Tmin21_quintile")

# ---- T9 correlation of candidate exposures within region ---------------------------------------------
vv <- c("Tmin_f", "Tmax_f", "Tmean_f", "Tafternoon_f", "DTR_f", "RH07_f", "Td07_f")
cl <- lapply(split(panel, by = "Region"), function(d) suppressWarnings(cor(d[, ..vv], use = "pairwise.complete.obs")))
cm <- apply(simplify2array(cl), c(1, 2), mean, na.rm = TRUE)                    # mean over regions with data (VJI has no RH)
write_tab(as.data.table(round(cm, 3), keep.rownames = "variable"), STAGE, "T9_exposure_correlations")

# ---- figures (direct labels: colours repeat across the sex / race / age panels) --------------------------
fam <- c(All = "All", Female = "Sex", Male = "Sex", Black = "Race", White = "Race", Under50 = "Age", Age50_74 = "Age", Age75plus = "Age")
famf <- function(d) { d[, `:=`(fam = factor(fam[as.character(Group)], c("Sex", "Race", "Age")), lab = GROUP_LABEL[as.character(Group)])]; d }
endlab <- function(d, x) d[get(x) == max(get(x))]
lines_by_group <- function(d, x, y, title, subtitle, xlab = NULL, ylab = NULL, ref = NULL) {
  d <- famf(copy(d[Group != "All"])); e <- endlab(d, x)
  ggplot(d, aes(.data[[x]], .data[[y]], colour = Group)) + { if (!is.null(ref)) geom_hline(yintercept = ref, colour = COL[["grey"]]) } +
    geom_line(linewidth = .8) + geom_point(size = 1.1) + geom_text(data = e, aes(label = lab), hjust = -.15, size = 3.1, show.legend = FALSE) +
    facet_wrap(~fam) + scale_colour_manual(values = GROUP_COL, guide = "none") + scale_x_continuous(expand = expansion(mult = c(.03, .22))) +
    labs(title = title, subtitle = subtitle, x = xlab, y = ylab)
}
save_fig(lines_by_group(t6, "year", "crude_rate_per100k", "Crude death rate per 100,000 by year",
         "Statewide deaths (12 workbooks) / Census July population. 2020 includes COVID-19.") + facet_wrap(~fam, scales = "free_y"), STAGE, "F1_annual_rates", 10, 3.8)
save_fig(lines_by_group(t7, "month", "relative_to_annual", "Seasonality: deaths per day relative to the group's annual mean",
         "Winter excess is largest for 75+ and smallest for under-50s.", ref = 1) +
         scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c("Jan", "Apr", "Jul", "Oct"), expand = expansion(mult = c(.03, .22))), STAGE, "F2_seasonality", 10, 3.8)
g3 <- ggplot(panel, aes(Tmin_f, reorder(Region, Tmin_f, median))) + geom_boxplot(outlier.size = .2, colour = COL[["ink2"]], fill = COL[["lgrey"]], coef = 1.5) +
  labs(title = "Daily minimum temperature by region, 2005–2020", subtitle = "Stations ordered by median Tmin; gap-filled series.", x = "Tmin (°C)", y = NULL)
save_fig(g3, STAGE, "F3_Tmin_by_region", 8, 4.5)
save_fig(lines_by_group(t8, "Tmin_q", "relative_deaths", "Raw association: colder-than-usual weeks, same region, year and month",
         "Deaths / stratum mean (region × year × month), by quintile of the mean Tmin over lags 0–21 within the stratum. Model-free; 1 = coldest.",
         xlab = "Quintile of 21-day mean Tmin within the stratum", ylab = "Relative deaths", ref = 1) +
         scale_x_continuous(breaks = 1:5, expand = expansion(mult = c(.03, .25))), STAGE, "F4_raw_association_Tmin", 10, 3.8)
g5 <- ggplot(t4, aes(SMR, reorder(Region, SMR))) + geom_vline(xintercept = 1, colour = COL[["grey"]]) + geom_point(colour = COL[["blue"]], size = 2.6) +
  geom_text(aes(label = sprintf("%.0f", crude_rate_per100k)), hjust = -.4, size = 3, colour = COL[["ink2"]]) +
  labs(title = "Regional mortality: indirectly age-standardised ratio (crude rate per 100,000 as label)", x = "SMR (state = 1)", y = NULL)
save_fig(g5, STAGE, "F5_regional_SMR", 8, 4.3)
print(t1[match(GROUPS, Group), .(Group, deaths, share = round(share_of_deaths, 3), rate = round(crude_rate_per100k))]); print(t2)
