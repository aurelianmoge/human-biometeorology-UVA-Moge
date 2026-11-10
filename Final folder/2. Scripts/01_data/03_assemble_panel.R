# =============================================================================
# 03_assemble_panel.R - regional daily calendars from the extracted workbook cells.
# (Data-preparation logic inherited unchanged from the earlier audited pipeline.)
#   * complete 12 x 5,844 region-day calendar, 2005-2020;
#   * a duplicated outcome date makes the whole region-date missing (never summed);
#   * physically impossible / inverted Tmin-Tmax and unsupported humidity set to missing;
#   * check: identical to the frozen reference panel (1. Data/reference/regional_panel.csv).
# Input : 1. Data/derived/workbook_cells.csv (script 01)
# Output: 1. Data/derived/regional_panel.csv, aggregate_daily.csv, dtr_daily.csv;
#         3. Outputs/tables/01_data/ (missingness, quality flags, margin reconciliation)
# =============================================================================
source("2. Scripts/00_setup/config.R")
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })
STAGE <- "01_data"
DATES <- seq(as.Date("2005-01-01"), as.Date("2020-12-31"), by = "day")
CELL_GROUPS <- c("All", "Black", "White", "Female", "Male", "Under50", "Age50_74", "Age75plus", "Respiratory", "Circulatory")
EXPOSURES <- c("Tmin", "Tmax", "Tmean_midrange", "DTR")
write_table <- function(data, name, theme) write_tab(as.data.table(data), STAGE, name)
read_derived <- function(name) { d <- read.csv(file.path(DATA, "derived", paste0(name, ".csv")), check.names = FALSE)
  if ("Date" %in% names(d)) d$Date <- as.Date(d$Date); d }
theme <- STAGE
raw <- read_derived("workbook_cells")
weather <- c("Tmin", "Tmax", "T07", "Td07", "RH07", "Tafternoon", "Wind", "Pressure", "PM", "Ozone")

###################################
### Resolve the regional calendar ###
###################################

# Duplicate outcomes are ambiguous. Their entire region-date is missing.
# Weather duplicates are averaged; no duplicated deaths are silently summed.
regional <- raw |>
  group_by(Region, Date) |>
  summarise(
    Source_rows = n(),
    across(all_of(CELL_GROUPS), function(x) if (n() > 1) NA_real_ else x[1]),
    across(all_of(weather), function(x) if (all(!is.finite(x))) NA_real_ else mean(x[is.finite(x)])),
    .groups = "drop"
  ) |>
  group_by(Region) |>
  complete(Date = DATES) |>
  arrange(Date, .by_group = TRUE) |>
  ungroup()

######################################
### Weather quality and definitions ###
######################################

inverted <- with(regional, is.finite(Tmin) & is.finite(Tmax) & Tmax < Tmin)
impossible <- with(regional, (is.finite(Tmin) & (Tmin < -50 | Tmin > 45)) |
                           (is.finite(Tmax) & (Tmax < -40 | Tmax > 55)))
humidity_invalid <- with(regional, is.finite(RH07) & (RH07 < 0 | RH07 > 100))
dewpoint_conflict <- with(regional, is.finite(T07) & is.finite(Td07) & Td07 > T07 + 0.5)
calculated_rh <- with(regional, 100 * exp(17.625 * Td07 / (243.04 + Td07) - 17.625 * T07 / (243.04 + T07)))
humidity_unsupported <- humidity_invalid | dewpoint_conflict |
  !is.finite(calculated_rh) | (is.finite(regional$RH07) & abs(regional$RH07 - calculated_rh) > 15)
flags <- regional |> select(Region, Date) |>
  mutate(Inverted = inverted, Impossible = impossible, Humidity_unsupported = humidity_unsupported)
regional$Tmin[inverted | impossible] <- NA_real_
regional$Tmax[inverted | impossible] <- NA_real_
regional$RH07[humidity_unsupported] <- NA_real_
regional <- regional |>
  group_by(Region) |>
  mutate(DTR = Tmax - Tmin, Tmean_midrange = (Tmax + Tmin) / 2,
         DeltaTmean = Tmean_midrange - lag(Tmean_midrange)) |>
  ungroup()

# Complete station sets are fixed across dates; missing stations never change weights.
aggregate_panel <- function(data) {
  data |> group_by(Date) |>
    summarise(across(all_of(CELL_GROUPS), function(x) sum(x, na.rm = FALSE)),
              across(all_of(c(EXPOSURES, "RH07", "PM", "Ozone", "DeltaTmean")),
                     function(x) if (all(is.finite(x))) mean(x) else NA_real_), .groups = "drop") |>
    arrange(Date)
}
humidity_support <- regional |> group_by(Region) |>
  summarise(Observed = sum(is.finite(RH07)), .groups = "drop")
dtr_regions <- humidity_support$Region[humidity_support$Observed > 1000]
aggregate <- aggregate_panel(regional)
dtr <- aggregate_panel(filter(regional, Region %in% dtr_regions))
for (entry in list(c("regional_panel", "regional"), c("aggregate_daily", "aggregate"), c("dtr_daily", "dtr"))) {
  write.csv(get(entry[2]), file.path(DATA, "derived", paste0(entry[1], ".csv")), row.names = FALSE, na = "")
}
write_table(flags, "weather_quality_flags", theme)
write_table(data.frame(Region = dtr_regions, Rule = "more than 1000 physically supported RH07 days"), "dtr_fixed_footprint", theme)
write_table(raw |> count(Region, Date) |> filter(n > 1), "ambiguous_outcome_dates", theme)
write_table(raw |> mutate(Age_residual = All - Under50 - Age50_74 - Age75plus,
                         Sex_residual = All - Female - Male,
                         Race_residual = All - Black - White) |>
              group_by(Region) |> summarise(Age_failures = sum(Age_residual != 0, na.rm = TRUE),
              Sex_unallocated = sum(Sex_residual, na.rm = TRUE),
              Other_or_unknown_race = sum(Race_residual, na.rm = TRUE), .groups = "drop"), "margin_reconciliation", theme)

long <- regional |> select(Region, Date, all_of(c(EXPOSURES, "RH07", CELL_GROUPS))) |>
  pivot_longer(-c(Region, Date), names_to = "Variable", values_to = "Value")
coverage <- long |> group_by(Region, Variable) |>
  summarise(Calendar_days = n(), Observed = sum(is.finite(Value)), Missing = sum(!is.finite(Value)),
            Missing_percent = 100 * mean(!is.finite(Value)), .groups = "drop")
write_table(coverage, "missingness", theme)
write_table(long |> mutate(Year = as.integer(format(Date, "%Y"))) |> group_by(Region, Variable, Year) |>
              summarise(Observed = sum(is.finite(Value)), Missing = sum(!is.finite(Value)), .groups = "drop"), "annual_missingness", theme)

# Compare independent assembly with the previous corrected panel, rather than
# accepting an old intermediate as a new raw observation.
frozen <- file.path(DATA, "reference", "regional_panel.csv")
if (file.exists(frozen)) {
  reference <- read.csv(frozen)
  reference$Date <- as.Date(reference$Date)
  compared <- inner_join(regional, reference, by = c("Region", "Date"), suffix = c("_new", "_old"))
  checks <- lapply(c(EXPOSURES, CELL_GROUPS, "RH07"), function(variable) {
    new <- compared[[paste0(variable, "_new")]]
    old <- compared[[paste0(variable, "_old")]]
    both <- is.finite(new) & is.finite(old)
    data.frame(Variable = variable, Rows = nrow(compared), Missingness_disagreements = sum(is.finite(new) != is.finite(old)),
               Maximum_absolute_difference = if (any(both)) max(abs(new[both] - old[both])) else NA_real_)
  }) |> bind_rows()
  write_table(checks, "independent_panel_reproduction", theme)
  stopifnot(all(checks$Missingness_disagreements == 0), all(checks$Maximum_absolute_difference < 1e-7, na.rm = TRUE))
}

