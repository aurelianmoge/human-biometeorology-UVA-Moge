# Allocate Census population to classified days under a shared weather proxy.
source("2. Scripts/00_shared/setup.R")
theme <- THEMES[4]
annual <- read_derived("statewide_population_margins")
flags <- read.csv(file.path(OUTPUT, theme, "csv/classification_date_flags_RR105.csv"))
flags$Date <- as.Date(flags$Date)
allocations <- list()
for (group in unique(flags$Group)) {
  population <- annual[annual$Group == group, ]
  knots <- as.numeric(as.Date(paste0(population$Year, "-07-01")))
  daily <- approx(knots, population$Population, xout = as.numeric(DATES), rule = 2)$y
  selected <- flags[flags$Group == group, ]
  selected$Population <- daily[match(selected$Date, DATES)]
  allocations[[group]] <- selected
}
allocations <- bind_rows(allocations)
summary <- allocations |> group_by(Exposure, Group, Window) |>
  summarise(Assessed_days = n(), Denominator_population_days = sum(Population),
            Pooled_population_days = sum(Population * Pooled), Own_population_days = sum(Population * Own),
            Own_only_population_days = sum(Population * Own * !Pooled),
            Pooled_only_population_days = sum(Population * Pooled * !Own), .groups = "drop") |>
  mutate(Change_population_days = Own_population_days - Pooled_population_days,
         Change_percentage_points = 100 * Change_population_days / Denominator_population_days,
         Status = "full-Virginia allocation scenario under a shared proxy; unverified catchments; not unique people or attributable deaths")
stopifnot(max(abs(summary$Change_population_days - summary$Own_only_population_days + summary$Pooled_only_population_days)) < .01)
write_table(summary, "population_day_allocation_scenarios", theme)
write_table(annual, "annual_population_descriptions", THEMES[2])
save_plot(ggplot(annual |> filter(Group %in% c("Female", "Male", "Black", "White", "Under50", "Age75plus")),
  aes(Year, Population / 1e6, color = Group)) + geom_line() +
  labs(title = "Corrected Census population descriptions", y = "Million residents", caption = "Recorded demographic margins overlap across race, age and sex; do not sum these lines."),
  "population_demographic_trajectories", THEMES[2])
