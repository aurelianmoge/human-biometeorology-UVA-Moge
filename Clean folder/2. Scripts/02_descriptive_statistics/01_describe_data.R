# Descriptive associations use observed records, with no causal interpretation.
source("2. Scripts/00_shared/setup.R")
theme <- THEMES[2]
regional <- read_derived("regional_panel")
aggregate <- read_derived("aggregate_daily")

############################
### Weather and mortality ###
############################

long <- regional |> select(Region, Date, all_of(c(EXPOSURES, GROUPS))) |>
  pivot_longer(-c(Region, Date), names_to = "Variable", values_to = "Value")
summary <- long |> group_by(Region, Variable) |>
  summarise(N = sum(is.finite(Value)), Mean = mean(Value, na.rm = TRUE), SD = sd(Value, na.rm = TRUE),
            P01 = quantile(Value, .01, na.rm = TRUE), P10 = quantile(Value, .1, na.rm = TRUE),
            Median = median(Value, na.rm = TRUE), P90 = quantile(Value, .9, na.rm = TRUE),
            P99 = quantile(Value, .99, na.rm = TRUE), .groups = "drop")
write_table(summary, "regional_summary", theme)
monthly <- long |> mutate(Month = as.integer(format(Date, "%m"))) |>
  group_by(Region, Variable, Month) |>
  summarise(Observed = sum(is.finite(Value)), Mean = mean(Value, na.rm = TRUE), .groups = "drop")
write_table(monthly, "seasonal_summary", theme)
annual <- regional |> mutate(Year = as.integer(format(Date, "%Y"))) |> group_by(Region, Year) |>
  summarise(across(all_of(GROUPS), function(x) sum(x, na.rm = TRUE)),
            Observed_outcome_days = sum(is.finite(All)), .groups = "drop")
write_table(annual, "annual_mortality_observed_days", theme)

save_plot(ggplot(monthly |> filter(Variable %in% EXPOSURES), aes(Month, Mean, color = Region)) +
            geom_line() + facet_wrap(~Variable, scales = "free_y") +
            labs(title = "Seasonal weather profiles by workbook region", y = "Degrees Celsius", x = "Month"),
          "seasonal_weather", theme, 12, 8)
save_plot(ggplot(annual, aes(Year, All, color = Region)) + geom_line() +
            labs(title = "Observed annual deaths by region", y = "Death count", caption = "Incomplete regional calendars affect annual totals; these are not mortality rates."),
          "annual_deaths", theme)

cloud <- aggregate |> select(Date, All, all_of(EXPOSURES)) |> pivot_longer(all_of(EXPOSURES), names_to = "Exposure", values_to = "Value")
write_table(cloud, "mortality_weather_cloud_data", theme)
save_plot(ggplot(cloud, aes(Value, All)) + geom_point(alpha = .12, size = .6, color = "#236B8E") +
            geom_smooth(method = "loess", se = FALSE, color = "#AF4528") + facet_wrap(~Exposure, scales = "free_x") +
            labs(title = "Daily aggregate mortality and station weather", x = "Degrees Celsius", y = "Supplied deaths", caption = "Raw association: season, trends, population composition and lag histories are not controlled."),
          "mortality_temperature_clouds", theme, 12, 9)
correlation <- cor(regional[c(EXPOSURES, "All", "RH07", "Wind", "Pressure", "PM", "Ozone")], use = "pairwise.complete.obs")
correlation <- as.data.frame(as.table(correlation))
names(correlation) <- c("First", "Second", "Correlation")
write_table(correlation, "weather_mortality_correlations", theme)
save_plot(ggplot(correlation, aes(First, Second, fill = Correlation)) + geom_tile() +
            scale_fill_gradient2(low = "#236B8E", mid = "white", high = "#AF4528", limits = c(-1, 1)) +
            labs(title = "Pooled regional correlations", x = NULL, y = NULL,
                 caption = "Pairwise observed rows; unequal missingness and between-region differences affect comparisons.") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)), "correlation_heatmap", theme, 10, 8)

########################################
### County socioeconomic descriptions ###
########################################

county <- read.csv(file.path(DATA, "reference", "county_baseline_analysis_data.csv"))
write_table(county, "county_baseline_panel", theme)
county_variables <- c("Crude_rate", "Age_adjusted_rate", "poverty_pct", "median_household_income", "Black_share_pct", "Age75plus_pct", "Population")
write_table(county |> summarise(across(all_of(county_variables), list(Mean = ~mean(.x, na.rm = TRUE), SD = ~sd(.x, na.rm = TRUE)))), "county_summary", theme)
save_plot(ggplot(county, aes(poverty_pct, Age_adjusted_rate, color = factor(Year))) +
            geom_point(alpha = .6) + geom_smooth(method = "lm", se = FALSE) +
            labs(title = "County poverty and annual age-adjusted mortality", x = "Poverty (%)", y = "Published age-adjusted rate per 100,000", color = "Year",
                 caption = "Ecological annual comparison; it does not estimate modification of acute temperature response."),
          "county_poverty_mortality_cloud", theme)
