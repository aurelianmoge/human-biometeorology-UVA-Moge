# Create an interpretation-oriented atlas, including null and uncertain results.
source("2. Scripts/00_shared/setup.R")
theme <- THEMES[4]
curves <- read.csv(file.path(OUTPUT, THEMES[3], "csv/response_curves.csv"))
estimates <- read.csv(file.path(OUTPUT, THEMES[3], "csv/response_estimates.csv"))

for (exposure in EXPOSURES) {
  selected <- curves |> filter(Exposure == exposure, Specification == "primary",
    Outcome %in% c("All", "Female", "Male", "Black", "White", "Under50", "Age75plus"),
    Window %in% c("Acute_0_2", "Cumulative_0_21"))
  save_plot(ggplot(selected, aes(Exposure_value, RR, color = Outcome)) +
    geom_hline(yintercept = 1, linetype = 2) + geom_line() + facet_wrap(~Window, scales = "free_y") +
    labs(title = paste(exposure, "responses by demographic margin"), x = "Degrees Celsius", y = "RR vs common metric median",
      caption = "Differences between separate curves require a direct contrast; no age standardization of race/sex effects."),
    paste0("subgroup_curves_", exposure), theme, 12, 6)
}

# Sensitivity panels expose the dependence on calendar adjustment and count family.
selected <- estimates |> filter(Outcome == "All", Window == "Cumulative_0_21")
save_plot(ggplot(selected, aes(RR, Specification, color = factor(round(Exposure_value, 1)))) +
  geom_vline(xintercept = 1, linetype = 2) + geom_point(position = position_dodge(.5)) +
  geom_errorbar(aes(xmin = Lower, xmax = Upper), width = .2, position = position_dodge(.5), orientation = "y") +
  facet_wrap(~Exposure, scales = "free") + labs(title = "Model specification sensitivity", x = "Cumulative RR", y = NULL, color = "Exposure °C",
    caption = "References and weather basis remain fixed. Omitting 2020 changes the fitted calendar; nominal intervals are diagnostics."),
  "specification_sensitivity_forest", theme, 14, 10)

causes <- estimates |> filter(Outcome %in% c("Respiratory", "Circulatory"), Specification == "primary", Window == "Cumulative_0_21")
write_table(causes, "cause_specific_response_summary", theme)
save_plot(ggplot(causes, aes(RR, Outcome, color = factor(round(Exposure_value, 1)))) +
  geom_vline(xintercept = 1, linetype = 2) + geom_point(position = position_dodge(.5)) +
  geom_errorbar(aes(xmin = Lower, xmax = Upper), width = .2, position = position_dodge(.5), orientation = "y") +
  facet_wrap(~Exposure, scales = "free_x") + labs(title = "Cause-specific cumulative responses", x = "RR", y = NULL, color = "Exposure °C",
    caption = "Cause margins are supplied workbook classifications. Separate intervals do not test a difference between causes."),
  "cause_specific_forest", theme, 13, 7)

lag <- read.csv(file.path(OUTPUT, THEMES[3], "csv/lag_profiles.csv")) |>
  group_by(Exposure) |> mutate(Anchor = ifelse(Exposure_value == min(Exposure_value), "Lower anchor", "Upper anchor")) |> ungroup()
save_plot(ggplot(lag, aes(Lag, Exposure, fill = Log_RR)) + geom_tile() + facet_wrap(~Anchor) +
  scale_fill_gradient2(low = "#236B8E", mid = "white", high = "#AF4528", name = "log RR") +
  labs(title = "Exposure-specific lag timing", y = NULL, caption = "P10/P90 for temperature; P95/P99 for DTR. Tiles show effect magnitude, not significance."),
  "lag_effect_heatmap", theme, 12, 6)
residual <- read.csv(file.path(OUTPUT, theme, "csv/residual_calendar_autocorrelation.csv"))
save_plot(ggplot(residual, aes(Lag, ACF, color = Exposure)) + geom_hline(yintercept = 0, linetype = 2) + geom_line() +
  labs(title = "Residual dependence after temporal adjustment", y = "Calendar-pair residual correlation", caption = "Missing residuals are not compressed into consecutive days; no iid confidence band is asserted."),
  "residual_autocorrelation", theme)
