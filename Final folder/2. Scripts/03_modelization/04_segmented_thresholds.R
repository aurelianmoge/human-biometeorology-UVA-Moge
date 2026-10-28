# Explore a continuous two-slope same-day threshold model.
# The search is descriptive: uncertainty conditional on a selected breakpoint
# would ignore selection. No breakpoint significance test is reported.
source("2. Scripts/00_shared/distributed_lag.R")
theme <- THEMES[3]
profiles <- list()
predictions <- list()
for (exposure in c("Tmax", "DTR")) {
  data <- read_derived(if (exposure == "DTR") "dtr_daily" else "aggregate_daily")
  base <- make_design(data, exposure)
  # Use the cumulative-model calendar for comparisons; the hinge effect itself
  # remains same-day and must not be interpreted as a cumulative DLNM threshold.
  eligible <- base$eligible
  nuisance <- base$design[, -base$target, drop = FALSE]
  value <- data[[exposure]]
  candidates <- unique(as.numeric(quantile(value[eligible], seq(.15, .85, by = .025))))
  for (threshold in candidates) {
    design <- cbind(nuisance, Slope = value, Hinge = pmax(value - threshold, 0))
    fit <- glm.fit(design[eligible, , drop = FALSE], data$All[eligible], family = quasipoisson())
    if (!fit$converged || any(!is.finite(fit$coefficients))) stop("Invalid segmented fit.")
    coefficients <- tail(fit$coefficients, 2)
    profiles[[length(profiles) + 1]] <- data.frame(Exposure = exposure, Threshold = threshold,
      Poisson_deviance = fit$deviance, Before_slope = coefficients[1], After_slope = sum(coefficients),
      Days = sum(eligible), Interpretation = "retrospective same-day hinge search; no selection-adjusted uncertainty")
  }
}
profiles <- bind_rows(profiles) |> group_by(Exposure) |>
  mutate(Deviance_difference = Poisson_deviance - min(Poisson_deviance), Selected = Poisson_deviance == min(Poisson_deviance)) |>
  ungroup()
write_table(profiles, "segmented_threshold_profile", theme)
save_plot(ggplot(profiles, aes(Threshold, Deviance_difference)) + geom_line(color = "#236B8E") +
  geom_point(data = filter(profiles, Selected), color = "#AF4528", size = 2.5) + facet_wrap(~Exposure, scales = "free") +
  labs(title = "Exploratory same-day breakpoint profiles", x = "Candidate threshold (°C)", y = "Poisson deviance above best fit",
       caption = "Same calendar within each metric; in-sample profile, not a confidence interval or confirmed biological threshold."),
  "segmented_threshold_profiles", theme)
