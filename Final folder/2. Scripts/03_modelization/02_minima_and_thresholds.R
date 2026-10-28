# Minima, RR crossings and published DTR event rules are separate estimands.
source("2. Scripts/00_shared/distributed_lag.R")
theme <- THEMES[3]
fits <- readRDS(file.path(OUTPUT, theme, "models", "model_index.rds"))
minima <- list()
classifications <- list()
uncertainty <- list()
for (id in names(fits)) {
  fit <- fits[[id]]
  if (fit$settings$Specification != "primary") next
  grid <- seq(fit$basis$quantiles[1], fit$basis$quantiles[6], length.out = 1001)
  contrast <- model_contrast(fit, grid, 0:21)
  estimate <- drop(contrast %*% fit$beta)
  minimum_index <- which.min(estimate)
  minimum <- grid[minimum_index]
  # The finite grid and boundary flags make minimizer instability visible.
  qualifying <- exp(estimate) >= 1.05
  upper <- which(qualifying & grid > fit$basis$reference)
  minima[[id]] <- data.frame(Model_ID = id, Exposure = fit$settings$Exposure, Outcome = fit$settings$Outcome,
    Minimum = minimum, Boundary_minimum = minimum_index %in% c(1, length(grid)),
    First_upper_RR105 = if (length(upper)) grid[upper[1]] else NA_real_,
    Upper_crossing_status = if (length(upper)) "reached within domain" else "not reached within domain",
    Domain_lower = min(grid), Domain_upper = max(grid))
  if (fit$settings$Outcome == "All") {
    covariance <- eigen(fit$covariance, symmetric = TRUE)
    if (min(covariance$values) < -1e-8) stop("Materially indefinite covariance.")
    root <- covariance$vectors %*% diag(sqrt(pmax(0, covariance$values)), nrow = length(fit$beta))
    draws <- matrix(rnorm(500 * length(fit$beta)), 500) %*% t(root)
    draws <- sweep(draws, 2, fit$beta, `+`)
    drawn_minimum <- apply(contrast %*% t(draws), 2, function(values) grid[which.min(values)])
    uncertainty[[id]] <- data.frame(Model_ID = id, Draw = seq_along(drawn_minimum), Minimum = drawn_minimum,
      Status = "conditional coefficient-draw minimum; not a calibrated confidence interval")
  }
}
write_table(bind_rows(minima), "fitted_minima_and_crossings", theme)
write_table(bind_rows(uncertainty), "minimum_coefficient_draws", theme)

# Reproduce fixed article rules against audited local extrema, retaining incomplete
# observations and unique event versus forward-lag outcome dates separately.
regional <- read_derived("regional_panel")
rules <- read.csv(file.path(DATA, "reference/metadata/reference_rules.csv"))
events <- list()
counts <- list()
for (row in seq_len(nrow(rules))) {
  rule <- rules[row, ]
  data <- regional[regional$Region == rule$Region, ]
  observed <- is.finite(data$DTR)
  event <- observed & data$DTR >= rule$lower_C
  if (is.finite(rule$upper_C)) event <- event & data$DTR <= rule$upper_C
  outcome <- rep(FALSE, nrow(data))
  # ROA's published lag varies by DTR band; lower bounds are inclusive.
  for (index in which(event)) {
    maximum_lag <- rule$max_lag
    if (rule$Region == "ROA") maximum_lag <- if (data$DTR[index] < 18) 1 else if (data$DTR[index] < 20) 2 else if (data$DTR[index] < 21) 3 else if (data$DTR[index] < 22) 4 else 5
    outcome[index:min(nrow(data), index + maximum_lag)] <- TRUE
  }
  events[[row]] <- data.frame(Region = rule$Region, Date = data$Date, Observed_DTR = observed, DTR = data$DTR,
                              Exposure_event = event, Outcome_union = outcome)
  counts[[row]] <- data.frame(Region = rule$Region, Observed_days = sum(observed), Missing_days = sum(!observed),
    Event_dates = sum(event), Unique_outcome_dates = sum(outcome), Published_event_dates = rule$Published_event_dates,
    Missing_event_lower = sum(event), Missing_event_upper = sum(event) + sum(!observed),
    Interpretation = "Workbook fixed-rule benchmark; source/catchment and weather-window discrepancies remain")
}
write_table(bind_rows(events), "published_rule_date_flags", theme)
write_table(bind_rows(counts), "published_rule_benchmark", theme)

curves <- read.csv(file.path(OUTPUT, theme, "csv/response_curves.csv"))
save_plot(ggplot(curves |> filter(Outcome == "All", Specification == "primary", Window == "Cumulative_0_21"),
  aes(Exposure_value, RR)) + geom_hline(yintercept = 1, linetype = 2) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "#236B8E", alpha = .2) + geom_line(color = "#236B8E") +
  facet_wrap(~Exposure, scales = "free") + labs(title = "Cumulative exposure–mortality responses", x = "Degrees Celsius", y = "RR relative to weather median",
  caption = "Lag 0–21; nominal pointwise intervals. DTR has a different fixed footprint and additional controls."),
  "cumulative_response_curves", theme, 12, 9)
lag <- read.csv(file.path(OUTPUT, theme, "csv/lag_profiles.csv"))
save_plot(ggplot(lag, aes(Lag, RR, color = factor(round(Exposure_value, 1)))) +
  geom_hline(yintercept = 1, linetype = 2) + geom_line() + facet_wrap(~Exposure, scales = "free_y") +
  labs(title = "Timing of modeled responses", color = "Exposure °C", y = "Single-lag relative response",
       caption = "Selected physical anchors differ by metric. These are single-lag contrasts, not cumulative effects."), "lag_response_profiles", theme, 12, 8)
