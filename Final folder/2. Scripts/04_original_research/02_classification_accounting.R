# Compare pooled and subgroup classification rules on identical observed dates.
# Coefficient thresholds classify dates; they are not identified individual risks.
source("2. Scripts/00_shared/distributed_lag.R")
theme <- THEMES[4]
fits <- readRDS(file.path(OUTPUT, THEMES[3], "models", "model_index.rds"))
aggregate <- read_derived("aggregate_daily")
dtr <- read_derived("dtr_daily")
rows <- list()
date_rows <- list()
for (exposure in EXPOSURES) {
  data <- if (exposure == "DTR") dtr else aggregate
  pooled <- fits[[paste(exposure, "All_primary", sep = "_")]]
  values <- data[[exposure]]
  supported <- is.finite(values) & values >= pooled$basis$quantiles[1] & values <= pooled$basis$quantiles[6]
  # Assessment uses observed weather within a common physical domain; it is
  # explicitly different from the lag-history-dependent model fitting calendar.
  eligible <- which(supported)
  for (window in c("Acute_0_2", "Cumulative_0_21")) {
    lags <- if (window == "Acute_0_2") 0:2 else 0:21
    pooled_rr <- exp(drop(model_contrast(pooled, values[eligible], lags) %*% pooled$beta))
    for (group in c("Female", "Male", "Black", "White", "Age75plus", "Under50")) {
      own <- fits[[paste(exposure, group, "primary", sep = "_")]]
      own_rr <- exp(drop(model_contrast(own, values[eligible], lags) %*% own$beta))
      for (cutoff in c(1.02, 1.05, 1.10)) {
        pooled_flag <- pooled_rr >= cutoff
        own_flag <- own_rr >= cutoff
        counts <- data.frame(Exposure = exposure, Group = group, Window = window, Rule = paste0("RR", cutoff),
          Assessed = length(eligible), Pooled_dates = sum(pooled_flag), Own_dates = sum(own_flag),
          Own_only = sum(own_flag & !pooled_flag), Pooled_only = sum(pooled_flag & !own_flag),
          Both = sum(own_flag & pooled_flag), Neither = sum(!own_flag & !pooled_flag),
          Net_dates = sum(own_flag) - sum(pooled_flag),
          Change_percentage_points = 100 * mean(as.integer(own_flag) - as.integer(pooled_flag)))
        stopifnot(counts$Own_only - counts$Pooled_only == counts$Net_dates,
                  counts$Own_only + counts$Pooled_only + counts$Both + counts$Neither == counts$Assessed)
        rows[[length(rows) + 1]] <- counts
        if (cutoff == 1.05) date_rows[[length(date_rows) + 1]] <- data.frame(Exposure = exposure, Group = group,
          Window = window, Date = data$Date[eligible], Value = values[eligible], Pooled = pooled_flag, Own = own_flag)
      }
    }
  }
}
accounting <- bind_rows(rows)
write_table(accounting, "pooled_vs_own_date_classification", theme)
write_table(bind_rows(date_rows), "classification_date_flags_RR105", theme)
write_table(bind_rows(date_rows) |> mutate(Year = format(Date, "%Y"), Month = format(Date, "%m")) |>
  group_by(Exposure, Group, Window, Year, Month) |>
  summarise(Assessed = n(), Pooled = sum(Pooled), Own = sum(Own), .groups = "drop"), "classification_monthly", theme)
save_plot(ggplot(accounting |> filter(Rule == "RR1.05"), aes(Exposure, Group, fill = Change_percentage_points)) +
  geom_tile() + geom_text(aes(label = sprintf("%.1f", Change_percentage_points)), size = 3) + facet_wrap(~Window) +
  scale_fill_gradient2(low = "#236B8E", mid = "white", high = "#AF4528", name = "Change pp") +
  labs(title = "Changing from pooled to subgroup RR ≥ 1.05 rules", x = NULL, y = NULL,
       caption = "Signed change in observed date fraction, full qualifying set. Uncertainty in thresholding is not calibrated."),
  "classification_change_heatmap", theme, 12, 7)
