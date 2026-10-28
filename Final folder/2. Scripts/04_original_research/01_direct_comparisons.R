# Compare fitted responses directly, retaining shared-date dependence.
# Nominal p-values reproduce the old calculations; later calibration concerns
# prevent us from presenting these reruns as newly validated formal discoveries.
source("2. Scripts/00_shared/distributed_lag.R")
theme <- THEMES[4]
fits <- readRDS(file.path(OUTPUT, THEMES[3], "models", "model_index.rds"))
pairs <- list(Female_Male = c("Female", "Male"), Black_White = c("Black", "White"),
              Age75plus_Under50 = c("Age75plus", "Under50"))
direct <- list()
for (pair in names(pairs)) {
  for (exposure in c("Tmax", "DTR")) {
    first <- fits[[paste(exposure, pairs[[pair]][1], "primary", sep = "_")]]
    second <- fits[[paste(exposure, pairs[[pair]][2], "primary", sep = "_")]]
    targets <- if (exposure == "Tmax") first$basis$quantiles[c(2, 4)] else rep(first$basis$quantiles[5], 2)
    windows <- if (exposure == "Tmax") c("Cumulative_0_21", "Acute_0_2") else c("Acute_0_2", "Cumulative_0_21")
    for (index in 1:2) {
      lags <- if (windows[index] == "Acute_0_2") 0:2 else 0:21
      result <- paired_contrast(first, second, model_contrast(first, targets[index], lags), model_contrast(second, targets[index], lags))
      direct[[length(direct) + 1]] <- cbind(data.frame(Comparison = pair, Exposure = exposure,
        Window = windows[index], Exposure_value = targets[index], Days = sum(first$eligible)), result)
    }
  }
}
direct <- bind_rows(direct)
direct$Historical_Holm12_recalculated <- p.adjust(direct$Nominal_P, "holm", n = 12)
direct$Inference_status <- "nominal diagnostic; calibration concerns; historical family arithmetic reproduced"
write_table(direct, "direct_demographic_contrasts", theme)

# Compare metric-specific percentile responses on exactly identical dates.
metric <- list()
for (pair in combn(EXPOSURES[1:3], 2, simplify = FALSE)) {
  first <- fits[[paste(pair[1], "All_primary", sep = "_")]]
  second <- fits[[paste(pair[2], "All_primary", sep = "_")]]
  for (quantile in c(2, 4)) {
    result <- paired_contrast(first, second,
      model_contrast(first, first$basis$quantiles[quantile], 0:21),
      model_contrast(second, second$basis$quantiles[quantile], 0:21))
    metric[[length(metric) + 1]] <- cbind(data.frame(Comparison = paste(pair, collapse = "/"),
      Tail = if (quantile == 2) "P10" else "P90"), result)
  }
}
write_table(bind_rows(metric), "paired_temperature_metrics", theme)

# Independent checks against full-precision previous exports.
reference <- read.csv(file.path(DATA, "reference/results/csv/demographic_direct_contrasts.csv"))
check <- inner_join(direct, reference, by = c("Comparison", "Exposure", "Window"), suffix = c("_new", "_old")) |>
  transmute(Comparison, Exposure, Window, Difference_log_RRR = Log_RR - Log_estimate,
            Difference_SE = SE_new - SE_old,
            Difference_Holm = Historical_Holm12_recalculated - Adjusted_P)
write_table(check, "demographic_reproduction_checks", theme)
stopifnot(nrow(check) == 12, max(abs(check$Difference_log_RRR)) < 1e-7,
          max(abs(check$Difference_SE)) < 1e-7, max(abs(check$Difference_Holm)) < 1e-6)

bandwidth <- list()
residual <- list()
for (exposure in EXPOSURES) {
  fit <- fits[[paste(exposure, "All_primary", sep = "_")]]
  targets <- fit$basis$quantiles[if (exposure == "DTR") c(5, 6) else c(2, 4)]
  for (bandwidth_days in c(14, 28, 56)) {
    bandwidth[[length(bandwidth) + 1]] <- cbind(data.frame(Exposure = exposure, Bandwidth = bandwidth_days, Exposure_value = targets),
      summarize_contrast(fit, model_contrast(fit, targets, 0:21), bandwidth_days))
  }
  for (lag in 1:60) {
    index <- seq_len(length(fit$residual) - lag)
    complete <- is.finite(fit$residual[index]) & is.finite(fit$residual[index + lag])
    residual[[length(residual) + 1]] <- data.frame(Exposure = exposure, Lag = lag, Pairs = sum(complete),
      ACF = cor(fit$residual[index[complete]], fit$residual[(index + lag)[complete]]))
  }
}
write_table(bind_rows(bandwidth), "hac_bandwidth_sensitivity", theme)
write_table(bind_rows(residual), "residual_calendar_autocorrelation", theme)

save_plot(ggplot(direct, aes(RR, Comparison, color = Exposure)) +
  geom_vline(xintercept = 1, linetype = 2) + geom_point(position = position_dodge(width = .45)) +
  geom_errorbar(aes(xmin = Lower, xmax = Upper), width = .15, position = position_dodge(width = .45), orientation = "y") +
  facet_wrap(~Window) + labs(title = "Direct demographic response ratios", x = "Ratio of relative responses", y = NULL,
    caption = "Nominal pointwise intervals; marginal groups; calibration and catchment limitations apply."), "demographic_contrast_forest", theme, 12, 6)
