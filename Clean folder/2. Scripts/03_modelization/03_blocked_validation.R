# Buffered internal holdouts assess prediction, not causal validity or calibration.
source("2. Scripts/00_shared/distributed_lag.R")
theme <- THEMES[3]
results <- list()
stopped <- list()
for (exposure in c("Tmax", "DTR")) {
  data <- read_derived(if (exposure == "DTR") "dtr_daily" else "aggregate_daily")
  design <- make_design(data, exposure)
  for (start_text in c("2008-04-01", "2011-07-01", "2014-10-01", "2017-01-01")) {
    start <- as.Date(start_text)
    test <- design$eligible & data$Date >= start & data$Date < start + 28
    train <- design$eligible & (data$Date < start - 42 | data$Date >= start + 28 + 42)
    if (sum(test) < 5) {
      stopped[[length(stopped) + 1]] <- data.frame(Exposure = exposure, Block = start_text, Reason = "fewer than 5 eligible test days")
      next
    }
    # Fit temporal spline knots on training dates only; outcome data never select bases.
    training_spline <- ns(as.numeric(data$Date[train]), df = 112)
    x <- design$design
    x[, 2:113] <- predict(training_spline, as.numeric(data$Date))
    for (candidate in c("exposure", "nuisance_only")) {
      columns <- if (candidate == "exposure") seq_len(ncol(x)) else setdiff(seq_len(ncol(x)), design$target)
      train_x <- x[train, columns, drop = FALSE]
      if (qr(train_x)$rank < ncol(train_x)) {
        stopped[[length(stopped) + 1]] <- data.frame(Exposure = exposure, Block = start_text, Reason = paste(candidate, "rank deficient"))
        next
      }
      fit <- glm.fit(train_x, data$All[train], family = poisson())
      if (!fit$converged || any(!is.finite(fit$coefficients))) stop("Invalid held-out fit.")
      mean <- exp(drop(x[test, columns, drop = FALSE] %*% fit$coefficients))
      observed <- data$All[test]
      deviance <- 2 * (ifelse(observed == 0, 0, observed * log(pmax(observed, 1e-20) / mean)) - observed + mean)
      results[[length(results) + 1]] <- data.frame(Exposure = exposure, Candidate = candidate, Block = start_text,
        Date = data$Date[test], Observed = observed, Predicted = mean, Poisson_deviance = deviance,
        Absolute_error = abs(observed - mean), Training_days = sum(train), Buffer_days = 42)
    }
  }
}
results <- bind_rows(results)
write_table(results, "blocked_validation_daily", theme)
write_table(results |> group_by(Exposure, Candidate) |>
  summarise(Test_observations = n(), Mean_deviance = mean(Poisson_deviance), MAE = mean(Absolute_error), .groups = "drop"),
  "blocked_validation_scores", theme)
if (length(stopped)) write_table(bind_rows(stopped), "blocked_validation_omissions", theme)
