# Fit an explicit retrospective model registry and retain every requested fit.
source("2. Scripts/00_shared/distributed_lag.R")
theme <- THEMES[3]
aggregate <- read_derived("aggregate_daily")
dtr <- read_derived("dtr_daily")
registry <- list()
add_model <- function(exposure, outcome = "All", specification = "primary", time_df = 112,
                      controls = "full", interaction = FALSE, negative_binomial = FALSE) {
  registry[[length(registry) + 1]] <<- data.frame(
    Model_ID = paste(exposure, outcome, specification, sep = "_"), Exposure = exposure, Outcome = outcome,
    Specification = specification, Time_df = time_df, Controls = controls,
    Interaction = interaction, Negative_binomial = negative_binomial)
}
for (exposure in EXPOSURES) {
  for (outcome in c("All", "Female", "Male", "Black", "White", "Under50", "Age75plus", "Respiratory", "Circulatory")) {
    add_model(exposure, outcome)
  }
  for (df in c(48, 80)) add_model(exposure, specification = paste0("time_df_", df), time_df = df)
  add_model(exposure, specification = "exclude_2020", time_df = 105)
  add_model(exposure, specification = "negative_binomial", negative_binomial = TRUE)
}
for (control in c("temporal_only", "same_day", "temperature")) {
  add_model("DTR", specification = paste0("control_", control), controls = control)
}
add_model("DTR", specification = "temperature_interaction", interaction = TRUE)
registry <- bind_rows(registry)
write_table(registry, "model_registry", theme)

# Hash every fitting dependency. Changing settings, data or code invalidates caches.
dependencies <- c(file.path(DATA, "derived", c("aggregate_daily.csv", "dtr_daily.csv")),
                  "2. Scripts/00_shared/distributed_lag.R", "2. Scripts/00_shared/setup.R",
                  "2. Scripts/03_modelization/01_fit_models.R")
signature <- tools::md5sum(dependencies)
fits <- list()
statuses <- list()
estimates <- list()
curves <- list()
diagnostics <- list()
calendars <- list()
lag_profiles <- list()
windows <- list(Acute_0_2 = 0:2, Days_3_7 = 3:7, Cumulative_0_21 = 0:21)
for (row in seq_len(nrow(registry))) {
  settings <- registry[row, ]
  id <- settings$Model_ID
  message(row, "/", nrow(registry), " ", id)
  path <- file.path(OUTPUT, theme, "models", paste0(id, ".rds"))
  fit <- NULL
  if (file.exists(path)) {
    cached <- readRDS(path)
    if (identical(cached$signature, signature)) fit <- cached$fit
  }
  reused <- !is.null(fit)
  if (is.null(fit)) {
    fit <- tryCatch({
      data <- if (settings$Exposure == "DTR") dtr else aggregate
      design <- make_design(data, settings$Exposure, settings$Outcome, settings$Time_df,
                            settings$Controls, settings$Specification == "exclude_2020", settings$Interaction)
      fit_count_model(design, settings$Negative_binomial)
    }, error = function(error) error)
    if (inherits(fit, "error")) {
      # A previous cache for this request is obsolete and must not look usable.
      if (file.exists(path)) unlink(path)
      statuses[[id]] <- cbind(settings, Status = if (settings$Negative_binomial) "stopped" else "failed", Reason = conditionMessage(fit), Cached = FALSE)
      next
    }
    fit$settings <- settings
    saveRDS(list(signature = signature, fit = fit), path)
  }
  fits[[id]] <- fit
  statuses[[id]] <- cbind(settings, Status = "completed", Reason = "", Cached = reused)
  diagnostics[[id]] <- cbind(settings, Observations = fit$observations, Deaths = fit$deaths,
                             Parameters = fit$parameters, Dispersion = fit$dispersion, Converged = fit$converged)
  calendars[[id]] <- data.frame(Model_ID = id, Date = fit$dates, Eligible = fit$eligible)
  targets <- if (settings$Exposure == "DTR") fit$basis$quantiles[c(5, 6)] else fit$basis$quantiles[c(2, 4)]
  grid <- sort(unique(c(seq(fit$basis$quantiles[1], fit$basis$quantiles[6], length.out = 151), fit$basis$reference)))
  for (window in names(windows)) {
    estimates[[length(estimates) + 1]] <- cbind(settings[rep(1, length(targets)), , drop = FALSE], Window = window, Exposure_value = targets,
       Reference = fit$basis$reference, summarize_contrast(fit, model_contrast(fit, targets, windows[[window]])))
    curves[[length(curves) + 1]] <- cbind(settings[rep(1, length(grid)), , drop = FALSE], Window = window, Exposure_value = grid,
       Reference = fit$basis$reference, summarize_contrast(fit, model_contrast(fit, grid, windows[[window]])))
  }
  if (settings$Outcome == "All" && settings$Specification == "primary") {
    for (lag in 0:21) lag_profiles[[length(lag_profiles) + 1]] <- cbind(settings[rep(1, length(targets)), , drop = FALSE], Lag = lag, Exposure_value = targets,
      summarize_contrast(fit, model_contrast(fit, targets, lag)))
  }
  write_table(bind_rows(statuses), "model_status", theme)
}
write_table(bind_rows(statuses), "model_status", theme)
write_table(bind_rows(estimates), "response_estimates", theme)
write_table(bind_rows(curves), "response_curves", theme)
write_table(bind_rows(diagnostics), "model_diagnostics", theme)
write_table(bind_rows(calendars), "model_calendars", theme)
write_table(bind_rows(lag_profiles), "lag_profiles", theme)
saveRDS(fits, file.path(OUTPUT, theme, "models", "model_index.rds"))
writeLines(capture.output(sessionInfo()), file.path(OUTPUT, theme, "logs", "model_session.txt"))
if (any(vapply(statuses, function(status) status$Status == "failed", logical(1)))) stop("One or more registered count fits failed; inspect model_status.csv.")
