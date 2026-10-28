# Numerical verification checks are scientific identities and independent
# checkpoints; passing them does not resolve source provenance or identification.
source("2. Scripts/00_shared/distributed_lag.R")
checks <- list()
record <- function(name, condition, detail = "") {
  checks[[length(checks) + 1]] <<- data.frame(Check = name, Passed = isTRUE(condition), Detail = detail)
}
regional <- read_derived("regional_panel")
record("Twelve complete regional calendars", nrow(regional) == 12 * length(DATES) && length(unique(regional$Region)) == 12)
record("Unique region-date keys", !anyDuplicated(regional[c("Region", "Date")]))
record("Derived DTR identity", max(abs(regional$DTR - regional$Tmax + regional$Tmin), na.rm = TRUE) < 1e-8)
record("Nonnegative observed DTR", all(regional$DTR[is.finite(regional$DTR)] >= 0))
record("Fourteen corroborated date repairs", nrow(read_derived("date_repairs")) == 14)
raw <- read_derived("workbook_cells")
source_flags <- raw |> mutate(Comparable = is.finite(DTR_supplied) & is.finite(Tmin) & is.finite(Tmax),
  Difference = DTR_supplied - (Tmax - Tmin)) |> group_by(Region) |>
  summarise(Comparable = sum(Comparable), Supplied_DTR_disagreements = sum(abs(Difference) > 1e-7, na.rm = TRUE), .groups = "drop")
write_table(source_flags, "supplied_DTR_identity_audit", THEMES[1])
record("PHF supplied-DTR issue independently recovered", source_flags$Supplied_DTR_disagreements[source_flags$Region == "PHF"] == 4236)

population <- read_derived("statewide_population_margins")
for (year in unique(population$Year)) {
  values <- population[population$Year == year, ]
  total <- values$Population[values$Group == "All"]
  for (partition in list(c("Female", "Male"), c("Black", "White", "Residual"), c("Under50", "Age50_74", "Age75plus"))) {
    record(paste("Census partition", year, paste(partition, collapse = "+")), sum(values$Population[values$Group %in% partition]) == total)
  }
}
reference_population <- read.csv(file.path(OUTPUT, THEMES[4], "historical_evidence/3. Results/Stage 3 composition inference/statewide_joint_populations.csv"))
reference_total <- reference_population |> group_by(Year) |> summarise(Reference = sum(Population), .groups = "drop")
compared <- left_join(filter(population, Group == "All"), reference_total, by = "Year")
record("Corrected statewide Census July mapping reproduced", all(compared$Population == compared$Reference) && nrow(compared) == 16)

# Compare every common primary curve anchor and uncertainty calculation.
estimates <- read.csv(file.path(OUTPUT, THEMES[3], "csv/response_estimates.csv")) |>
  filter(Specification == "primary")
reference <- read.csv(file.path(DATA, "reference/results/csv/estimates.csv")) |>
  filter(Spec == "primary")
compared <- inner_join(estimates, reference,
  by = c("Exposure", "Outcome", "Window", "Exposure_value" = "Contrast_value"), suffix = c("_new", "_old"))
# Some CSV decimal representations differ at the final digit; join physical
# values through rounded keys as well, never silently drop a missing checkpoint.
estimates$Physical_key <- round(estimates$Exposure_value, 7)
reference$Physical_key <- round(reference$Contrast_value, 7)
compared <- inner_join(estimates, reference, by = c("Exposure", "Outcome", "Window", "Physical_key"), suffix = c("_new", "_old")) |>
  transmute(Exposure, Outcome, Window, Physical_key, Log_RR_difference = Log_RR - Log_estimate,
            SE_difference = SE_new - SE_old, RR_difference = RR - Estimate)
write_table(compared, "primary_model_reproduction_checks", THEMES[3])
record("Primary model checkpoints present", nrow(compared) >= 100, paste(nrow(compared), "anchor/window checks"))
record("Primary log-response reproduction", max(abs(compared$Log_RR_difference)) < 1e-7)
record("Primary covariance reproduction", max(abs(compared$SE_difference)) < 1e-7)

county <- read.csv(file.path(OUTPUT, THEMES[4], "csv/county_social_gradient_regressions.csv"))
county_reference <- read.csv(file.path(DATA, "reference/results/csv/county_poverty_checkpoint.csv"))
selected <- filter(county, Outcome == "Age_adjusted_rate", Specification == "poverty", Term == "Poverty5")
record("County poverty coefficient reproduction", abs(selected$Beta - county_reference$Beta) < 1e-8)
record("County poverty clustered covariance reproduction", abs(selected$SE - county_reference$SE) < 1e-8)

fits <- readRDS(file.path(OUTPUT, THEMES[3], "models/model_index.rds"))
for (id in names(fits)) {
  fit <- fits[[id]]
  record(paste("Valid count fit", id), fit$converged && all(is.finite(fit$beta)) && fit$observations > fit$parameters)
  record(paste("Covariance PSD", id), min(eigen(fit$covariance, symmetric = TRUE, only.values = TRUE)$values) > -1e-8)
  contrast <- model_contrast(fit, fit$basis$reference, 0:21)
  record(paste("RR reference identity", id), max(abs(contrast %*% fit$beta)) < 1e-10)
}
status <- read.csv(file.path(OUTPUT, THEMES[3], "csv/model_status.csv"))
record("All requested count fits accounted for", nrow(status) == 56 && !any(status$Status == "failed"))
record("Stopped fits explicitly labeled", all(status$Negative_binomial[status$Status == "stopped"]))

classification <- read.csv(file.path(OUTPUT, THEMES[4], "csv/pooled_vs_own_date_classification.csv"))
record("Classification date partition", all(with(classification, Own_only + Pooled_only + Both + Neither == Assessed)))
record("Classification signed net identity", all(with(classification, Own_only - Pooled_only == Net_dates)))
record("Classification covers every registered combination", nrow(classification) == 4 * 6 * 2 * 3)

# Recompute the archived reliability frequencies from the original replicate rows.
path <- file.path(OUTPUT, THEMES[4], "historical_evidence/population_risk_comparison_v1/method_reliability/results/simulation_replicates.csv")
simulations <- read.csv(path)
calibration <- simulations |> filter(Status == "completed") |> group_by(Scenario) |>
  summarise(Successful = n(), Coverage = mean(Coverage), Excludes_zero = mean(Excludes_zero),
            Mean_width = mean(Width), Stable_beyond_1pp = mean(Stable_beyond_1pp), .groups = "drop")
write_table(calibration, "archived_calibration_recomputed", THEMES[4])
record("Archived equal-scenario undercoverage reproduced", abs(calibration$Coverage[calibration$Scenario == "Equal"] - .885) < 1e-10)

result <- bind_rows(checks)
write.csv(result, file.path(DATA, "metadata", "verification_checks.csv"), row.names = FALSE)
print(result |> count(Passed))
if (!all(result$Passed)) stop("Verification failure: ", paste(result$Check[!result$Passed], collapse = "; "))
