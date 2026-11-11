# =============================================================================
# run_all.R - runs the whole pipeline in order, from the "Final folder" root.
#   Rscript "2. Scripts/run_all.R"                    # everything except the null simulation (~1 h 10 min;
#                                                     #   stage 06 thresholds ~40 min of it)
#   Rscript "2. Scripts/run_all.R" --with-simulation  # + 300-replicate null simulation (~100 min more)
#   Rscript "2. Scripts/run_all.R" --from 03          # restart at a stage (01 ... 08)
# Each step's console output goes to 3. Outputs/logs/<step>.log; the run stops at the first error.
# =============================================================================
if (!dir.exists("1. Data/raw")) stop("Run from the 'Final folder' root.")
args <- commandArgs(TRUE)
from <- if ("--from" %in% args) args[which(args == "--from") + 1] else "01"
steps <- c(
  "01_data/01_extract_workbooks.py", "01_data/02_census_populations.py", "01_data/03_assemble_panel.R",
  "01_data/04_fetch_ghcn.py", "01_data/05_fill_weather.R", "01_data/06_county_crosswalk.R", "01_data/07_analysis_panel.R",
  "02_descriptives/01_descriptives.R", "02_descriptives/02_exposure_choice.R", "02_descriptives/03_exposure_by_group.R",
  "03_temperature/01_two_stage_dlnm.R", "03_temperature/02_group_contrasts.R", "03_temperature/03_sensitivity_grid.R",
  "04_dtr/01_two_stage_dtr.R",
  "05_group_differences/01_published_design.R", "05_group_differences/02_conditional_poisson.R", "05_group_differences/03_attributable_burden.R",
  "05_group_differences/04_burden_decomposition.R",
  "06_thresholds/01_mmt_and_significance.R", "06_thresholds/02_equal_size.R", "06_thresholds/03_threshold_model.R", "06_thresholds/04_risk_level.R",
  if ("--with-simulation" %in% args) "07_methods_supplement/01_null_simulation.R", "07_methods_supplement/02_real_vs_null.R", "07_methods_supplement/03_metric_power.R",
  "08_figures/01_main_figures.R", "08_figures/02_paper_assets.R")
steps <- steps[substr(steps, 1, 2) >= from]
dir.create("3. Outputs/logs", recursive = TRUE, showWarnings = FALSE)
for (s in steps) {
  log <- file.path("3. Outputs/logs", sub("\\.(R|py)$", ".log", basename(s))); t0 <- Sys.time()
  cmd <- if (grepl("\\.py$", s)) "python3" else "Rscript"
  status <- system2(cmd, c(if (cmd == "python3") "-I", shQuote(file.path("2. Scripts", s))), stdout = log, stderr = log)
  cat(sprintf("%-50s %s  (%.1f min)\n", s, if (status == 0) "ok" else "FAILED", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  if (status != 0) stop("Step failed: ", s, " - see ", log)
}
writeLines(capture.output(sessionInfo()), "3. Outputs/logs/sessionInfo.txt")
