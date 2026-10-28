# Shared paths and explicit project settings.
# Run the scripts from the clean project root, or use run.py.

if (!file.exists("2. Scripts/00_shared/setup.R") || !dir.exists("1. Data/raw")) {
  stop("Run from the Clean folder root, or use its 2. Scripts/run.py entry point.")
}

required_packages <- c("dplyr", "tidyr", "splines", "ggplot2", "gridExtra", "mgcv", "MASS")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "), ". Install deliberately before running; this pipeline never upgrades packages.")
}
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(splines)
  library(ggplot2)
})

DATA <- "1. Data"
OUTPUT <- "3. Outputs"
THEMES <- c("01_data_description", "02_descriptive_statistics", "03_modelization", "04_original_research")
DATES <- seq(as.Date("2005-01-01"), as.Date("2020-12-31"), by = "day")
GROUPS <- c("All", "Black", "White", "Female", "Male", "Under50", "Age50_74", "Age75plus", "Respiratory", "Circulatory")
EXPOSURES <- c("Tmin", "Tmax", "Tmean_midrange", "DTR")
SEED <- 20261028L
set.seed(SEED)

for (theme in THEMES) {
  for (format in c("csv", "png", "tables_png", "models", "logs")) {
    dir.create(file.path(OUTPUT, theme, format), recursive = TRUE, showWarnings = FALSE)
  }
}
write_table <- function(data, name, theme) {
  write.csv(data, file.path(OUTPUT, theme, "csv", paste0(name, ".csv")), row.names = FALSE, na = "")
}
save_plot <- function(plot, name, theme, width = 10, height = 6) {
  ggsave(file.path(OUTPUT, theme, "png", paste0(name, ".png")), plot,
         width = width, height = height, dpi = 220, bg = "white", limitsize = FALSE)
}
read_derived <- function(name) {
  data <- read.csv(file.path(DATA, "derived", paste0(name, ".csv")), check.names = FALSE)
  if ("Date" %in% names(data)) data$Date <- as.Date(data$Date)
  data
}
plot_theme <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.caption = element_text(hjust = 0), legend.position = "bottom")
theme_set(plot_theme)
