# =============================================================================
# config.R - settings shared by every R script of the pipeline.
# Every script is run from the "Final folder" root (run_all.R does this).
# =============================================================================
if (!dir.exists("1. Data/raw") || !file.exists("2. Scripts/00_setup/config.R"))
  stop("Run from the 'Final folder' root, e.g.  Rscript '2. Scripts/run_all.R'")

pkgs <- c("data.table", "dlnm", "mixmeta", "splines", "MASS", "parallel", "ggplot2", "fixest", "dplyr", "tidyr")
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Missing R packages: ", paste(miss, collapse = ", "))
suppressPackageStartupMessages({
  library(data.table); library(dlnm); library(mixmeta); library(splines); library(MASS); library(parallel); library(ggplot2)
})

# ---- paths ---------------------------------------------------------------------
DATA    <- "1. Data"
DERIVED <- file.path(DATA, "derived")
META    <- file.path(DATA, "reference", "metadata")
OUT     <- "3. Outputs"
out_path <- function(kind, stage, file) {            # kind: tables | figures | models | logs
  d <- file.path(OUT, kind, stage); dir.create(d, recursive = TRUE, showWarnings = FALSE); file.path(d, file)
}
write_tab  <- function(x, stage, name) fwrite(x, out_path("tables", stage, paste0(name, ".csv")))
save_model <- function(x, name) saveRDS(x, out_path("models", "", paste0(name, ".rds")))
read_model <- function(name) readRDS(file.path(OUT, "models", paste0(name, ".rds")))
save_fig   <- function(g, stage, name, w = 9, h = 4.5) ggsave(out_path("figures", stage, paste0(name, ".png")), g, width = w, height = h, dpi = 200)

# ---- study constants ----------------------------------------------------------------
NCORES  <- max(1L, min(8L, parallel::detectCores() - 1L))
YEARS   <- 2005:2020
GROUPS  <- c("All", "Female", "Male", "Black", "White", "Under50", "Age50_74", "Age75plus")
GROUP_LABEL <- c(All = "All deaths", Female = "Female", Male = "Male", Black = "Black", White = "White",
                 Under50 = "Under 50", Age50_74 = "50–74", Age75plus = "75+")
PAIRS   <- list(c("Female", "Male"), c("Black", "White"), c("Age75plus", "Under50"), c("Age75plus", "Age50_74"), c("Age50_74", "Under50"))
MAX_LAG <- 21
DF_YEAR <- 7                                         # seasonal/long-term spline df per year
PCT     <- c(seq(0, 1, .1), 2:98, seq(99, 100, .1)) / 100   # percentile grid for curves
pidx    <- function(p) which.min(abs(PCT - p))
INRANGE <- which(PCT >= .01 & PCT <= .99)            # MMT searched on P1-P99 (Tobias et al. 2017)
SEED    <- 20261110L

# ---- figure style (colour-blind-safe categorical palette) ------------------------------
COL <- c(blue = "#2a78d6", orange = "#eb6834", aqua = "#1baf7a", grey = "#8a8984", lgrey = "#e6e5e1", ink2 = "#52514e")
GROUP_COL <- c(All = "#0b0b0b", Female = COL[["orange"]], Male = COL[["blue"]], Black = COL[["orange"]], White = COL[["blue"]],
               Under50 = COL[["aqua"]], Age50_74 = COL[["blue"]], Age75plus = COL[["orange"]])
theme_set(theme_minimal(base_size = 11) + theme(
  panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = COL[["lgrey"]], linewidth = .3),
  plot.background = element_rect(fill = "#fcfcfb", colour = NA), legend.position = "top",
  plot.title = element_text(face = "bold", size = 12), plot.subtitle = element_text(colour = COL[["ink2"]], size = 9.5),
  plot.title.position = "plot", strip.text = element_text(face = "bold", hjust = 0)))
