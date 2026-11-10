# =============================================================================
# 03_sensitivity_grid.R - one-factor-at-a-time sensitivity of the main model.
#   seasonal df/yr {4, 7*, 10}; maximum lag {7, 14, 21*, 28}; exposure {Tmin*, Tmean, Tmax,
#   1 pm temperature}; exposure basis {ns P10/75/90*, quadratic B-spline (Davis / DeGuzman),
#   ns P10/50/90}; no weather gap-filling; excluding 2020 (COVID-19); excluding ROA; excluding SHD.
# For each specification: pooled cumulative RR at P1 / P2.5 / P99 (reference = all-cause MMT),
# MMT by margin, and for each pair the RRR at P1 / P99 and the whole-curve test.
# Output: tables/03_temperature/sensitivity_RR.csv, sensitivity_contrasts.csv
# =============================================================================
source("2. Scripts/00_setup/config.R"); source("2. Scripts/00_setup/dlnm_tools.R")
STAGE <- "03_temperature"
panel <- load_panel()
`%or%` <- function(a, b) if (is.null(a)) b else a
specs <- list(base = list(), df4 = list(df_year = 4), df10 = list(df_year = 10), lag7 = list(lag = 7), lag14 = list(lag = 14), lag28 = list(lag = 28),
              Tmean = list(expo = "Tmean_f"), Tmax = list(expo = "Tmax_f"), Afternoon = list(expo = "Tafternoon_f"),
              bs2 = list(vtype = "bs2"), ns_median = list(vtype = "nsm"), no_fill = list(expo = "Tmin"),
              no2020 = list(sub = quote(year < 2020)), noROA = list(sub = quote(Region != "ROA")),
              noSHD = list(sub = quote(Region != "SHD")))                  # SHD workbook Tmin departs from GHCN by ~5 C
rows <- list(); tests <- list()
for (nm in names(specs)) {
  sp <- specs[[nm]]; vt <- sp$vtype %or% "ns3"
  dat <- if (is.null(sp$sub)) panel else panel[eval(sp$sub)]
  s1 <- run_stage1(dat, sp$expo %or% "Tmin_f", lag = sp$lag %or% MAX_LAG, df_year = sp$df_year %or% DF_YEAR, vtype = vt, windows = list(cum = c(0, NA)))
  pl <- setNames(lapply(GROUPS, function(g) pool(s1, g, "cum", vt)), GROUPS); ref <- pl$All$mmt
  for (g in GROUPS) for (q in c(.01, .025, .99)) { l <- lrr(pl[[g]], q, ref)
    rows[[length(rows) + 1]] <- data.table(spec = nm, group = g, pct = 100 * q, RR = exp(l[1]), lo = exp(l[1] - 1.96 * l[2]), hi = exp(l[1] + 1.96 * l[2]),
                                           MMT = 100 * pl[[g]]$mmt, I2 = pl[[g]]$i2) }
  for (pr in PAIRS) { ct <- curve_test(s1, pr)
    for (q in c(.01, .99)) { r <- rr_ratio(pl[[pr[1]]], pl[[pr[2]]], q, ref)
      tests[[length(tests) + 1]] <- data.table(spec = nm, pair = paste(pr, collapse = "/"), pct = 100 * q, RRR = r[["RRR"]], lo = r[["lo"]], hi = r[["hi"]],
        p_point = r[["p"]], degenerate = q == .99 && ref >= .99, p_curve = ct[["p"]]) } }
  cat(nm, "done\n")
}
write_tab(rbindlist(rows), STAGE, "sensitivity_RR"); write_tab(rbindlist(tests), STAGE, "sensitivity_contrasts")
print(dcast(rbindlist(rows)[group == "All"], spec ~ pct, value.var = "RR")[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
