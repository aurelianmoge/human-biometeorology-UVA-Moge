# =============================================================================
# 05_fill_weather.R - fill the missing station Tmin / Tmax (deaths are never imputed).
#
# 544 of 70,128 station-days lack Tmin (0.8%): one 244-day block at OKV (Dec 2012 -
# Aug 2013), 198 scattered EZF days (mostly 2017-18), a 17-day PHF block (2009), and
# <= 18 days elsewhere.  Davis et al. (2023) filled gaps from GHCN / PRISM; we do the same:
#   M3  same station in GHCN-Daily, if the workbook series IS that station
#       (mean |difference| < 0.1 C and r > 0.998 on common days);
#   M2  nearby GHCN station(s), month-specific linear calibration;
#   M1  mean of the other 11 Virginia stations, month-specific linear calibration.
# Choice per station and variable: M3 when available, otherwise the better of M1 / M2
# in a masked validation on OBSERVED days (mean RMSE over two masks: 2,000 random days
# and one hidden 244-day block, the size of the OKV gap).
# Input : 1. Data/derived/regional_panel.csv; 1. Data/raw/ghcn/ (script 04)
# Output: 1. Data/derived/weather_filled.csv (Tmin, Tmax + source flags);
#         3. Outputs/tables/01_data/weather_station_identity.csv, weather_fill_validation.csv,
#         weather_fill_sources.csv
# =============================================================================
source("2. Scripts/00_setup/config.R")
STAGE <- "01_data"; set.seed(SEED)
panel <- fread(file.path(DERIVED, "regional_panel.csv"))[, .(Region, Date = as.IDate(Date), Tmin, Tmax)]
gdir <- file.path(DATA, "raw", "ghcn"); gs <- fread(file.path(gdir, "ghcn_stations.csv"))
G <- rbindlist(lapply(unique(gs$Station), function(s) { d <- fread(file.path(gdir, paste0(s, ".csv")))
  if (!nrow(d)) return(NULL); d[, .(Station = STATION, Date = as.IDate(DATE), Tmin = as.numeric(TMIN), Tmax = as.numeric(TMAX))] }), fill = TRUE)
regions <- sort(unique(panel$Region))
W <- list(Tmin = dcast(panel, Date ~ Region, value.var = "Tmin"), Tmax = dcast(panel, Date ~ Region, value.var = "Tmax"))
dates <- W$Tmin$Date; mon <- factor(month(dates))
ghcn <- function(st, var) { g <- G[Station == st]; g[[var]][match(dates, g$Date)] }

# ---- (1) is the workbook series the GHCN candidate station? --------------------------------
ident <- rbindlist(lapply(regions, function(r) { st <- gs[Region == r & Role == "candidate", Station]
  y <- W$Tmin[[r]]; g <- ghcn(st, "Tmin"); ok <- is.finite(y) & is.finite(g)
  data.table(Region = r, Station = st, common_days = sum(ok), mean_abs_diff = mean(abs(y[ok] - g[ok])), correlation = cor(y[ok], g[ok]),
             workbook_missing = sum(!is.finite(y)), covered_by_station = sum(!is.finite(y) & is.finite(g))) }))
ident[, same_station := mean_abs_diff < .1 & correlation > .998]
write_tab(ident, STAGE, "weather_station_identity")

# ---- (2) candidate fills -----------------------------------------------------------------------
calib <- function(y, x, train) { ok <- train & is.finite(y) & is.finite(x); p <- rep(NA_real_, length(y))
  if (sum(ok) < 200) return(p); f <- lm(y ~ mon * x, subset = ok); i <- is.finite(x)
  p[i] <- predict(f, data.frame(mon = mon, x = x)[i, ]); p }
fills <- function(r, var, train) {
  y <- W[[var]][[r]]
  oth <- rowMeans(as.matrix(W[[var]][, setdiff(regions, r), with = FALSE]), na.rm = TRUE)
  M1 <- calib(y, oth, train)
  don <- lapply(gs[Region == r & Role == "donor", Station], ghcn, var = var)
  cr <- sapply(don, function(x) { ok <- train & is.finite(y) & is.finite(x); if (sum(ok) < 200) NA_real_ else cor(y[ok], x[ok]) })
  M2 <- rep(NA_real_, length(y))
  for (j in order(-cr)[seq_len(sum(!is.na(cr)))]) { p <- calib(y, don[[j]], train); M2[is.na(M2)] <- p[is.na(M2)] }
  list(y = y, M1 = M1, M2 = M2, M3 = ghcn(gs[Region == r & Role == "candidate", Station], var))
}

# ---- (3) masked validation on observed days ---------------------------------------------------------
val <- rbindlist(lapply(regions, function(r) rbindlist(lapply(c("Tmin", "Tmax"), function(var) {
  y <- W[[var]][[r]]; obs <- which(is.finite(y))
  masks <- list(random = sample(obs, min(2000, length(obs))), block = { s <- sample(seq(400, length(y) - 650), 1); intersect(s:(s + 243), obs) })
  rbindlist(lapply(names(masks), function(mk) { h <- masks[[mk]]; tr <- rep(TRUE, length(y)); tr[h] <- FALSE; M <- fills(r, var, tr)
    rbindlist(lapply(c("M1", "M2", "M3"), function(k) { e <- M[[k]][h] - y[h]
      data.table(Region = r, variable = var, mask = mk, method = k, coverage = mean(is.finite(e)),
                 rmse = sqrt(mean(e^2, na.rm = TRUE)), bias = mean(e, na.rm = TRUE)) })) })) }))))
write_tab(val, STAGE, "weather_fill_validation")

# ---- (4) final filled series ------------------------------------------------------------------------
lab <- c(M1 = "other Virginia stations (calibrated)", M2 = "nearby GHCN station (calibrated)", M3 = "same station, GHCN-Daily")
out <- rbindlist(lapply(regions, function(r) {
  one <- lapply(c("Tmin", "Tmax"), function(var) {
    M <- fills(r, var, rep(TRUE, length(dates))); v <- M$y; src <- ifelse(is.finite(v), "observed", NA_character_)
    best <- val[Region == r & variable == var & method %in% c("M1", "M2"), .(rmse = mean(rmse)), by = method][order(rmse), method]
    for (k in c(if (ident[Region == r, same_station]) "M3", best)) { i <- !is.finite(v) & is.finite(M[[k]]); v[i] <- M[[k]][i]; src[i] <- lab[[k]] }
    list(v = v, src = src) })
  data.table(Region = r, Date = dates, Tmin = one[[1]]$v, Tmax = one[[2]]$v, Tmin_source = one[[1]]$src, Tmax_source = one[[2]]$src) }))
swap <- out[, Tmax < Tmin & (Tmin_source != "observed" | Tmax_source != "observed")]     # filled pair crossing: order them
out[which(swap), c("Tmin", "Tmax") := .(pmin(Tmin, Tmax), pmax(Tmin, Tmax))]
fwrite(out, file.path(DERIVED, "weather_filled.csv"))
src <- dcast(out[Tmin_source != "observed", .N, by = .(Region, Tmin_source)], Region ~ Tmin_source, value.var = "N", fill = 0)
write_tab(src, STAGE, "weather_fill_sources")
cat("Tmin sources:\n"); print(out[, .N, by = Tmin_source]); print(val[mask == "block" & variable == "Tmin", .(rmse = round(mean(rmse, na.rm = TRUE), 2)), by = method])
