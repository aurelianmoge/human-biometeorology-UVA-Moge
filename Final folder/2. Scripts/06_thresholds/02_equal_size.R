# =============================================================================
# 02_equal_size.R - curve thresholds compared at EQUAL SAMPLE SIZE.
# Significance thresholds depend on precision, hence on a group's number of deaths. The larger
# group's daily deaths are binomially thinned to the smaller group's total (thinning keeps the
# group's true curve and only lowers precision), the two-stage model is refitted, and the MMT and
# significance thresholds are recomputed; R replicates of the thinning.
#   Sex : Female and Male totals already equal -> no thinning
#   Race: White thinned to the Black total;  Age: 50-74 and 75+ thinned to the Under-50 total
#   Reference: all deaths thinned to the Black and to the Under-50 totals
# Size-matched test: p = share of replicates in which the thinned larger group's threshold is at
# least as far from its median as the smaller group's actual threshold.
# Output: tables/06_thresholds/equal_size_replicates.csv, equal_size_thresholds.csv,
#         equal_size_partition_totals.csv
# =============================================================================
source("2. Scripts/00_setup/config.R"); source("2. Scripts/00_setup/dlnm_tools.R"); source("2. Scripts/00_setup/threshold_tools.R")
STAGE <- "06_thresholds"; R <- 100
panel <- load_panel(); setorder(panel, Region, Date); PL <- read_model("pooled_Tmin")
rp <- fread(file.path(DERIVED, "regional_populations.csv"))[Year %in% YEARS]; EX <- exposure_tools(panel, rp)
TOTD <- sapply(GROUPS, function(g) sum(panel[[g]], na.rm = TRUE))
THIN <- data.table(col = c("White_eq", "Age50_74_eq", "Age75plus_eq", "All_eqBlack", "All_eqUnder50"),
                   src = c("White", "Age50_74", "Age75plus", "All", "All"), target = c("Black", "Under50", "Under50", "Black", "Under50"))[, f := TOTD[target] / TOTD[src]]
TT <- c("cold_mmt", "cold_sig", "heat_mmt", "heat_sig")
FULL <- rbindlist(lapply(GROUPS, function(g) as.data.table(as.list(thr_pair(PL[[paste(g, "cum", sep = "|")]], PL[[paste(g, "acute", sep = "|")]])))[, group := g]))

REP <- rbindlist(lapply(seq_len(R), function(k) { set.seed(SEED + 600 + k); d <- copy(panel)
  for (i in seq_len(nrow(THIN))) { y <- d[[THIN$src[i]]]; d[, (THIN$col[i]) := fifelse(is.na(y), NA_integer_, rbinom(.N, fifelse(is.na(y), 0L, as.integer(y)), THIN$f[i]))] }
  s1 <- run_stage1(d, "Tmin_f", groups = THIN$col)
  if (k %% 10 == 0) cat("replicate", k, "\n")
  rbindlist(lapply(THIN$col, function(g) as.data.table(as.list(thr_pair(pool(s1, g, "cum"), pool(s1, g, "acute"))))[, `:=`(group = g, rep = k)])) }))
write_tab(REP, STAGE, "equal_size_replicates")

q <- function(v, p) if (all(is.na(v))) NA_real_ else quantile(v, p, na.rm = TRUE)
SUMM <- rbindlist(lapply(THIN$col, function(g) rbindlist(lapply(TT, function(t) { v <- REP[group == g][[t]]; src <- THIN[col == g]$src; tgt <- THIN[col == g]$target
  data.table(group = g, source = src, thinned_to = tgt, threshold = t, full_size = FULL[group == src][[t]], equal_size_median = q(v, .5),
             equal_size_lo = q(v, .025), equal_size_hi = q(v, .975), share_with_threshold = mean(!is.na(v)), smaller_group_actual = FULL[group == tgt][[t]]) }))))
SUMM[, p_size_matched := mapply(function(g, t, x, src) { if (src == "All") return(NA_real_); v <- REP[group == g][[t]]; fill <- if (startsWith(t, "cold")) 0 else 1
  if (is.na(x)) x <- fill; v[is.na(v)] <- fill; m <- median(v); mean(abs(v - m) >= abs(x - m)) }, group, threshold, smaller_group_actual, source)]
write_tab(SUMM, STAGE, "equal_size_thresholds")

CA <- FULL[group == "All"]
PRT <- rbindlist(lapply(list(Race = list(keep = "Black", thin = c(White = "White_eq")), Age = list(keep = "Under50", thin = c(Age50_74 = "Age50_74_eq", Age75plus = "Age75plus_eq"))), function(pt)
  rbindlist(lapply(TT, function(t) { side <- substr(t, 1, 4); gs <- c(pt$keep, names(pt$thin))
    common <- sum(sapply(gs, EX$pd_pct, p = CA[[t]], side = side)); full <- sum(sapply(gs, function(g) EX$pd_pct(g, FULL[group == g][[t]], side)))
    eq <- sapply(seq_len(R), function(k) EX$pd_pct(pt$keep, FULL[group == pt$keep][[t]], side) + sum(sapply(names(pt$thin), function(g) EX$pd_pct(g, REP[group == pt$thin[[g]] & rep == k][[t]], side))))
    data.table(threshold = t, E_common_millions = common / 1e6, delta_full_size_pct = 100 * (full / common - 1), delta_equal_size_median = 100 * (median(eq) / common - 1),
               delta_equal_size_lo = 100 * (quantile(eq, .025) / common - 1), delta_equal_size_hi = 100 * (quantile(eq, .975) / common - 1)) }))), idcol = "partition")
write_tab(PRT, STAGE, "equal_size_partition_totals")
print(SUMM[, .(group, threshold, full = 100 * full_size, eq = 100 * equal_size_median, lo = 100 * equal_size_lo, hi = 100 * equal_size_hi, has = share_with_threshold, smaller = 100 * smaller_group_actual, p = round(p_size_matched, 2))], nrows = 30)
