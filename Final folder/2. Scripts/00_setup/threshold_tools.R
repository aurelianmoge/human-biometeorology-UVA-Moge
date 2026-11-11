# =============================================================================
# threshold_tools.R - helpers shared by the threshold stage (06_thresholds).
# Sourced after config.R and dlnm_tools.R.
#   (1) thresholds read off a pooled DLNM curve: MMT and significance thresholds
#   (2) person-days beyond a regional Tmin percentile or a regional temperature, by group
#   (3) the hockey-stick threshold model (conditional Poisson, fixest): bases, fit, profile,
#       likelihood-weighted threshold draws, likelihood-ratio test of a common threshold
# =============================================================================
suppressPackageStartupMessages(library(fixest)); setFixest_nthreads(1); setFixest_notes(FALSE)

# ---- (1) thresholds on a pooled DLNM curve ---------------------------------------------------------------------
# log RR (and SE) of every grid percentile relative to grid index r
curve_vs <- function(B, cf, V, r) { X <- sweep(B, 2, B[r, ], "-"); list(l = drop(X %*% cf), se = sqrt(pmax(rowSums((X %*% V) * X), 0))) }
mmt_of <- function(B, cf) INRANGE[which.min(drop(B %*% cf)[INRANGE])]          # grid index of the MMT (P1-P99)
# significance threshold: from the extreme (P1 / P99) towards the MMT, the last percentile of the
# uninterrupted run where the lower 95% bound of RR (vs the MMT) exceeds 1; NA if not significant at P1 / P99
sig_thr <- function(l, se, r, side) { s <- (l - 1.96 * se) > 0
  idx <- if (side == "cold") INRANGE[INRANGE < r] else rev(INRANGE[INRANGE > r])
  if (!length(idx) || !s[idx[1]]) return(NA_real_)
  k <- which(!s[idx])[1]; PCT[if (is.na(k)) idx[length(idx)] else idx[k - 1]] }
thr_set <- function(B, cf, V, side) { r <- mmt_of(B, cf); cv <- curve_vs(B, cf, V, r); c(mmt = PCT[r], sig = sig_thr(cv$l, cv$se, r, side)) }
# MMT and significance thresholds of a pair of pooled models (cold: cumulative curve, heat: lag 0-2 curve)
thr_pair <- function(pc, pa) c(cold_mmt = thr_set(pc$basis, pc$coef, pc$vcov, "cold")[["mmt"]], cold_sig = thr_set(pc$basis, pc$coef, pc$vcov, "cold")[["sig"]],
                               heat_mmt = thr_set(pa$basis, pa$coef, pa$vcov, "heat")[["mmt"]], heat_sig = thr_set(pa$basis, pa$coef, pa$vcov, "heat")[["sig"]])

# ---- (2) person-days at risk --------------------------------------------------------------------------------------
# returns a closure-based toolkit for one panel: day-level flags weighted by each group's regional population
exposure_tools <- function(panel, rp) {
  regs <- sort(unique(panel$Region)); d <- panel[, .(Region, year, Tmin_f)]
  d[, pr := ecdf(Tmin_f)(Tmin_f), by = Region]
  pop <- sapply(GROUPS, function(g) merge(d[, .(Region, year)], rp[Group == g, .(Region, year = Year, Population)], by = c("Region", "year"), all.x = TRUE, sort = FALSE)$Population)
  pop[is.na(pop)] <- 0; pop[is.na(d$Tmin_f), ] <- 0
  ridx <- match(d$Region, regs)
  Q <- sapply(regs, function(r) quantile(d[Region == r]$Tmin_f, PCT, na.rm = TRUE))           # PCT x region, deg C
  wg <- sapply(GROUPS, function(g) { x <- rp[Group == g, .(P = sum(Population)), by = Region]; x$P[match(regs, x$Region)] / sum(x$P) })
  list(regs = regs, PY = colSums(pop),                                                         # person-days per group
       # person-days beyond a regional percentile p (NA = no threshold = 0 days)
       pd_pct = function(g, p, side) if (is.na(p)) 0 else sum(pop[, g] * (if (side == "cold") d$pr < p else d$pr > p), na.rm = TRUE),
       # person-days beyond region-specific temperatures thr (deg C, one per region; +-Inf allowed)
       pd_temp = function(g, thr, side, flag = FALSE) { t <- thr[ridx]; f <- !is.na(d$Tmin_f) & (if (side == "cold") d$Tmin_f < t else d$Tmin_f > t); if (flag) f else sum(pop[f, g]) },
       pop = pop,
       pct_C = function(g, p) if (is.na(p)) NA_real_ else sum(wg[, g] * Q[pidx(p), ]),      # percentile -> group's population-weighted deg C
       temp_C = function(g, thr) sum(wg[, g] * thr))                                          # regional temperatures -> population-weighted deg C
}
PARTITIONS <- list(Sex = c("Female", "Male"), Race = c("Black", "White"), Age = c("Under50", "Age50_74", "Age75plus"))

# ---- (3) hockey-stick threshold model ---------------------------------------------------------------------------------
# log E[deaths] = stratum + b_c * sum_l w_l (tau_c - Tmin_{t-l})_+  (lags 0-21)  +  b_h * sum_l v_l (Tmin_{t-l} - tau_h)_+  (lags 0-2)
# strata region x year x month x day of week; tau = a regional Tmin percentile q common to all regions.
tm_setup <- function(panel, QC = seq(.30, .99, .01), QH = seq(.50, .99, .01)) {
  regs <- unique(panel$Region)                                                                 # panel must be sorted by Region, Date
  mk <- function(q, side) { out <- lapply(regs, function(r) { x <- panel[Region == r]$Tmin_f
    cb <- crossbasis(x, lag = if (side == "l") MAX_LAG else 2, argvar = list(fun = "thr", thr.value = quantile(x, q, na.rm = TRUE), side = side),
                     arglag = if (side == "l") list(knots = logknots(MAX_LAG, 3)) else list(fun = "integer"))
    list(m = matrix(unclass(cb), nrow = length(x)), lag = attr(cb, "arglag")) })
    structure(do.call(rbind, lapply(out, `[[`, "m")), arglag = out[[1]]$lag) }
  BC <- lapply(setNames(QC, QC), mk, side = "l"); BH <- lapply(setNames(QH, QH), mk, side = "h")
  cumw <- function(B, L) colSums(do.call(onebasis, c(list(x = 0:L), attr(B, "arglag"))))
  QALL <- sort(unique(round(c(QC, QH), 2)))
  TAUM <- sapply(regs, function(r) quantile(panel[Region == r]$Tmin_f, QALL, na.rm = TRUE)); rownames(TAUM) <- sprintf("%.2f", QALL)
  list(QC = QC, QH = QH, BC = BC, BH = BH, WC = cumw(BC[[1]], MAX_LAG), WH = cumw(BH[[1]], 2), regs = regs, TAUM = TAUM,
       ok = complete.cases(BC[[1]], BH[[1]], panel$Tmin_f), st = paste(panel$Region, panel$year, month(panel$Date), wday(panel$Date)),
       clu = paste(panel$Region, panel$year, month(panel$Date))) }
TAU <- function(S, q) S$TAUM[sprintf("%.2f", q), ]                                            # regional thresholds (deg C)
# one fit at thresholds (qc, qh); se = TRUE adds the dispersion and the cumulative log slopes with clustered SEs
tm_fit <- function(S, y, qc, qh, se = FALSE) {
  X <- cbind(S$BC[[as.character(qc)]], S$BH[[as.character(qh)]]); colnames(X) <- c(paste0("c", seq_len(ncol(S$BC[[1]]))), paste0("h", seq_len(ncol(S$BH[[1]]))))
  d <- data.table(y = y, st = S$st, clu = S$clu, X)[S$ok & !is.na(y)]
  m <- fepois(as.formula(paste("y ~", paste(colnames(X), collapse = " + "), "| st")), d)
  out <- list(dev = deviance(m))
  if (se) { mu <- fitted(m); out$phi <- sum((d$y[obs(m)] - mu)^2 / mu) / (m$nobs - sum(m$fixef_sizes) - length(coef(m)))
    V <- vcov(m, cluster = ~clu); b <- coef(m); ic <- grep("^c", names(b)); ih <- grep("^h", names(b))
    out$bc <- sum(S$WC * b[ic]); out$sc <- sqrt(drop(t(S$WC) %*% V[ic, ic] %*% S$WC))
    out$bh <- sum(S$WH * b[ih]); out$sh <- sqrt(drop(t(S$WH) %*% V[ih, ih] %*% S$WH)) }
  out }
# profile likelihood over the grids: cold | heat at P90, heat | best cold, cold | best heat
tm_profile <- function(S, y) {
  pc0 <- sapply(S$QC, function(q) tm_fit(S, y, q, .90)$dev); qc <- S$QC[which.min(pc0)]
  ph <- sapply(S$QH, function(q) tm_fit(S, y, qc, q)$dev); qh <- S$QH[which.min(ph)]
  pc <- sapply(S$QC, function(q) tm_fit(S, y, q, qh)$dev); qc <- S$QC[which.min(pc)]
  f <- tm_fit(S, y, qc, qh, se = TRUE)
  list(qc = qc, qh = qh, pc = pc, ph = ph, phi = f$phi, bc = f$bc, sc = f$sc, bh = f$bh, sh = f$sh) }
# slope (and SE) at every grid threshold, the other side held at its best value
tm_grid_slopes <- function(S, y, p) {
  p$gc <- rbindlist(lapply(S$QC, function(q) { f <- tm_fit(S, y, q, p$qh, se = TRUE); data.table(q = q, b = f$bc, s = f$sc) }))
  p$gh <- rbindlist(lapply(S$QH, function(q) { f <- tm_fit(S, y, p$qc, q, se = TRUE); data.table(q = q, b = f$bh, s = f$sh) })); p }
tm_ci <- function(prof, Q, phi) range(Q[(prof - min(prof)) / phi <= qchisq(.95, 1)])               # profile-likelihood 95% CI
tm_lr_common <- function(pa, pb, phi) { st <- (min(pa + pb) - min(pa) - min(pb)) / phi; c(LR = st, p = pchisq(st, 1, lower.tail = FALSE)) }
tm_draws <- function(prof, Q, phi, n) sample(Q, n, replace = TRUE, prob = exp(-(prof - min(prof)) / (2 * phi)))   # likelihood-weighted thresholds
# joint draws of (threshold, slope): q ~ likelihood weights, b | q ~ N(b_hat(q), se(q)^2)
tm_draws_qb <- function(p, side, n) { prof <- if (side == "cold") p$pc else p$ph; G <- if (side == "cold") p$gc else p$gh
  i <- sample(seq_len(nrow(G)), n, replace = TRUE, prob = exp(-(prof - min(prof)) / (2 * p$phi))); data.table(q = G$q[i], b = rnorm(n, G$b[i], G$s[i])) }
# regional temperatures at which excess risk k is reached (b <= 0: never; cold -Inf, heat +Inf)
risk_level_temp <- function(S, q, b, k, side) { d <- if (b > 0) log(1 + k) / b else Inf; if (side == "cold") TAU(S, q) - d else TAU(S, q) + d }
# binomial split of a count series into two pseudo-groups with the same true curve
split_null <- function(y, prob) { a <- fifelse(is.na(y), NA_integer_, rbinom(length(y), fifelse(is.na(y), 0L, as.integer(y)), prob)); list(a = a, b = y - a) }
