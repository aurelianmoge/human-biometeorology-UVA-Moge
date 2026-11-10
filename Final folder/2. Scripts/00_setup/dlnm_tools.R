# =============================================================================
# dlnm_tools.R - two-stage DLNM toolkit used by stages 03-07.
#   Stage 1 (fit_region / run_stage1): quasi-Poisson DLNM per region x group,
#     reduced to overall cumulative curves (Gasparrini & Armstrong 2013).
#   Stage 2 (pool): multivariate REML meta-analysis (mixmeta) on the percentile
#     scale; pooled curve evaluated on the average-region temperature
#     distribution (Gasparrini et al. 2015).
# =============================================================================

load_panel <- function() readRDS(file.path(DERIVED, "analysis_panel.rds"))

# exposure-response basis; knots at the region's own percentiles
argvar_spec <- function(x, type = "ns3") switch(type,
  ns3 = list(fun = "ns", knots = quantile(x, c(.10, .75, .90), na.rm = TRUE), Boundary.knots = range(x, na.rm = TRUE)),
  nsm = list(fun = "ns", knots = quantile(x, c(.10, .50, .90), na.rm = TRUE), Boundary.knots = range(x, na.rm = TRUE)),
  bs2 = list(fun = "bs", degree = 2, knots = seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = 4)[2:3],
             Boundary.knots = range(x, na.rm = TRUE)))

# Stage 1 for one region/outcome. `control`: optional second cross-basis (e.g. Tmean in DTR models);
# `extra`: optional formula terms (character) added to the model.
fit_region <- function(dat, expo, outcome, lag = MAX_LAG, df_year = DF_YEAR, vtype = "ns3", control = NULL, extra = NULL,
                       windows = list(cum = c(0, NA), acute = c(0, 2), delayed = c(3, NA)), lag_at = NULL) {
  x  <- dat[[expo]]
  cb <- crossbasis(x, lag = lag, argvar = argvar_spec(x, vtype), arglag = list(knots = logknots(lag, 3)))
  rhs <- c("cb", if (!is.null(control)) "cc", if (!is.null(extra)) extra, sprintf("ns(time, df = %d)", df_year * length(unique(dat$year))), "dow")
  if (!is.null(control)) { z <- dat[[control]]; cc <- crossbasis(z, lag = lag, argvar = argvar_spec(z, "ns3"), arglag = list(knots = logknots(lag, 3))) }
  m <- glm(as.formula(paste("y ~", paste(rhs, collapse = " + "))), data = cbind(dat, y = dat[[outcome]]),
           family = quasipoisson, na.action = "na.exclude")
  cen <- median(x, na.rm = TRUE)
  red <- lapply(windows, function(w) { w[is.na(w)] <- lag; if (w[1] > lag) return(NULL)
    r <- crossreduce(cb, m, type = "overall", lag = w, cen = cen); list(coef = coef(r), vcov = vcov(r)) })
  if (!is.null(lag_at)) { r <- crossreduce(cb, m, type = "var", value = quantile(x, lag_at, na.rm = TRUE), cen = cen)
    red$lag_response <- list(coef = coef(r), vcov = vcov(r)) }
  list(red = red, pct = quantile(x, PCT, na.rm = TRUE), n = sum(!is.na(fitted(m))), dispersion = summary(m)$dispersion)
}

# Stage 1 for all regions x groups (parallel)
run_stage1 <- function(panel, expo, groups = GROUPS, regions = sort(unique(panel$Region)), ...) {
  jobs <- CJ(region = regions, group = groups)
  res <- mclapply(seq_len(nrow(jobs)), function(i) {
    d <- panel[Region == jobs$region[i]][order(Date)]
    out <- tryCatch(fit_region(d, expo, jobs$group[i], ...), error = function(e) NULL)
    c(list(region = jobs$region[i], group = jobs$group[i]), out)
  }, mc.cores = NCORES)
  Filter(function(z) !is.null(z$red), res)
}

# Stage 2 for one group/window: pooled curve, MMT and optional coefficient draws
pool <- function(s1, g, window = "cum", vtype = "ns3", nsim = 0, boot_mmt = FALSE) {
  s <- Filter(function(z) z$group == g && !is.null(z$red[[window]]), s1)
  co <- t(sapply(s, function(z) z$red[[window]]$coef)); S <- lapply(s, function(z) z$red[[window]]$vcov)
  meta <- mixmeta(co ~ 1, S = S, method = "reml")
  tavg <- rowMeans(sapply(s, `[[`, "pct"))
  av <- switch(vtype,
    ns3 = list(fun = "ns", knots = tavg[sapply(c(.10, .75, .90), pidx)], Boundary.knots = range(tavg)),
    nsm = list(fun = "ns", knots = tavg[sapply(c(.10, .50, .90), pidx)], Boundary.knots = range(tavg)),
    bs2 = argvar_spec(tavg, "bs2"))
  b <- do.call(onebasis, c(list(x = tavg), av))
  lin <- drop(b %*% coef(meta))
  p <- list(coef = coef(meta), vcov = vcov(meta), basis = b, tavg = tavg, regions = sapply(s, `[[`, "region"),
            mmt = PCT[INRANGE[which.min(lin[INRANGE])]], i2 = summary(meta)$i2stat[1],
            blup = setNames(lapply(blup(meta, vcov = TRUE), function(z) z), sapply(s, `[[`, "region")))
  if (nsim > 0 || boot_mmt) p$sims <- mvrnorm(max(nsim, 1000), p$coef, p$vcov)
  if (boot_mmt) p$mmt_boot <- apply(p$sims, 1, function(cf) PCT[INRANGE[which.min((b %*% cf)[INRANGE])]])   # Tobias et al. 2017
  p
}

# log RR (and SE) at percentile q relative to percentile ref
lrr <- function(p, q, ref) { x <- p$basis[pidx(q), ] - p$basis[pidx(ref), ]
  c(est = sum(x * p$coef), se = sqrt(drop(t(x) %*% p$vcov %*% x))) }
rr_curve <- function(p, ref) { X <- sweep(p$basis, 2, p$basis[pidx(ref), ], "-")
  data.table(pct = PCT, temp = p$tavg, logrr = drop(X %*% p$coef), se = sqrt(rowSums((X %*% p$vcov) * X))) }
# ratio of RR between two pooled groups at percentile q (common reference); groups treated as independent
rr_ratio <- function(pa, pb, q, ref) { a <- lrr(pa, q, ref); b <- lrr(pb, q, ref); d <- unname(a[1] - b[1]); se <- unname(sqrt(a[2]^2 + b[2]^2))
  c(RRR = exp(d), lo = exp(d - 1.96 * se), hi = exp(d + 1.96 * se), p = 2 * pnorm(-abs(d / se)), se = se) }

# Whole-curve test between two disjoint groups: multivariate meta-regression coef ~ group (Gasparrini et al. 2012)
curve_test <- function(s1, pair, window = "cum") {
  s <- Filter(function(z) z$group %in% pair && !is.null(z$red[[window]]), s1)
  co <- t(sapply(s, function(z) z$red[[window]]$coef)); S <- lapply(s, function(z) z$red[[window]]$vcov)
  md <- data.frame(region = sapply(s, `[[`, "region"), grp = factor(sapply(s, `[[`, "group"), levels = rev(pair)))
  m <- tryCatch(mixmeta(co ~ grp, S = S, random = ~ 1 | region/grp, data = md, method = "reml"), error = function(e) NULL)
  if (is.null(m)) return(c(wald = NA, df = NA, p = NA))
  k <- grep("grp", names(coef(m))); b <- coef(m)[k]; V <- vcov(m)[k, k]; w <- drop(t(b) %*% solve(V) %*% b)
  c(wald = w, df = length(k), p = pchisq(w, length(k), lower.tail = FALSE))
}
