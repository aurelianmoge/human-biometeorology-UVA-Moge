# =============================================================================
# 01_null_simulation.R - what do the group-comparison methods report when NO group differs?
# Truth: every pseudo-group has exactly the all-cause temperature-mortality curve.  For each
# region, the all-cause stage-1 model gives fitted daily means mu_t; each replicate draws
#   total_t ~ NegBin(mu_t, theta)  (theta from the quasi-Poisson dispersion), then splits it
#   binomially / multinomially with the REAL margin shares (Black 20.2%, sex 50/50, age
#   9.7 / 35.5 / 54.8%), so pseudo-groups share day-level shocks like the real margins and
#   differ only in size.  The full pipeline (stage 1, mixmeta, MMT, curve test, RRR,
#   significance-based at-risk days) is applied unchanged.
# Long run (~100 min for 300 replicates on 8 cores); run_all.R runs it only with --with-simulation.
# Usage: Rscript "2. Scripts/07_methods_supplement/01_null_simulation.R" [replicates = 300]
# Output: tables/07_methods_supplement/null_groups.csv, null_pairs.csv
# =============================================================================
source("2. Scripts/00_setup/config.R"); source("2. Scripts/00_setup/dlnm_tools.R")
STAGE <- "07_methods_supplement"
args <- commandArgs(TRUE); NREP <- if (length(args)) as.integer(args[1]) else 300L
RNGkind("L'Ecuyer-CMRG"); set.seed(SEED)
panel <- load_panel(); regions <- sort(unique(panel$Region))
truth <- rbindlist(lapply(regions, function(r) { d <- panel[Region == r][order(Date)]; x <- d$Tmin_f
  cb <- crossbasis(x, lag = MAX_LAG, argvar = argvar_spec(x, "ns3"), arglag = list(knots = logknots(MAX_LAG, 3)))
  m <- glm(d$All ~ cb + ns(time, df = DF_YEAR * 16) + dow, data = d, family = quasipoisson, na.action = "na.exclude")
  phi <- summary(m)$dispersion; mu <- fitted(m)
  data.table(Region = r, Date = d$Date, mu = mu, theta = if (phi > 1) mean(mu, na.rm = TRUE) / (phi - 1) else 1e6) }))
panel <- merge(panel, truth, by = c("Region", "Date")); setorder(panel, Region, Date)
true_mmt <- read_model("pooled_Tmin")[["All|cum"]]$mmt
SIMG <- c("sB", "sW", "sF", "sM", "sU50", "s5074", "s75")
SIZE <- c(sB = .202, sW = .798, sF = .5, sM = .5, sU50 = .097, s5074 = .355, s75 = .548)
SPAIRS <- list(race = c("sB", "sW"), sex = c("sF", "sM"), age_old_young = c("s75", "sU50"), age_old_mid = c("s75", "s5074"))
w_pct <- c(diff(PCT), 0)
sig_share <- function(p) { rc <- rr_curve(p, p$mmt); lo <- rc$logrr - 1.96 * rc$se
  c(sig_cold = sum(w_pct[lo > 0 & PCT < p$mmt]), sig_heat = sum(w_pct[lo > 0 & PCT > p$mmt])) }
one_rep <- function(k) {
  d <- copy(panel[, .(Region, Date, year, time, dow, Tmin_f, mu, theta)]); ok <- !is.na(d$mu)
  tot <- rep(0L, nrow(d)); tot[ok] <- rnbinom(sum(ok), mu = d$mu[ok], size = d$theta[ok])
  d[, sB := rbinom(.N, tot, .202)][, sW := tot - sB][, sF := rbinom(.N, tot, .5)][, sM := tot - sF]
  d[, sU50 := rbinom(.N, tot, .097)][, s5074 := rbinom(.N, tot - sU50, .355 / (1 - .097))][, s75 := tot - sU50 - s5074]
  s1 <- run_stage1(d, "Tmin_f", groups = SIMG, windows = list(cum = c(0, NA)))
  pl <- setNames(lapply(SIMG, function(g) pool(s1, g, "cum")), SIMG)
  grp <- rbindlist(lapply(SIMG, function(g) { p <- pl[[g]]; ss <- sig_share(p); c1 <- lrr(p, .01, true_mmt); c99 <- lrr(p, .99, true_mmt)
    data.table(rep = k, group = g, size = SIZE[[g]], mmt = p$mmt, hot_days = 1 - p$mmt, lrr_p1 = c1[[1]], se_p1 = c1[[2]], lrr_p99 = c99[[1]], se_p99 = c99[[2]],
               own_p1 = lrr(p, .01, p$mmt)[[1]], own_p99 = lrr(p, .99, p$mmt)[[1]], sig_cold = ss[["sig_cold"]], sig_heat = ss[["sig_heat"]]) }))
  prs <- rbindlist(lapply(names(SPAIRS), function(nm) { pr <- SPAIRS[[nm]]; ct <- curve_test(s1, pr); a <- grp[group == pr[1]]; b <- grp[group == pr[2]]
    z1 <- (a$lrr_p1 - b$lrr_p1) / sqrt(a$se_p1^2 + b$se_p1^2); z99 <- (a$lrr_p99 - b$lrr_p99) / sqrt(a$se_p99^2 + b$se_p99^2)
    data.table(rep = k, pair = nm, p_curve = ct[["p"]], rrr_p1 = exp(a$lrr_p1 - b$lrr_p1), p_p1 = 2 * pnorm(-abs(z1)), rrr_p99 = exp(a$lrr_p99 - b$lrr_p99),
               p_p99 = 2 * pnorm(-abs(z99)), own_rrr_p1 = exp(a$own_p1 - b$own_p1), own_rrr_p99 = exp(a$own_p99 - b$own_p99),
               dmmt = a$mmt - b$mmt, dhot = a$hot_days - b$hot_days, dsig_cold = a$sig_cold - b$sig_cold) }))
  list(grp = grp, prs = prs)
}
t0 <- Sys.time(); res <- vector("list", NREP)
for (k in seq_len(NREP)) { res[[k]] <- tryCatch(one_rep(k), error = function(e) { message("replicate ", k, ": ", conditionMessage(e)); NULL })
  if (k %% 10 == 0) cat(k, "replicates,", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n") }
res <- Filter(Negate(is.null), res)
write_tab(rbindlist(lapply(res, `[[`, "grp")), STAGE, "null_groups"); write_tab(rbindlist(lapply(res, `[[`, "prs")), STAGE, "null_pairs")
cat(length(res), "replicates written\n")
