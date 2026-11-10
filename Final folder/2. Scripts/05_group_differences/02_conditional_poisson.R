# =============================================================================
# 02_conditional_poisson.R - interpretable slopes and group-equality tests.
# Time-stratified conditional Poisson (fixest::fepois), strata region x year x month x day
# of week (x group in stacked models): the count analogue of a time-stratified
# case-crossover (Armstrong et al. 2014).  SEs clustered by region x year-month; in stacked
# models the clusters are shared, so the covariance between groups' slopes is estimated.
# Terms (region-specific thresholds at the regional P90 of Tmin, close to the MMT):
#   cold   = max(0, P90 - mean Tmin lag 0-21)    per 1 C
#   heat   = max(0, mean Tmin lag 0-2 - P90)     per 1 C
#   DTR02  = mean DTR lag 0-2 / 5;  DTR21 = mean DTR lag 0-21 / 5   (per 5 C)
#   lowDTR = share of lag 3-21 days with DTR < regional P10, per 10 pp (joint sex test)
# Models: M1 all deaths; M2 sex; M3 race; M4 age (stacked, group-specific slopes);
#         J  joint 3-df test of equal cold / heat / lowDTR slopes for each pair.
# Output: tables/05_group_differences/cpoisson_*.csv
# =============================================================================
source("2. Scripts/00_setup/config.R"); source("2. Scripts/00_setup/dlnm_tools.R")
suppressPackageStartupMessages(library(fixest))
STAGE <- "05_group_differences"
panel <- load_panel(); setorder(panel, Region, Date)
ma <- function(x, k) frollmean(x, k, align = "right")
panel[, `:=`(t21 = ma(Tmin_f, 22), t02 = ma(Tmin_f, 3), d21 = ma(DTR_f, 22) / 5, d02 = ma(DTR_f, 3) / 5,
             low = as.numeric(DTR_f < quantile(DTR_f, .10, na.rm = TRUE)), p90 = quantile(Tmin_f, .90, na.rm = TRUE)), by = Region]
panel[, low321 := shift(frollmean(low, 19, align = "right"), 3) * 10, by = Region]
panel[, `:=`(cold = pmax(0, p90 - t21), heat = pmax(0, t02 - p90), stratum = paste(Region, year, month, dow), clu = paste(Region, year, month))]
pct <- function(b) 100 * (exp(b) - 1)
tidy <- function(m, model) { ct <- coeftable(m)
  data.table(model = model, term = rownames(ct), b = ct[, 1], se = ct[, 2], p = ct[, 4], effect_pct = pct(ct[, 1]),
             lo_pct = pct(ct[, 1] - 1.96 * ct[, 2]), hi_pct = pct(ct[, 1] + 1.96 * ct[, 2]), n = nobs(m)) }
stack <- function(groups, terms) { l <- melt(panel[, c("stratum", "clu", terms, groups), with = FALSE], measure.vars = groups, variable.name = "grp", value.name = "y")[!is.na(y)]
  l[, grp := factor(grp, groups)]
  fepois(as.formula(paste("y ~", paste0("i(grp, ", terms, ")", collapse = " + "), "| stratum^grp")), l, cluster = ~clu) }
eq_test <- function(m, terms, groups) {                                        # Wald test of equal slopes across groups
  nm <- unlist(lapply(terms, function(t) paste0("grp::", groups, ":", t))); b <- coef(m)[nm]; V <- vcov(m)[nm, nm]
  k <- length(groups); blk <- if (k == 2) matrix(c(1, -1), 1) else cbind(1, -diag(k - 1))
  C <- kronecker(diag(length(terms)), blk); d <- C %*% b; w <- drop(t(d) %*% solve(C %*% V %*% t(C)) %*% d)
  data.table(terms = paste(terms, collapse = "+"), groups = paste(groups, collapse = "/"), wald = w, df = nrow(C), p = pchisq(w, nrow(C), lower.tail = FALSE)) }

# ---- M1: all deaths ---------------------------------------------------------------------------
a <- panel[!is.na(All)]
m1 <- fepois(All ~ cold + heat + d21 + d02 | stratum, a, cluster = ~clu)
# ---- M2-M4: stacked margins ----------------------------------------------------------------------
T4 <- c("cold", "heat", "d21", "d02")
M <- list(sex = c("Female", "Male"), race = c("Black", "White"), age = c("Under50", "Age50_74", "Age75plus"))
fits <- lapply(M, stack, terms = T4)
tab <- rbindlist(c(list(tidy(m1, "M1_all")), Map(tidy, fits, paste0("M", 2:4, "_", names(M)))))
tests <- rbindlist(lapply(names(M), function(k) rbindlist(lapply(T4, function(t) eq_test(fits[[k]], t, M[[k]])))))
# ---- covariance between the Female and Male slopes (why overlapping CIs can differ significantly) --------
V <- vcov(fits$sex); cov_fm <- rbindlist(lapply(c("cold", "heat"), function(t) { f <- paste0("grp::Female:", t); m <- paste0("grp::Male:", t)
  d <- coef(fits$sex)[f] - coef(fits$sex)[m]; s_ind <- sqrt(V[f, f] + V[m, m]); s <- sqrt(V[f, f] + V[m, m] - 2 * V[f, m])
  data.table(term = t, corr = V[f, m] / sqrt(V[f, f] * V[m, m]), diff = d, se_if_independent = s_ind, p_if_independent = 2 * pnorm(-abs(d / s_ind)),
             se_with_covariance = s, p_with_covariance = 2 * pnorm(-abs(d / s))) }))
# ---- J: joint tests of cold / heat / lowDTR slopes ------------------------------------------------------
TJ <- c("cold", "heat", "low321", "d02")
joint <- rbindlist(lapply(list(c("Female", "Male"), c("Black", "White"), c("Age75plus", "Age50_74")), function(g) {
  m <- stack(g, TJ); rbind(eq_test(m, c("cold", "heat", "low321"), g), tidy(m, paste(g, collapse = "/"))[, .(terms = term, groups = model, wald = effect_pct, df = NA, p)], fill = TRUE) }))
write_tab(tab, STAGE, "cpoisson_slopes"); write_tab(tests, STAGE, "cpoisson_equality_tests")
write_tab(cov_fm, STAGE, "cpoisson_sex_covariance"); write_tab(joint, STAGE, "cpoisson_joint_tests")
options(width = 200)
print(tab[, .(model, term, effect = round(effect_pct, 2), ci = sprintf("[%.2f, %.2f]", lo_pct, hi_pct), p = signif(p, 2))])
print(tests); print(joint[!is.na(df)])
