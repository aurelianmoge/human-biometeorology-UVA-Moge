# =============================================================================
# 03_attributable_burden.R - absolute inequality: cold- and heat-attributable deaths per
# 100,000 person-years, by margin, Virginia 2005-2020.
#   AN_g = sum over regions and days of (1 - 1/RR_g,r(x_t)) x mean deaths over t..t+21
#          (forward attributable number, Gasparrini & Leone 2014), with the region-specific
#          BLUP curve of margin g centred at the COMMON all-cause MMT (same counterfactual
#          for everyone) and, as sensitivity, at the margin's own MMT;
#   rate = AN_g / statewide person-years_g (the 12 regions cover the whole state);
#   decomposition: attributable-rate ratio = baseline-rate ratio x attributable-fraction ratio;
#   race, indirectly age-standardised (Anderson & Rosenberg 1998): expected_k = sum_a rate_a PY_k,a
#     with Census race x age person-years, for all-cause and cold-attributable deaths.
# Uncertainty: 300 draws of the regional BLUP coefficients.
# Output: tables/05_group_differences/burden_*.csv
# =============================================================================
source("2. Scripts/00_setup/config.R"); source("2. Scripts/00_setup/dlnm_tools.R")
STAGE <- "05_group_differences"; set.seed(SEED); NSIM <- 300
panel <- load_panel(); PL <- read_model("pooled_Tmin")
spop <- fread(file.path(DERIVED, "statewide_population_margins.csv"))[Year %in% YEARS]
cpop <- fread(file.path(DERIVED, "county_joint_populations.csv"))[Year %in% YEARS]
mmt_all <- PL[["All|cum"]]$mmt
rb <- function(x, at) onebasis(at, fun = "ns", knots = quantile(x, c(.1, .75, .9), na.rm = TRUE), Boundary.knots = range(x, na.rm = TRUE))
reg <- lapply(split(panel, by = "Region"), function(d) { d <- d[order(Date)]; list(d = d, B = rb(d$Tmin_f, d$Tmin_f)) })
an <- function(g, cf, mmt_p) { tot <- c(cold = 0, heat = 0, deaths = 0)
  for (r in names(cf)) { d <- reg[[r]]$d; x <- d$Tmin_f; y <- d[[g]]; cen <- quantile(x, mmt_p)
    l <- drop(sweep(reg[[r]]$B, 2, rb(x, cen), "-") %*% cf[[r]]); a <- (1 - exp(-l)) * frollmean(y, MAX_LAG + 1, align = "left", na.rm = TRUE)
    tot <- tot + c(sum(a[x < cen], na.rm = TRUE), sum(a[x > cen], na.rm = TRUE), sum(y, na.rm = TRUE)) }
  tot }
py <- spop[, .(person_years = sum(Population)), by = Group]
out <- list(); draws <- list()
for (g in GROUPS) { p <- PL[[paste(g, "cum", sep = "|")]]
  bl <- lapply(p$blup, `[[`, "blup"); bv <- lapply(p$blup, `[[`, "vcov")
  for (cen in c("common", "own")) { mp <- if (cen == "common") mmt_all else p$mmt; pt <- an(g, bl, mp)
    sims <- do.call(rbind, mclapply(seq_len(NSIM), function(s) { set.seed(SEED + s); an(g, lapply(setNames(names(bl), names(bl)), function(r) mvrnorm(1, bl[[r]], bv[[r]])), mp) }, mc.cores = NCORES))
    pyg <- py[Group == g, person_years]
    out[[length(out) + 1]] <- data.table(group = g, centre = cen, deaths = pt[["deaths"]], person_years = pyg, baseline_rate = 1e5 * pt[["deaths"]] / pyg,
      AF_cold = pt[["cold"]] / pt[["deaths"]], AF_heat = pt[["heat"]] / pt[["deaths"]],
      cold_rate = 1e5 * pt[["cold"]] / pyg, cold_lo = 1e5 * quantile(sims[, "cold"], .025) / pyg, cold_hi = 1e5 * quantile(sims[, "cold"], .975) / pyg,
      heat_rate = 1e5 * pt[["heat"]] / pyg, heat_lo = 1e5 * quantile(sims[, "heat"], .025) / pyg, heat_hi = 1e5 * quantile(sims[, "heat"], .975) / pyg)
    draws[[paste(g, cen)]] <- data.table(group = g, centre = cen, sim = seq_len(NSIM), cold = sims[, "cold"], heat = sims[, "heat"], person_years = pyg) }
  cat(g, "done\n") }
out <- rbindlist(out); draws <- rbindlist(draws); write_tab(out, STAGE, "burden_rates")
dec <- rbindlist(lapply(PAIRS, function(pr) { a <- out[group == pr[1] & centre == "common"]; b <- out[group == pr[2] & centre == "common"]
  da <- draws[group == pr[1] & centre == "common"]; db <- draws[group == pr[2] & centre == "common"]; s <- (da$cold / da$person_years) / (db$cold / db$person_years)
  data.table(pair = paste(pr, collapse = "/"), baseline_rate_ratio = a$baseline_rate / b$baseline_rate, AF_ratio = a$AF_cold / b$AF_cold,
             cold_rate_ratio = a$cold_rate / b$cold_rate, lo = quantile(s, .025), hi = quantile(s, .975), cold_rate_difference_per100k = a$cold_rate - b$cold_rate) }))
write_tab(dec, STAGE, "burden_rate_ratio_decomposition")

# ---- race, indirectly age-standardised ------------------------------------------------------------
ages <- c("Under50", "Age50_74", "Age75plus")
PY <- cpop[Race %in% c("Black", "White"), .(PY = sum(Population)), by = .(Race, Age)]
br <- out[centre == "common"]
ex <- merge(PY, br[group %in% ages, .(Age = group, all_rate = deaths / person_years, cold_rate = cold_rate / 1e5)], by = "Age")[, .(exp_all = sum(PY * all_rate), exp_cold = sum(PY * cold_rate)), by = Race]
ob <- br[group %in% c("Black", "White"), .(Race = group, obs_all = deaths, obs_cold = cold_rate * person_years / 1e5)]
smr <- merge(ob, ex, by = "Race")[, `:=`(SMR_all = obs_all / exp_all, SMR_cold = obs_cold / exp_cold)]
dc <- draws[centre == "common"]
sim <- sapply(seq_len(NSIM), function(s) { e <- merge(PY, dc[sim == s & group %in% ages, .(Age = group, cr = cold / person_years)], by = "Age")[, .(e = sum(PY * cr)), by = Race]
  o <- dc[sim == s & group %in% c("Black", "White"), .(Race = group, o = cold)]; m <- merge(o, e); m[Race == "Black", o / e] / m[Race == "White", o / e] })
std <- data.table(measure = c("all-cause SMR ratio (Black/White)", "cold-attributable SMR ratio (Black/White)"),
                  estimate = c(smr[Race == "Black", SMR_all] / smr[Race == "White", SMR_all], smr[Race == "Black", SMR_cold] / smr[Race == "White", SMR_cold]),
                  lo = c(NA, quantile(sim, .025)), hi = c(NA, quantile(sim, .975)))
write_tab(smr, STAGE, "burden_race_SMR"); write_tab(std, STAGE, "burden_race_standardised_ratio")
print(out[centre == "common", .(group, AF_cold = round(AF_cold, 4), cold_rate = round(cold_rate, 1), cold_lo = round(cold_lo, 1), cold_hi = round(cold_hi, 1))]); print(std)
