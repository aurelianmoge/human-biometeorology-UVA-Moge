# =============================================================================
# 01_days_at_risk.R - Cao et al. (2026) comparison: days and person-days "at risk" counted
# with each margin's OWN pooled curve vs the OVERALL all-cause curve.
#   rules, Tmin : colder / warmer than the MMT; RR >= 1.03, 1.05, 1.07, 1.10 (cold and heat side)
#   rules, DTR  : below / above the curve minimum; RR >= 1.03 (minimum poorly identified)
#   share of region-days at risk (no denominators needed);
#   person-days PD_g = sum_{region, year} population_g x at-risk days (regional populations,
#     script 01/06); Delta_g = (PD_own - PD_common) / PD_common  (Cao's under-estimation);
#     partition totals E_specific = sum_g PD_g(own) vs E_common;
#   uncertainty: 1,000 joint coefficient draws; noise-only benchmark: the margin's curve drawn
#     around the ALL-CAUSE curve with the margin's own covariance (identical true curves).
# Input : models/pooled_Tmin.rds, pooled_DTR.rds; derived/regional_populations.csv
# Output: tables/06_days_at_risk/days_at_risk.csv, partition_totals.csv
# =============================================================================
source("2. Scripts/00_setup/config.R"); source("2. Scripts/00_setup/dlnm_tools.R")
STAGE <- "06_days_at_risk"; set.seed(SEED); NS <- 1000
panel <- load_panel()
rp <- fread(file.path(DERIVED, "regional_populations.csv"))
PLT <- read_model("pooled_Tmin"); Pt <- setNames(lapply(GROUPS, function(g) PLT[[paste(g, "cum", sep = "|")]]), GROUPS)
Pd <- read_model("pooled_DTR")

days <- panel[!is.na(Tmin_f) & !is.na(DTR_f), .(Region, year, Tmin_f, DTR_f)]
days[, `:=`(ip_T = sapply(ecdf(Tmin_f)(Tmin_f), pidx), ip_D = sapply(ecdf(DTR_f)(DTR_f), pidx)), by = Region]
w_days <- function(ip) { v <- numeric(length(PCT)); a <- days[, .N, by = .(i = get(ip))]; v[a$i] <- a$N; v }
w_pd <- function(ip, g) { a <- merge(days[, .(n = .N), by = .(Region, year, i = get(ip))], rp[Group == g], by.x = c("Region", "year"), by.y = c("Region", "Year"))[, .(pd = sum(n * Population)), by = i]
  v <- numeric(length(PCT)); v[a$i] <- a$pd; v }
TH <- c(1.03, 1.05, 1.07, 1.10)
flags <- function(basis, cf, expo) {
  lin <- drop(basis %*% cf); r <- INRANGE[which.min(lin[INRANGE])]; l <- lin - lin[r]; lo <- seq_along(PCT) < r; hi <- seq_along(PCT) > r
  if (expo == "Tmin") cbind(MMT_cold = lo, MMT_heat = hi, sapply(setNames(TH, paste0("RR", TH, "_cold")), function(t) lo & l >= log(t)),
                            sapply(setNames(TH, paste0("RR", TH, "_heat")), function(t) hi & l >= log(t)))
  else cbind(MIN_low = lo, MIN_high = hi, RR1.03_low = lo & l >= log(1.03), RR1.03_high = hi & l >= log(1.03)) }

res <- list()
for (expo in c("Tmin", "DTR")) {
  P <- if (expo == "Tmin") Pt else Pd; ip <- if (expo == "Tmin") "ip_T" else "ip_D"; dw <- w_days(ip)
  Fc <- flags(P$All$basis, P$All$coef, expo); Fc_d <- lapply(seq_len(NS), function(s) flags(P$All$basis, P$All$sims[s, ], expo))
  for (g in setdiff(GROUPS, "All")) { p <- P[[g]]; pw <- w_pd(ip, g); Fo <- flags(p$basis, p$coef, expo)
    own_d <- sapply(seq_len(NS), function(s) colSums(flags(p$basis, p$sims[s, ], expo) * pw))
    com_d <- sapply(seq_len(NS), function(s) colSums(Fc_d[[s]] * pw))
    nul <- mvrnorm(NS, P$All$coef, p$vcov); nul_d <- sapply(seq_len(NS), function(s) colSums(flags(P$All$basis, nul[s, ], expo) * pw))
    sh_d <- sapply(seq_len(NS), function(s) colSums(flags(p$basis, p$sims[s, ], expo) * dw)) / sum(dw)
    PDo <- colSums(Fo * pw); PDc <- colSums(Fc * pw); rel <- (own_d - com_d) / pmax(com_d, 1); reln <- (nul_d - PDc) / pmax(PDc, 1)
    res[[length(res) + 1]] <- data.table(exposure = expo, group = g, rule = colnames(Fo),
      share_own = colSums(Fo * dw) / sum(dw), share_own_lo = apply(sh_d, 1, quantile, .025), share_own_hi = apply(sh_d, 1, quantile, .975),
      share_common = colSums(Fc * dw) / sum(dw), PD_own_millions = PDo / 1e6, PD_common_millions = PDc / 1e6,
      delta_pct = 100 * (PDo - PDc) / pmax(PDc, 1), delta_lo = 100 * apply(rel, 1, quantile, .025), delta_hi = 100 * apply(rel, 1, quantile, .975),
      noise_lo = 100 * apply(reln, 1, quantile, .025), noise_hi = 100 * apply(reln, 1, quantile, .975),
      p_vs_noise = sapply(seq_len(nrow(reln)), function(j) { d <- (PDo[j] - PDc[j]) / max(PDc[j], 1); m <- median(reln[j, ]); mean(abs(reln[j, ] - m) >= abs(d - m)) })) }
  cat(expo, "done\n") }
res <- rbindlist(res); write_tab(res, STAGE, "days_at_risk")
part <- list(Sex = c("Female", "Male"), Race = c("Black", "White"), Age = c("Under50", "Age50_74", "Age75plus"))
tot <- rbindlist(lapply(names(part), function(nm) res[group %in% part[[nm]], .(E_specific_millions = sum(PD_own_millions), E_common_millions = sum(PD_common_millions)), by = .(exposure, rule)][, partition := nm]))
tot[, delta_pct := fifelse(E_common_millions > 0, 100 * (E_specific_millions - E_common_millions) / E_common_millions, NA_real_)]
write_tab(tot, STAGE, "partition_totals")
print(res[exposure == "Tmin" & rule %in% c("MMT_heat", "MMT_cold", "RR1.05_cold"), .(rule, group, own = round(100 * share_own, 1), common = round(100 * share_common, 1),
          delta = round(delta_pct), p = round(p_vs_noise, 2))])
