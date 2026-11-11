# =============================================================================
# 02_paper_assets.R - tables (Markdown) and figures of the paper, built from saved outputs only.
#   T1  sample and exposure by group             (02_descriptives/03_exposure_by_group)
#   T2  temperature-mortality by group: exact DeGuzman model, direct test, regional model
#   T2b conditional Poisson cold / heat slopes
#   T3  DTR contrasts (regional model, variants); T3b DTR, statewide design
#   T4  decomposition of attributable-rate gaps; T4b with the nearest-station crosswalk
#   T5  regions; T5b poverty meta-regression; T5c high- vs low-poverty tercile decomposition
#   T6a-g thresholds (06_thresholds): curve thresholds, equal size, threshold model, risk levels, tests,
#       days at risk, known-null comparison
#   F1 exposure by group, F3 DTR curves, F4 decomposition, F5 regions, F6 risk levels, F7 missed deaths
#   (+ copies of 01_main_figures)
# Output: 3. Outputs/paper/tables_md/*.md, 3. Outputs/paper/figures/*.png
# =============================================================================
source("2. Scripts/00_setup/config.R"); source("2. Scripts/00_setup/dlnm_tools.R")
tab <- function(stage, name) fread(out_path("tables", stage, paste0(name, ".csv")))
P_DIR <- file.path(OUT, "paper"); for (d in c("tables_md", "figures")) dir.create(file.path(P_DIR, d), recursive = TRUE, showWarnings = FALSE)
md <- function(d, name) { d <- copy(as.data.table(d)); for (j in names(d)) if (is.numeric(d[[j]])) set(d, j = j, value = formatC(d[[j]], format = "fg", digits = 3))
  writeLines(c(paste("|", paste(names(d), collapse = " | "), "|"), paste("|", paste(rep("---", ncol(d)), collapse = " | "), "|"),
               apply(d, 1, function(r) paste("|", paste(r, collapse = " | "), "|"))), file.path(P_DIR, "tables_md", paste0(name, ".md"))) }
fg <- function(g, n, w = 10, h = 4.5) ggsave(file.path(P_DIR, "figures", paste0(n, ".png")), g, width = w, height = h, dpi = 170)
ci <- function(e, l, h, d = 2) sprintf(paste0("%.", d, "f (%.", d, "f–%.", d, "f)"), e, l, h)
GL <- GROUP_LABEL
PLB <- c("Female/Male" = "Female vs Male", "Black/White" = "Black vs White", "Age75plus/Age50_74" = "75+ vs 50–74",
         "Age75plus/Under50" = "75+ vs Under 50", "Age50_74/Under50" = "50–74 vs Under 50")
FAM <- c(All = "All", Female = "Sex", Male = "Sex", Black = "Race", White = "Race", Under50 = "Age", Age50_74 = "Age", Age75plus = "Age")

# ---- T1, F1: sample and exposure ------------------------------------------------------------------------------
t1d <- tab("02_descriptives", "T1_deaths_rates_statewide"); E <- tab("02_descriptives", "exposure_by_group")
t1 <- merge(t1d[, .(Group, Deaths = deaths, `Share of deaths` = sprintf("%.1f%%", 100 * share_of_deaths), `Crude death rate per 100,000` = round(crude_rate_per100k))],
            E[weights == "residents, rebuilt crosswalk", .(Group, `Mean Tmin (°C)` = round(mean_Tmin, 2), `Frost days / yr` = round(frost_days, 1),
              `Days Tmin ≤ −10 °C / yr` = round(very_cold_days, 1), `Warm nights ≥ 22 °C / yr` = round(warm_nights, 1),
              `Days heat index ≥ 100 °F / yr` = round(heat_index_100F_days, 1), `Mean DTR (°C)` = round(mean_DTR, 2), `Days DTR ≥ 15 °C / yr` = round(high_DTR_days, 1))],
            by = "Group")[match(names(GL), Group)][, Group := GL[Group]]
md(t1, "T1_sample_exposure")
geo <- tab("02_descriptives", "residential_distribution"); pc <- function(v) sprintf("%.1f%%", 100 * v)
md(geo[, .(Region, Black = pc(Black), White = pc(White), Female = pc(Female), Male = pc(Male), `75+` = pc(Age75plus))], "T1b_residence")
EL <- c(frost_days = "Frost days (Tmin < 0 °C)", very_cold_days = "Very cold days (Tmin ≤ −10 °C)", warm_nights = "Warm nights (Tmin ≥ 22 °C)",
        heat_index_100F_days = "Heat index ≥ 100 °F", high_DTR_days = "DTR ≥ 15 °C")
ev <- melt(E[Group != "All"], id.vars = c("Group", "weights"), measure.vars = names(EL))[, `:=`(v = factor(EL[variable], EL), g = factor(GL[Group], rev(GL)))]
fg(ggplot(ev, aes(value, g, colour = Group, shape = weights)) + geom_point(size = 2.4, position = position_dodge(.6)) + facet_wrap(~v, scales = "free_x", nrow = 1) +
     scale_colour_manual(values = GROUP_COL, guide = "none") + scale_shape_manual(values = c(16, 1, 4), name = "weighted by") +
     labs(title = "Days per year experienced by each group, given where it lives",
          subtitle = "Regional station series weighted by each group's person-years (two county–region crosswalks) or by where its deaths occur.", x = "Days per year", y = NULL),
   "F1_exposure_by_group", 13, 4)

# ---- T2: temperature-mortality ---------------------------------------------------------------------------------------
EX <- tab("05_group_differences", "deguzman_exact_RR"); EC <- tab("05_group_differences", "deguzman_exact_contrasts")
PW <- tab("03_temperature", "pointwise_RRR")[exposure == "Tmin" & window == "cum" & centre == "MMT_all"]
t2 <- rbindlist(lapply(PAIRS[c(1, 2, 4, 3)], function(pr) { pn <- paste(pr, collapse = "/")
  a <- EX[group == pr[1] & temp == -10]; b <- EX[group == pr[2] & temp == -10]; k <- EC[pair == pn & temp == -10]
  q1 <- PW[pair == pn & pct == 1]; q99 <- PW[pair == pn & pct == 99]
  data.table(Comparison = PLB[pn], `Statewide model, RR at −10 °C vs 5 °C (A; B)` = paste(ci(a$RR, a$lo, a$hi), ";", ci(b$RR, b$lo, b$hi)),
             `Excess-risk ratio` = sprintf("%.1f×", k$excess_ratio), `Direct test: RR ratio (p)` = sprintf("%s, p = %.2f", ci(k$RRR, k$RRR_lo, k$RRR_hi), k$p_RRR),
             `Regional model, RR ratio at P1 (p)` = sprintf("%.2f, p = %.2f", q1$RRR, q1$p), `RR ratio at P99 (p)` = sprintf("%.2f, p = %.3f", q99$RRR, q99$p)) }))
md(t2, "T2_temperature_groups")
md(EX[group %in% c("Black", "White") & !is.na(published_RR), .(Group = group, `Tmin (°C)` = temp, `RR, this study` = ci(RR, lo, hi), `RR, DeGuzman et al. 2026` = published_RR)], "T2a_deguzman_replication")
cs <- tab("05_group_differences", "cpoisson_slopes")[grepl("^M[234]", model) & grepl(":(cold|heat)$", term)]
md(cs[, .(Group = GL[sub("grp::(\\w+):.*", "\\1", term)], Side = sub(".*:", "", term), `% change in deaths per °C (95% CI)` = ci(effect_pct, lo_pct, hi_pct))], "T2b_slopes")

# ---- T3, F3: DTR --------------------------------------------------------------------------------------------------------
dc <- tab("04_dtr", "dtr_contrasts")[pair %in% names(PLB)[-4] & pct %in% c(10, 90) & window == "cum" & variant %in% c("base", "humidity", "warm", "cold")]
md(dcast(dc[, .(pair = PLB[pair], pct = paste0("P", pct), variant, v = sprintf("%.2f (p = %.3f)", RRR, p_point))], pair + pct ~ variant, value.var = "v")[
  , .(Comparison = pair, `DTR percentile` = pct, `All months` = base, `+ humidity` = humidity, `May–Sep` = warm, `Nov–Mar` = cold)], "T3_dtr")
pdd <- tab("05_group_differences", "published_design_contrasts")[exposure == "DTR" & point %in% c("P1", "P99") & df_year == 3]
md(pdd[, .(Comparison = PLB[pair], Lag = fifelse(window == "cum", "0–21", "0–2"), Point = point, `RR A` = round(RR_a, 3), `RR B` = round(RR_b, 3),
           `RR ratio (p)` = sprintf("%.3f (p = %.3f)", RRR, p_RRR))], "T3b_dtr_statewide")
PD <- read_model("pooled_DTR")
cur <- rbindlist(lapply(GROUPS, function(g) rr_curve(PD[[g]], .5)[pct >= .01 & pct <= .99][, group := g]))
allc <- cur[group == "All"][, k := 1][data.table(k = 1, fam = c("Sex", "Race", "Age")), on = "k", allow.cartesian = TRUE]
cg <- cur[group != "All"][, fam := FAM[group]]
lab <- unique(cg[, .(group, fam)])[, `:=`(x = 2, y = 1.24 - .035 * (seq_len(.N) - 1)), by = fam]   # key stacked in each panel corner
fg(ggplot(cg, aes(100 * pct, exp(logrr), colour = group, fill = group)) + geom_hline(yintercept = 1, colour = COL[["grey"]]) +
     geom_line(data = allc, aes(100 * pct, exp(logrr)), inherit.aes = FALSE, colour = COL[["grey"]], linetype = 2) +
     geom_ribbon(aes(ymin = exp(logrr - 1.96 * se), ymax = exp(logrr + 1.96 * se)), alpha = .12, colour = NA) + geom_line(linewidth = .8) +
     geom_text(data = lab, aes(x = x, y = y, label = GL[group]), hjust = 0, size = 3.2, fontface = "bold", inherit.aes = FALSE, colour = GROUP_COL[lab$group]) +
     facet_wrap(~factor(fam, c("Sex", "Race", "Age"))) + scale_colour_manual(values = GROUP_COL, guide = "none") + scale_fill_manual(values = GROUP_COL, guide = "none") +
     labs(title = "DTR–mortality curves by group (lag 0–21), net of mean temperature",
       subtitle = "Two-stage DLNM pooled over 12 regions; reference = regional median DTR; dashed grey = all deaths.", x = "Regional DTR percentile", y = "Relative risk"),
   "F3_dtr_curves", 10.5, 4)

# ---- T4, F4: decomposition ---------------------------------------------------------------------------------------------------
DEC <- tab("05_group_differences", "decomposition"); CL <- c(gap = "Total gap", geography = "Geography", frailty = "Frailty", sensitivity = "Sensitivity")
wide <- function(d) dcast(d[, .(side, pair, component, v = sprintf("%.2f (%.2f to %.2f)", estimate, lo, hi))], side + pair ~ component, value.var = "v")
w4 <- wide(DEC[crosswalk == "rebuilt"])
md(w4[, .(Outcome = fifelse(side == "cold", "Cold", "Heat (lag 0–2)"), Comparison = PLB[pair], `Rate A` = rate_A, `Rate B` = rate_B, `Gap (A − B)` = gap,
          Geography = geography, Frailty = frailty, Sensitivity = sensitivity)][order(Outcome)], "T4_decomposition")
w4b <- wide(DEC[crosswalk == "nearest-station"])
md(w4b[, .(Outcome = fifelse(side == "cold", "Cold", "Heat (lag 0–2)"), Comparison = PLB[pair], Geography = geography, Frailty = frailty, Sensitivity = sensitivity)][order(Outcome)],
   "T4b_decomposition_nearest")
dp <- DEC[crosswalk == "rebuilt" & component %in% names(CL) & pair != "Age75plus/Under50"][, `:=`(comp = factor(CL[component], CL), pl = factor(PLB[pair], rev(PLB)),
  out = fifelse(side == "cold", "Cold-attributable deaths per 100,000", "Heat-attributable deaths per 100,000 (lag 0–2)"))]
fg(ggplot(dp, aes(estimate, pl, colour = comp)) + geom_vline(xintercept = 0, colour = COL[["grey"]]) +
     geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(.65), size = .3) + facet_wrap(~out, scales = "free_x") +
     scale_colour_manual(values = c("Total gap" = "#0b0b0b", Geography = COL[["aqua"]], Frailty = COL[["orange"]], Sensitivity = COL[["blue"]]), name = NULL) +
     labs(title = "Gaps in temperature-attributable deaths: geography, frailty or sensitivity?",
          subtitle = "Shapley decomposition of the gap (first minus second group), per 100,000 person-years. 95% intervals from 200 draws of the regional curves.",
          x = "Contribution to the gap (deaths per 100,000 person-years)", y = NULL), "F4_decomposition", 12, 4.6)

# ---- T5, F5: regions --------------------------------------------------------------------------------------------------------------
RG <- tab("05_group_differences", "decomposition_regions")
md(RG[order(poverty), .(Region, `Poverty (%)` = round(poverty, 1), `Baseline deaths per 100,000` = round(baseline_rate), `Cold-attributable fraction (%)` = ci(cold_af, cold_af - 1.96 * se_cold_af, cold_af + 1.96 * se_cold_af),
                        `Cold-attributable deaths per 100,000` = ci(cold_rate, cold_rate - 1.96 * se_cold_rate, cold_rate + 1.96 * se_cold_rate, 1))], "T5_regions")
MR <- tab("05_group_differences", "regional_poverty_metaregression")
md(MR[, .(Outcome = outcome, `Change per 10 points of poverty (95% CI)` = ci(slope_per_10pp_poverty, lo, hi), p = signif(p, 2), `I² (%)` = round(I2), Spearman = round(spearman, 2))], "T5b_poverty_metaregression")
TD <- tab("05_group_differences", "regional_tercile_decomposition")
md(dcast(TD[, .(side, component, v = ci(estimate, lo, hi, 1))], side ~ component, value.var = "v")[, .(Outcome = fifelse(side == "cold", "Cold", "Heat (lag 0–2)"),
   `High-poverty regions` = rate_high, `Low-poverty regions` = rate_low, Gap = gap, Frailty = frailty, `Temperature risk` = temperature_risk)], "T5c_tercile_decomposition")
fg(ggplot(RG, aes(poverty, baseline_rate)) + geom_point(aes(size = cold_af), colour = COL[["blue"]], alpha = .8) + geom_text(aes(label = Region), vjust = -1.1, size = 3) +
     scale_size_continuous(name = "cold-attributable\nfraction (%)", range = c(2, 6)) +
     labs(title = "Poorer regions: higher baseline mortality, similar cold-attributable fraction", x = "Population-weighted poverty rate (%)",
          y = "All-cause deaths per 100,000 person-years"), "F5_regions", 8.5, 4.4)

# ---- T6, F6, F7: thresholds and days at risk (06_thresholds) --------------------------------------------------------------------
G7 <- setdiff(GROUPS, "All"); ord <- function(d) d[order(match(group, GROUPS))]
CT <- tab("06_thresholds", "curve_thresholds")
fpc <- function(p, l, h, C) fifelse(is.na(p), sprintf("none (draws: P%.0f–P%.0f)", 100 * l, 100 * h), sprintf("P%.0f (P%.0f–P%.0f), %.1f °C", 100 * p, 100 * l, 100 * h, C))
md(ord(dcast(CT[, .(group, v = fpc(pct, lo, hi, degC), key = paste(side, threshold))], group ~ key, value.var = "v"))[
  , .(Group = GL[group], `Cold MMT` = `cold mmt`, `Cold significance threshold` = `cold sig`, `Heat MMT` = `heat mmt`, `Heat significance threshold` = `heat sig`)], "T6a_curve_thresholds")
ES <- tab("06_thresholds", "equal_size_thresholds")
md(ES[source != "All", .(`Larger group` = GL[source], `Thinned to` = GL[thinned_to], Threshold = c(cold_mmt = "Cold MMT", cold_sig = "Cold significance", heat_mmt = "Heat MMT", heat_sig = "Heat significance")[threshold],
         `Full size` = fifelse(is.na(full_size), "none", sprintf("P%.0f", 100 * full_size)), `Equal size, median (95% range)` = sprintf("P%.0f (P%.0f–P%.0f)", 100 * equal_size_median, 100 * equal_size_lo, 100 * equal_size_hi),
         `Replicates with a threshold` = sprintf("%.0f%%", 100 * share_with_threshold), `Smaller group` = fifelse(is.na(smaller_group_actual), "none", sprintf("P%.0f", 100 * smaller_group_actual)), `Size-matched p` = round(p_size_matched, 2))], "T6b_equal_size")
TG <- tab("06_thresholds", "tm_thresholds"); HS <- TG[p_heat_slope < .05, group]                  # groups with a significant heat effect
md(ord(TG)[, .(Group = GL[group], `Cold: risk starts below` = sprintf("%.1f °C (%.1f–%.1f)", cold_thr_C, cold_lo_C, cold_hi_C), `Cold slope, % per °C` = ci(cold_slope, cold_slope_lo, cold_slope_hi),
                `Heat: risk starts above` = fifelse(group %in% HS, sprintf("%.1f °C (%.1f–%.1f)", heat_thr_C, heat_lo_C, heat_hi_C), "no significant heat effect"),
                `Heat slope, % per °C` = sprintf("%s, p = %.3f", ci(heat_slope, heat_slope_lo, heat_slope_hi), p_heat_slope),
                `Safe band (°C)` = fifelse(group %in% HS, sprintf("%.1f", band_C), "—"))], "T6c_threshold_model")
RL <- tab("06_thresholds", "risk_level_temperatures"); cf <- function(e, l, h) sprintf("%.1f (%s–%s)", e, ifelse(is.finite(l), sprintf("%.1f", l), "never"), ifelse(is.finite(h), sprintf("%.1f", h), "never"))
md(ord(dcast(RL[, .(group, key = paste0(side, k), v = fifelse(side == "heat" & !group %in% HS, "—", cf(T_k, lo, hi)))], group ~ key, value.var = "v"))[
  , .(Group = GL[group], `Cold +5%` = cold5, `Cold +10%` = cold10, `Heat +5%` = heat5, `Heat +10%` = heat10)], "T6d_risk_level")
PT <- tab("06_thresholds", "risk_level_pair_tests"); TP <- tab("06_thresholds", "tm_pair_tests")
md(merge(PT, TP[, .(pair, cold_LR_p, heat_LR_p, cold_p_null, heat_p_null)], by = "pair")[order(match(pair, names(PLB)), side), .(Comparison = PLB[pair], Side = side,
   `Onset: p (LR; null)` = sprintf("%.2f; %.2f", fifelse(side == "cold", cold_LR_p, heat_LR_p), fifelse(side == "cold", cold_p_null, heat_p_null)),
   `Slope ratio at common threshold (p)` = sprintf("%.2f (p = %.3f)", slope_ratio, p_slope_common), `Slope test false positives (null)` = sprintf("%.0f%%", 100 * slope_test_false_positive_null),
   `Cold +10% temperature gap, °C (95% CI)` = fifelse(side == "cold", sprintf("%+.1f (%.1f to %s)", T10_diff, T10_diff_lo, fifelse(is.finite(T10_diff_hi), sprintf("%.1f", T10_diff_hi), "not reached")), "—"),
   `p (draws; null)` = fifelse(side == "cold", sprintf("%.2f; %.2f", p_T10_draws, p_T10_null), "—"))], "T6e_pair_tests")
DA <- tab("06_thresholds", "risk_level_days_at_risk")[(side == "cold" & k == 10) | (side == "heat" & k == 5)]
DA[, ns := side == "heat" & !group %in% HS]
md(DA[order(side, match(group, GROUPS)), .(Rule = fifelse(side == "cold", "Cold, +10%", "Heat, +5%"), Group = GL[group], `Own threshold (°C)` = fifelse(ns, "no significant heat effect", sprintf("%.1f", own_Tk)), `Population threshold (°C)` = sprintf("%.1f", common_Tk),
          `Person-days at risk: own / population rule` = fifelse(ns, sprintf("— / %.1f%%", share_days_common), sprintf("%.1f%% / %.1f%%", share_days_own, share_days_common)),
          `Missed days` = fifelse(ns, "—", sprintf("%.1f%%", missed_days_pct)), `Missed attributable deaths` = fifelse(ns, "—", sprintf("%.1f%%", missed_attr_pct)),
          `Over-flagged days` = fifelse(ns, "—", sprintf("%.1f%%", extra_days_pct)))], "T6f_days_at_risk")
NSm <- tab("06_thresholds", "null_split_summary")
md(dcast(NSm[, .(metric, prob = fifelse(prob == .5, "50/50 split", "20/80 split"), v = sprintf("%.0f", q95_abs_diff))], metric ~ prob, value.var = "v")[
  , .(Metric = metric, `Spurious difference, 95th pct (percentile points), 50/50` = `50/50 split`, `20/80` = `20/80 split`)], "T6g_null_split")
# F6: safe band, onset and +10% / +5% points by group
fb <- merge(TG[, .(group, cold_on = cold_thr_C, heat_on = heat_thr_C, p_heat_slope, p_cold_slope)],
            dcast(RL[(side == "cold" & k == 10) | (side == "heat" & k == 5), .(group, key = side, T_k)], group ~ key, value.var = "T_k"), by = "group")
fb[, `:=`(g = factor(GL[group], rev(GL)), heat_on = fifelse(p_heat_slope < .05, heat_on, NA_real_), heat = fifelse(p_heat_slope < .05 & is.finite(heat), heat, NA_real_))]
fg(ggplot(fb, aes(y = g)) + geom_segment(aes(x = cold_on, xend = fifelse(is.na(heat_on), 26, heat_on), yend = g), colour = COL[["lgrey"]], linewidth = 5) +
     geom_segment(aes(x = cold, xend = cold_on, yend = g), colour = COL[["blue"]], linewidth = 1) + geom_segment(aes(x = heat_on, xend = heat, yend = g), colour = COL[["orange"]], linewidth = 1, na.rm = TRUE) +
     geom_point(aes(x = cold_on), shape = 21, fill = "white", colour = COL[["blue"]], size = 2.6) + geom_point(aes(x = cold), colour = COL[["blue"]], size = 2.8) +
     geom_point(aes(x = heat_on), shape = 21, fill = "white", colour = COL[["orange"]], size = 2.6, na.rm = TRUE) + geom_point(aes(x = heat), colour = COL[["orange"]], size = 2.8, na.rm = TRUE) +
     labs(title = "Where temperature-related mortality starts, and how soon it becomes large",
          subtitle = "Grey: band without excess mortality. Open circle: estimated threshold (cold / heat).\nFilled: Tmin at which deaths are raised by 10% (cold) or 5% (heat). No heat marker = no significant heat effect.",
          x = "Daily minimum temperature (°C, population-weighted)", y = NULL), "F6_risk_levels", 10.5, 4.4)
# F7: missed attributable deaths under the population rule
fm <- DA[, .(group, ns, rule = fifelse(side == "cold", "Cold, +10% rule", "Heat, +5% rule"), missed_attr_pct = fifelse(ns, 0, missed_attr_pct))][, `:=`(g = factor(GL[group], rev(GL)), lab = fifelse(ns, "no significant heat effect", sprintf("%.1f%%", missed_attr_pct)))]
fg(ggplot(fm, aes(missed_attr_pct, g)) + geom_col(fill = COL[["orange"]], width = .6) + geom_text(aes(label = lab), hjust = -.15, size = 3, colour = COL[["ink2"]]) + facet_wrap(~rule) +
     scale_x_continuous(expand = expansion(mult = c(0, .2))) +
     labs(title = "What a single population threshold misses", subtitle = "Share of each group's temperature-attributable deaths occurring on days when its own risk already exceeds the level but the all-deaths rule does not flag the day.",
          x = "% of the group's attributable deaths missed", y = NULL), "F7_missed_deaths", 10, 4)

# ---- copies of the main figures used in the paper ---------------------------------------------------------------------------------
file.copy(file.path(OUT, "figures", "08_figures", c("Fig1_curves_by_margin.png", "Fig4_cpoisson_slopes.png", "Fig7a_dtr_sex_variants.png")),
          file.path(P_DIR, "figures"), overwrite = TRUE)
cat("paper assets written to", P_DIR, "\n")
