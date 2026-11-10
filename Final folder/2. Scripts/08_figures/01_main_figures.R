# =============================================================================
# 01_main_figures.R - figures of the paper, built from saved tables and models only.
#   Fig1 pooled cumulative Tmin-mortality curves by margin
#   Fig2 group contrasts at extreme cold / heat across designs (two-stage, published design)
#   Fig3 group MMTs with bootstrap intervals (+ null-simulation spread when available)
#   Fig4 conditional Poisson cold / heat slopes by margin
#   Fig5 cold-attributable deaths per 100,000 by margin
#   Fig6 Cao comparison: person-days at risk, own vs all-cause curve
#   Fig7 DTR: Female / Male ratio across variants and lag-response by sex
#   Fig8 sensitivity grid: Female / Male RRR at P1 and P99
# Output: 3. Outputs/figures/08_figures/
# =============================================================================
source("2. Scripts/00_setup/config.R"); source("2. Scripts/00_setup/dlnm_tools.R")
STAGE <- "08_figures"
tab <- function(stage, name) fread(out_path("tables", stage, paste0(name, ".csv")))
PL <- read_model("pooled_Tmin")
FAM <- c(All = "All", Female = "Sex", Male = "Sex", Black = "Race", White = "Race", Under50 = "Age", Age50_74 = "Age", Age75plus = "Age")
PAIR_LABEL <- c("Female/Male" = "Female / Male", "Black/White" = "Black / White", "Age75plus/Under50" = "75+ / Under 50",
                "Age75plus/Age50_74" = "75+ / 50–74", "Age50_74/Under50" = "50–74 / Under 50")

# ---- Fig1 ------------------------------------------------------------------------------------------------
ref <- PL[["All|cum"]]$mmt
cur <- rbindlist(lapply(GROUPS, function(g) rr_curve(PL[[paste(g, "cum", sep = "|")]], ref)[pct >= .01 & pct <= .99][, group := g]))
allc <- cur[group == "All"][, k := 1][data.table(k = 1, fam = c("Sex", "Race", "Age")), on = "k", allow.cartesian = TRUE]
cg <- cur[group != "All"][, fam := FAM[group]]
f1 <- ggplot(cg, aes(100 * pct, exp(logrr), colour = group, fill = group)) + geom_hline(yintercept = 1, colour = COL[["grey"]]) +
  geom_line(data = allc, aes(100 * pct, exp(logrr)), inherit.aes = FALSE, colour = COL[["grey"]], linetype = 2) +
  geom_ribbon(aes(ymin = exp(logrr - 1.96 * se), ymax = exp(logrr + 1.96 * se)), alpha = .12, colour = NA) + geom_line(linewidth = .8) +
  geom_text(data = cg[pct == .01], aes(x = 0, label = GROUP_LABEL[group]), hjust = 1.1, size = 3.1, show.legend = FALSE) +
  facet_wrap(~factor(fam, c("Sex", "Race", "Age"))) + scale_colour_manual(values = GROUP_COL, guide = "none") +
  scale_fill_manual(values = GROUP_COL, guide = "none") + scale_x_continuous(limits = c(-22, 100), breaks = c(0, 25, 50, 75, 100)) +
  labs(title = "Cumulative (lag 0–21) temperature–mortality curves by margin", subtitle = "Two-stage DLNM, 12 regions, Tmin percentiles P1–P99; reference = all-cause MMT (P89). Dashed = all deaths.",
       x = "Regional percentile of daily minimum temperature", y = "Relative risk")
save_fig(f1, STAGE, "Fig1_curves_by_margin", 10.5, 4)

# ---- Fig2 ------------------------------------------------------------------------------------------------
pr <- tab("03_temperature", "pointwise_RRR")[exposure == "Tmin" & window == "cum" & centre == "MMT_all" & pct %in% c(1, 99)]
pd <- tab("05_group_differences", "published_design_contrasts")[exposure == "Tmin" & df_year == 3 & point %in% c("-10C", "24C")]
d2 <- rbind(pr[, .(pair, end = ifelse(pct == 1, "Extreme cold", "Extreme heat"), design = "Two-stage, 12 regions (P1 / P99 vs MMT)", RRR, lo, hi)],
            pd[, .(pair, end = ifelse(point == "-10C", "Extreme cold", "Extreme heat"), design = "Published statewide design (−10 / 24 °C vs 5 °C)", RRR, lo = RRR_lo, hi = RRR_hi)])
d2[, pl := factor(PAIR_LABEL[pair], rev(PAIR_LABEL))]
f2 <- ggplot(d2, aes(RRR, pl, colour = design)) + geom_vline(xintercept = 1, colour = COL[["grey"]]) +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(.5), size = .3) + facet_wrap(~end, scales = "free_x") +
  scale_colour_manual(values = c(COL[["blue"]], COL[["orange"]]), name = NULL) + scale_x_log10(breaks = c(.8, .9, 1, 1.1, 1.2, 1.3)) +
  labs(title = "Ratio of relative risks between margins, two designs", subtitle = "95% CI; groups modelled as independent series.", x = "RR ratio (log scale)", y = NULL)
save_fig(f2, STAGE, "Fig2_group_contrasts", 10, 4.2)

# ---- Fig3 ------------------------------------------------------------------------------------------------
mm <- rbindlist(lapply(setdiff(GROUPS, "All"), function(g) { p <- PL[[paste(g, "cum", sep = "|")]]
  data.table(group = g, mmt = 100 * p$mmt, lo = 100 * quantile(p$mmt_boot, .025), hi = 100 * quantile(p$mmt_boot, .975)) }))
mm[, g := factor(GROUP_LABEL[group], rev(GROUP_LABEL[setdiff(GROUPS, "All")]))]
f3 <- ggplot(mm, aes(mmt, g)) + geom_vline(xintercept = 100 * ref, linetype = 2, colour = COL[["grey"]])
nf <- out_path("tables", "07_methods_supplement", "null_groups.csv")
if (file.exists(nf)) { ng <- fread(nf)[, group := c(sB = "Black", sW = "White", sF = "Female", sM = "Male", sU50 = "Under50", s5074 = "Age50_74", s75 = "Age75plus")[group]]
  nq <- ng[, .(lo = 100 * quantile(mmt, .05), hi = 100 * quantile(mmt, .95)), by = group][, g := factor(GROUP_LABEL[group], levels(mm$g))]
  f3 <- f3 + geom_rect(data = nq, aes(xmin = lo, xmax = hi, ymin = as.numeric(g) - .3, ymax = as.numeric(g) + .3), inherit.aes = FALSE, fill = COL[["lgrey"]]) }
f3 <- f3 + geom_errorbarh(aes(xmin = lo, xmax = hi), height = .2, colour = COL[["orange"]]) + geom_point(colour = COL[["orange"]], size = 2.6) +
  scale_x_continuous(limits = c(0, 100), labels = function(x) paste0("P", x)) +
  labs(title = "Group-specific MMTs are poorly identified", subtitle = "Orange: estimate and 1,000-draw bootstrap 95% interval. Grey: 5–95% range when every group has the all-cause curve (null simulation).",
       x = "MMT (percentile of regional Tmin)", y = NULL)
save_fig(f3, STAGE, "Fig3_MMT_uncertainty", 9, 4)

# ---- Fig4 ------------------------------------------------------------------------------------------------
cs <- tab("05_group_differences", "cpoisson_slopes")[grepl("^M[234]", model) & grepl(":(cold|heat)$", term)]
cs[, `:=`(group = sub("grp::(\\w+):.*", "\\1", term), what = ifelse(grepl("cold$", term), "Cold: per °C of 21-day Tmin below P90", "Heat: per °C of 3-day Tmin above P90"))]
cs[, g := factor(GROUP_LABEL[group], rev(GROUP_LABEL[setdiff(GROUPS, "All")]))]
f4 <- ggplot(cs, aes(effect_pct, g, colour = group)) + geom_vline(xintercept = 0, colour = COL[["grey"]]) + geom_pointrange(aes(xmin = lo_pct, xmax = hi_pct), size = .3) +
  facet_wrap(~what, scales = "free_x") + scale_colour_manual(values = GROUP_COL, guide = "none") +
  labs(title = "Cold and heat slopes by margin, time-stratified conditional Poisson", subtitle = "Strata region × year × month × weekday; SEs clustered by region-month.", x = "% change in daily deaths (95% CI)", y = NULL)
save_fig(f4, STAGE, "Fig4_cpoisson_slopes", 9.5, 4)

# ---- Fig5 ------------------------------------------------------------------------------------------------
br <- tab("05_group_differences", "burden_rates")[centre == "common" & group != "All"]
br[, `:=`(g = factor(GROUP_LABEL[group], rev(GROUP_LABEL[setdiff(GROUPS, "All")])), fam = factor(FAM[group], c("Sex", "Race", "Age")))]
f5 <- ggplot(br, aes(cold_rate, g)) + geom_errorbarh(aes(xmin = pmax(cold_lo, .5), xmax = cold_hi), height = .2, colour = COL[["grey"]]) +
  geom_point(colour = COL[["blue"]], size = 2.6) + geom_text(aes(label = sprintf("%.0f", cold_rate)), vjust = -1, size = 3, colour = COL[["ink2"]]) +
  facet_grid(fam ~ ., scales = "free_y", space = "free_y") + scale_x_log10() +
  labs(title = "Cold-attributable deaths per 100,000 person-years, 2005–2020", subtitle = "Counterfactual = common all-cause MMT; 95% simulation interval.", x = "Rate (log scale)", y = NULL)
save_fig(f5, STAGE, "Fig5_cold_burden", 8.5, 4.4)

# ---- Fig6 ------------------------------------------------------------------------------------------------
dr <- tab("06_days_at_risk", "days_at_risk")[exposure == "Tmin" & rule %in% c("MMT_heat", "MMT_cold", "RR1.05_cold")]
cl <- function(x) pmin(pmax(1 + x / 100, .05), 20)
dr[, `:=`(ratio = cl(delta_pct), lo = cl(delta_lo), hi = cl(delta_hi), nlo = cl(noise_lo), nhi = cl(noise_hi), g = factor(GROUP_LABEL[group], rev(GROUP_LABEL[setdiff(GROUPS, "All")])))]
dr[, rl := factor(c(MMT_heat = "Heat: days warmer than the MMT", MMT_cold = "Cold: days colder than the MMT", RR1.05_cold = "Cold: days with RR ≥ 1.05")[rule],
                  c("Heat: days warmer than the MMT", "Cold: days colder than the MMT", "Cold: days with RR ≥ 1.05"))]
f6 <- ggplot(dr, aes(ratio, g)) + geom_vline(xintercept = 1, colour = COL[["grey"]]) +
  geom_rect(aes(xmin = nlo, xmax = nhi, ymin = as.numeric(g) - .32, ymax = as.numeric(g) + .32), fill = COL[["lgrey"]]) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = .15, colour = COL[["blue"]], alpha = .7) + geom_point(aes(colour = p_vs_noise < .05), size = 2.4) +
  geom_text(aes(label = sprintf("%+.0f%%", delta_pct)), vjust = -1, size = 2.8, colour = COL[["ink2"]]) + facet_wrap(~rl) +
  scale_colour_manual(values = c(`FALSE` = COL[["blue"]], `TRUE` = COL[["orange"]]), guide = "none") +
  scale_x_log10(breaks = c(.1, .5, 1, 2, 10), labels = function(x) sprintf("%+.0f%%", 100 * (x - 1))) +
  labs(title = "Person-days at risk: own curve vs all-cause curve (Cao et al. 2026 comparison)",
       subtitle = "Blue bar: 95% interval. Grey band: range produced by noise alone (identical true curves). Orange: outside it.", x = "Change in person-days at risk (log scale)", y = NULL)
save_fig(f6, STAGE, "Fig6_days_at_risk", 11, 4.2)

# ---- Fig7 ------------------------------------------------------------------------------------------------
dc <- tab("04_dtr", "dtr_contrasts")[pair == "Female/Male" & pct == 10 & window == "cum"]
vl <- c(base = "Base (DTR + Tmean cross-basis)", humidity = "+ humidity, dew point", anomaly = "Seasonally centred DTR", lag7 = "Lag 0–7 only", cold = "Outcomes Nov–Mar", warm = "Outcomes May–Sep")
dc[, v := factor(vl[variant], rev(vl))]
f7a <- ggplot(dc, aes(RRR, v)) + geom_vline(xintercept = 1, colour = COL[["grey"]]) + geom_pointrange(aes(xmin = lo, xmax = hi), colour = COL[["blue"]], size = .3) +
  geom_text(aes(x = max(hi) * 1.02, label = sprintf("p = %.2f", p_curve)), hjust = 0, size = 2.9, colour = COL[["ink2"]]) +
  scale_x_continuous(expand = expansion(mult = c(.05, .25))) +
  labs(title = "DTR: Female / Male ratio at the regional P10 of DTR (vs median), lag 0–21", subtitle = "Two-stage DLNM; p = whole-curve test.", x = "Female / Male RR ratio", y = NULL)
save_fig(f7a, STAGE, "Fig7a_dtr_sex_variants", 8.5, 3.6)
lr <- tab("04_dtr", "dtr_lag_response_P10_by_sex")[group != "All"]
f7b <- ggplot(lr, aes(lag, RR, colour = group, fill = group)) + geom_hline(yintercept = 1, colour = COL[["grey"]]) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .12, colour = NA) + geom_line(linewidth = .8) +
  scale_colour_manual(values = GROUP_COL, name = NULL) + scale_fill_manual(values = GROUP_COL, guide = "none") +
  labs(title = "Lag-response of mortality to a low-DTR day (regional P10 vs median), by sex", x = "Lag (days)", y = "Relative risk")
save_fig(f7b, STAGE, "Fig7b_dtr_lag_by_sex", 8.5, 3.6)

# ---- Fig8 ------------------------------------------------------------------------------------------------
sg <- tab("03_temperature", "sensitivity_contrasts")[pair == "Female/Male" & !degenerate]
sg[, `:=`(end = ifelse(pct == 1, "Cold end (P1)", "Heat end (P99)"), s = factor(spec, rev(unique(spec))))]
f8 <- ggplot(sg, aes(RRR, s, colour = end)) + geom_vline(xintercept = 1, colour = COL[["grey"]]) + geom_pointrange(aes(xmin = lo, xmax = hi), size = .25) +
  facet_wrap(~end, scales = "free_x") + scale_colour_manual(values = c(COL[["blue"]], COL[["orange"]]), guide = "none") +
  labs(title = "Female / Male RR ratio across model specifications", subtitle = "One factor changed at a time; reference = all-cause MMT of each specification.", x = "Female / Male RR ratio (95% CI)", y = NULL)
save_fig(f8, STAGE, "Fig8_sensitivity_sex", 9.5, 4.6)
cat("figures written to", out_path("figures", STAGE, ""), "\n")
