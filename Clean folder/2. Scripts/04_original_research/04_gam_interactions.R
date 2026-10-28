# Penalized GAMs describe same-day nonlinearity and temperature–DTR interaction.
# They answer a different lag question from cumulative DLNMs.
source("2. Scripts/00_shared/setup.R")
library(mgcv)
theme <- THEMES[4]
data <- read_derived("dtr_daily")
data$Time <- as.numeric(data$Date)
data$Weekday <- factor(as.integer(format(data$Date, "%u")))
complete <- complete.cases(data[c("All", "Tmean_midrange", "DTR", "RH07", "Time", "Weekday")])
training <- data[complete, ]

# Tensor interaction preserves separate smooth main effects. Penalized smooths
# are selected by REML, rather than an arbitrary polynomial order.
formulas <- list(
  additive = All ~ s(Tmean_midrange, k = 8) + s(DTR, k = 8) + s(RH07, k = 5) + s(Time, k = 112) + Weekday,
  interaction = All ~ s(Tmean_midrange, k = 8) + s(DTR, k = 8) + ti(Tmean_midrange, DTR, k = c(6, 6)) + s(RH07, k = 5) + s(Time, k = 112) + Weekday
)
statistics <- list()
surfaces <- list()
checks <- list()
for (name in names(formulas)) {
  message("Fitting penalized GAM: ", name)
  fit <- gam(formulas[[name]], data = training, family = quasipoisson(), method = "REML")
  saveRDS(fit, file.path(OUTPUT, theme, "models", paste0("GAM_", name, ".rds")))
  smooths <- as.data.frame(summary(fit)$s.table)
  smooths$Term <- rownames(smooths)
  smooths$Model <- name
  smooths$Inference_status <- "approximate nominal smooth tests; serial dependence not corrected; exploratory diagnostic only"
  statistics[[name]] <- smooths
  checks[[name]] <- data.frame(Model = name, Observations = nobs(fit), Converged = fit$converged,
                               Deviance_explained = summary(fit)$dev.expl, Scale = summary(fit)$scale,
                               Maximum_concurvity = max(concurvity(fit, full = TRUE), na.rm = TRUE))
  grid <- expand.grid(Tmean_midrange = seq(quantile(training$Tmean_midrange, .05), quantile(training$Tmean_midrange, .95), length.out = 55),
                      DTR = seq(quantile(training$DTR, .05), quantile(training$DTR, .95), length.out = 55))
  grid$RH07 <- median(training$RH07)
  grid$Time <- median(training$Time)
  grid$Weekday <- factor(levels(training$Weekday)[1], levels = levels(training$Weekday))
  reference <- grid[1, ]
  reference$Tmean_midrange <- median(training$Tmean_midrange)
  reference$DTR <- median(training$DTR)
  x <- predict(fit, grid, type = "lpmatrix")
  x0 <- predict(fit, reference, type = "lpmatrix")
  # Only weather terms change; same-day calendar, weekday and humidity are held fixed.
  difference <- sweep(x, 2, as.numeric(x0), `-`)
  grid$RR <- exp(drop(difference %*% coef(fit)))
  # Local support filter: retain cells near an observed temperature–DTR pair.
  scaled_observed <- cbind(training$Tmean_midrange / sd(training$Tmean_midrange), training$DTR / sd(training$DTR))
  distances <- vapply(seq_len(nrow(grid)), function(i) {
    target <- c(grid$Tmean_midrange[i] / sd(training$Tmean_midrange), grid$DTR[i] / sd(training$DTR))
    min(rowSums(sweep(scaled_observed, 2, target, `-`)^2))
  }, numeric(1))
  grid$Locally_supported <- distances <= .15^2
  grid$RR[!grid$Locally_supported] <- NA_real_
  grid$Model <- name
  surfaces[[name]] <- grid
}
write_table(bind_rows(statistics), "GAM_smooth_diagnostics", theme)
write_table(bind_rows(checks), "GAM_model_diagnostics", theme)
write_table(bind_rows(surfaces), "GAM_surface_data", theme)
save_plot(ggplot(bind_rows(surfaces), aes(Tmean_midrange, DTR, fill = RR)) + geom_raster() +
  facet_wrap(~Model) + scale_fill_viridis_c(na.value = "grey95") +
  guides(fill = guide_colorbar(barwidth = grid::unit(7, "cm"), barheight = grid::unit(.4, "cm"))) +
  labs(title = "Same-day GAM temperature–DTR response surface", x = "Midrange temperature (°C)", y = "DTR (°C)",
       caption = "Relative to joint medians; poorly supported cells masked. Same-day exploratory model, not cumulative or causal effects."),
  "GAM_temperature_DTR_surfaces", theme, 13, 6)
