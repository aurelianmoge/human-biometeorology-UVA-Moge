# Annual county regressions are separate from daily exposure-response models.
source("2. Scripts/00_shared/setup.R")
theme <- THEMES[4]
county <- read.csv(file.path(DATA, "reference/county_baseline_analysis_data.csv")) |>
  mutate(Year = factor(Year), Poverty5 = poverty_pct / 5, Black10 = Black_share_pct / 10)

cluster_covariance <- function(fit, cluster) {
  design <- model.matrix(fit)
  residual <- residuals(fit)
  score <- rowsum(design * as.numeric(residual), cluster, reorder = FALSE)
  n <- nrow(design)
  p <- ncol(design)
  clusters <- nrow(score)
  bread <- solve(crossprod(design))
  correction <- clusters / (clusters - 1) * (n - 1) / (n - p)
  correction * bread %*% crossprod(score) %*% bread
}
results <- list()
for (outcome in c("Crude_rate", "Age_adjusted_rate")) {
  for (specification in c("poverty", "joint", "geography")) {
    right <- switch(specification, poverty = "Year + Poverty5", joint = "Year + Poverty5 + Black10",
                    geography = "Year + Poverty5 + Black10 + Longitude + Latitude + Independent_city")
    formula <- as.formula(paste0("log(", outcome, ") ~ ", right))
    fit <- lm(formula, data = county, na.action = na.fail)
    covariance <- cluster_covariance(fit, county$FIPS)
    se <- sqrt(diag(covariance))
    coefficients <- coef(fit)
    terms <- intersect(names(coefficients), c("Poverty5", "Black10"))
    critical <- qt(.975, length(unique(county$FIPS)) - 1)
    for (term in terms) results[[length(results) + 1]] <- data.frame(Outcome = outcome, Specification = specification,
      Term = term, Beta = coefficients[term], SE = se[term], Ratio = exp(coefficients[term]),
      Lower = exp(coefficients[term] - critical * se[term]), Upper = exp(coefficients[term] + critical * se[term]),
      Counties = length(unique(county$FIPS)), Observations = nrow(county),
      Interpretation = "annual ecological log-rate association; county-cluster covariance; not acute weather modification")
  }
}
results <- bind_rows(results)
write_table(results, "county_social_gradient_regressions", theme)
save_plot(ggplot(results, aes(Ratio, interaction(Outcome, Specification), color = Term)) +
  geom_vline(xintercept = 1, linetype = 2) + geom_point(position = position_dodge(.5)) +
  geom_errorbar(aes(xmin = Lower, xmax = Upper), width = .2, position = position_dodge(.5), orientation = "y") +
  labs(title = "Annual county mortality gradients", x = "Multiplicative rate index", y = NULL,
       caption = "Poverty5: +5 pp poverty; Black10: +10 pp Black share. Ecological covariate associations."),
  "county_regression_forest", theme, 12, 7)
