# Distributed lag models written explicitly with base-R spline matrices.
# The exposure-major, lag-minor ordering is checked against direct convolution.
# dlnm is not installed locally; no approximation of its installed code is claimed.
source("2. Scripts/00_shared/setup.R")

shift_days <- function(values, lag) {
  if (lag == 0) return(values)
  c(rep(NA_real_, lag), head(values, -lag))
}

make_crossbasis <- function(exposure, maximum_lag = 21L) {
  boundaries <- range(exposure, na.rm = TRUE)
  knots <- seq(boundaries[1], boundaries[2], length.out = 4)[2:3]
  exposure_basis <- bs(exposure, degree = 2, knots = knots, Boundary.knots = boundaries)
  lag_basis <- ns(0:maximum_lag,
                  knots = exp(seq(log(1), log(maximum_lag), length.out = 4)[2:3]),
                  Boundary.knots = c(0, maximum_lag), intercept = TRUE)
  columns <- list()
  for (exposure_column in seq_len(ncol(exposure_basis))) {
    for (lag_column in seq_len(ncol(lag_basis))) {
      contributions <- lapply(0:maximum_lag, function(lag) {
        shift_days(exposure_basis[, exposure_column], lag) * lag_basis[lag + 1, lag_column]
      })
      columns[[length(columns) + 1]] <- Reduce(`+`, contributions)
    }
  }
  design <- do.call(cbind, columns)
  for (row in head(which(complete.cases(design)), 3)) {
    direct <- unlist(lapply(seq_len(ncol(exposure_basis)), function(column) {
      colSums(lag_basis * exposure_basis[row - (0:maximum_lag), column])
    }))
    stopifnot(max(abs(direct - design[row, ])) < 1e-9)
  }
  list(design = design, lag_basis = lag_basis, boundaries = boundaries, knots = knots,
       reference = median(exposure, na.rm = TRUE), maximum_lag = maximum_lag,
       quantiles = quantile(exposure, c(.01, .1, .5, .9, .95, .99), na.rm = TRUE, names = FALSE))
}

prediction_contrast <- function(basis, values, lags = 0:21) {
  evaluated <- bs(values, degree = 2, knots = basis$knots, Boundary.knots = basis$boundaries)
  reference <- bs(basis$reference, degree = 2, knots = basis$knots, Boundary.knots = basis$boundaries)
  difference <- sweep(evaluated, 2, as.numeric(reference), `-`)
  lag_sums <- colSums(basis$lag_basis[lags + 1, , drop = FALSE])
  do.call(cbind, lapply(seq_len(ncol(difference)), function(column) {
    difference[, column, drop = FALSE] %*% t(lag_sums)
  }))
}

calendar_hac <- function(influence, bandwidth = 28L) {
  # Ineligible dates remain zero rows, so lag 7 always means seven calendar days.
  influence <- as.matrix(influence)
  covariance <- crossprod(influence)
  for (lag in seq_len(bandwidth)) {
    lag_product <- crossprod(influence[(lag + 1):nrow(influence), , drop = FALSE],
                             influence[1:(nrow(influence) - lag), , drop = FALSE])
    covariance <- covariance + (1 - lag / (bandwidth + 1)) * (lag_product + t(lag_product))
  }
  (covariance + t(covariance)) / 2
}

make_design <- function(data, exposure, outcome = "All", time_df = 112L,
                        controls = "full", exclude_2020 = FALSE, interaction = FALSE) {
  basis <- make_crossbasis(data[[exposure]])
  interest <- basis$design
  core_histories <- lapply(c("Tmin", "Tmax", "Tmean_midrange"), function(variable) make_crossbasis(data[[variable]])$design)
  eligible <- complete.cases(do.call(cbind, core_histories)) &
    complete.cases(data[c("All", "Female", "Male", "Black", "White", "Under50", "Age75plus", "Respiratory", "Circulatory")])
  calendar <- as.numeric(data$Date)
  weekday <- factor(as.integer(format(data$Date, "%u")))
  nuisance <- cbind(Intercept = 1, ns(calendar, df = time_df), model.matrix(~weekday)[, -1, drop = FALSE])
  if (exposure == "DTR") {
    eligible <- eligible & is.finite(data$RH07) & complete.cases(basis$design)
    if (controls == "same_day") nuisance <- cbind(nuisance, ns(data$Tmean_midrange, df = 3))
    if (controls %in% c("full", "temperature")) nuisance <- cbind(nuisance, make_crossbasis(data$Tmean_midrange)$design)
    if (controls %in% c("full", "same_day")) nuisance <- cbind(nuisance, ns(data$RH07, df = 3))
  }
  if (interaction) {
    # Centering preserves the main DTR interpretation at the median midrange.
    product <- (data$DTR - basis$reference) *
      (data$Tmean_midrange - median(data$Tmean_midrange, na.rm = TRUE)) / 10
    interest <- cbind(interest, product)
  }
  eligible <- eligible & complete.cases(nuisance) & complete.cases(interest) & is.finite(data[[outcome]])
  if (exclude_2020) eligible <- eligible & format(data$Date, "%Y") != "2020"
  if (sum(eligible) < 500) stop("Fewer than 500 eligible dates.")
  nuisance_qr <- qr(nuisance[eligible, , drop = FALSE], tol = 1e-9)
  retained <- sort(nuisance_qr$pivot[seq_len(nuisance_qr$rank)])
  design <- cbind(nuisance[, retained, drop = FALSE], interest)
  if (qr(design[eligible, , drop = FALSE], tol = 1e-9)$rank != ncol(design)) stop("Rank deficient exposure design.")
  list(data = data, design = design, eligible = eligible, basis = basis, outcome = outcome,
       target = (ncol(design) - ncol(interest) + 1):ncol(design), interaction = interaction)
}

fit_count_model <- function(design, negative_binomial = FALSE) {
  eligible <- design$eligible
  x <- design$design[eligible, , drop = FALSE]
  y <- design$data[[design$outcome]][eligible]
  if (any(y < 0) || any(abs(y - round(y)) > 1e-8)) stop("Outcomes must be nonnegative counts.")
  theta <- NA_real_
  if (negative_binomial) {
    optimizer_warnings <- character()
    fit <- withCallingHandlers(
      MASS::glm.nb(y ~ x - 1, control = glm.control(maxit = 100)),
      warning = function(warning) {
        optimizer_warnings <<- c(optimizer_warnings, conditionMessage(warning))
        invokeRestart("muffleWarning")
      }
    )
    if (any(grepl("limit|failed|NaN", optimizer_warnings, ignore.case = TRUE))) {
      stop("Negative-binomial optimizer unresolved: ", paste(unique(optimizer_warnings), collapse = "; "))
    }
    theta <- fit$theta
  } else {
    fit <- glm.fit(x, y, family = quasipoisson("log"), control = glm.control(maxit = 100))
  }
  if (!fit$converged || any(!is.finite(fit$coefficients))) stop("Nonconverged or aliased count model.")
  mean <- as.numeric(fit$fitted.values)
  weights <- if (is.finite(theta)) mean / (1 + mean / theta) else mean
  target <- design$target
  residualized <- qr.resid(qr(x[, -target, drop = FALSE] * sqrt(weights)),
                          x[, target, drop = FALSE] * sqrt(weights))
  bread <- solve(crossprod(residualized))
  score <- if (is.finite(theta)) (y - mean) / (1 + mean / theta) else y - mean
  influence <- matrix(0, nrow(design$data), length(target))
  influence[eligible, ] <- (residualized * (score / sqrt(weights))) %*% bread
  residual <- rep(NA_real_, nrow(design$data))
  residual[eligible] <- (y - mean) / sqrt(mean)
  list(beta = as.numeric(fit$coefficients[target]), covariance = calendar_hac(influence),
       influence = influence, basis = design$basis, eligible = eligible,
       residual = residual, fitted = mean, dates = design$data$Date,
       coefficients = fit$coefficients, dispersion = sum((y - mean)^2 / mean) / (length(y) - ncol(x)),
       observations = length(y), parameters = ncol(x), deaths = sum(y),
       converged = fit$converged, theta = theta, interaction = design$interaction)
}

model_contrast <- function(fit, values, lags) {
  contrast <- prediction_contrast(fit$basis, values, lags)
  if (fit$interaction) contrast <- cbind(contrast, 0)
  contrast
}

summarize_contrast <- function(fit, contrast, bandwidth = 28L) {
  estimate <- drop(contrast %*% fit$beta)
  covariance <- if (bandwidth == 28) fit$covariance else calendar_hac(fit$influence, bandwidth)
  standard_error <- sqrt(pmax(0, rowSums((contrast %*% covariance) * contrast)))
  data.frame(Log_RR = estimate, SE = standard_error, RR = exp(estimate),
             Lower = exp(estimate - 1.96 * standard_error), Upper = exp(estimate + 1.96 * standard_error),
             Nominal_P = ifelse(standard_error > 0, 2 * pnorm(-abs(estimate / standard_error)), NA_real_))
}

paired_contrast <- function(first, second, first_contrast, second_contrast) {
  stopifnot(identical(first$eligible, second$eligible), identical(first$dates, second$dates))
  estimate <- drop(first_contrast %*% first$beta - second_contrast %*% second$beta)
  influence <- first$influence %*% t(first_contrast) - second$influence %*% t(second_contrast)
  se <- sqrt(pmax(0, diag(calendar_hac(influence))))
  data.frame(Log_RR = estimate, SE = se, RR = exp(estimate), Lower = exp(estimate - 1.96 * se),
             Upper = exp(estimate + 1.96 * se), Nominal_P = 2 * pnorm(-abs(estimate / se)))
}
