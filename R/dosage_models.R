# Standalone dosage-response model functions.
#
# These functions were redesigned from exploratory analysis code so that every
# input is explicit, every fit returns a structured result, and no private
# workspace objects or files are required. Only packages included with R are
# used.

.rbind_fill <- function(frames) {
  frames <- Filter(function(x) !is.null(x) && NROW(x) > 0L, frames)
  if (length(frames) == 0L) {
    return(data.frame())
  }

  all_names <- unique(unlist(lapply(frames, names), use.names = FALSE))
  frames <- lapply(frames, function(frame) {
    missing_names <- setdiff(all_names, names(frame))
    for (name in missing_names) {
      frame[[name]] <- NA
    }
    frame[all_names]
  })
  rownames_result <- NULL
  result <- do.call(rbind, frames)
  rownames(result) <- rownames_result
  result
}

.assert_scalar_string <- function(value, name) {
  if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
    stop(sprintf("`%s` must be one non-empty character string.", name), call. = FALSE)
  }
}

.validate_columns <- function(data, columns) {
  missing_columns <- setdiff(columns, names(data))
  if (length(missing_columns) > 0L) {
    stop(
      sprintf("Missing required column(s): %s", paste(missing_columns, collapse = ", ")),
      call. = FALSE
    )
  }
}

.model_data <- function(
    data,
    dose_col,
    response_col,
    group_col = NULL,
    min_rows = 5L,
    require_positive_dose = FALSE,
    na_action = c("fail", "omit")) {
  na_action <- match.arg(na_action)
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  .assert_scalar_string(dose_col, "dose_col")
  .assert_scalar_string(response_col, "response_col")
  if (!is.null(group_col)) {
    .assert_scalar_string(group_col, "group_col")
  }

  required <- c(dose_col, response_col, group_col)
  .validate_columns(data, required)
  if (!is.numeric(data[[dose_col]]) || !is.numeric(data[[response_col]])) {
    stop("Dose and response columns must both be numeric.", call. = FALSE)
  }

  # Preserve original row identities for diagnostics; omitting incomplete values
  # explicitly changes the fitted analysis population.
  complete <- is.finite(data[[dose_col]]) & is.finite(data[[response_col]])
  if (!is.null(group_col)) {
    complete <- complete & !is.na(data[[group_col]])
  }
  if (any(!complete) && identical(na_action, "fail")) {
    stop(
      sprintf(
        "Found %d row(s) with missing or non-finite model values; use `na_action = \"omit\"` explicitly to omit them.",
        sum(!complete)
      ),
      call. = FALSE
    )
  }

  filtered <- data[complete, , drop = FALSE]
  if (NROW(filtered) < min_rows) {
    stop(sprintf("At least %d complete observations are required.", min_rows), call. = FALSE)
  }
  if (length(unique(filtered[[dose_col]])) < 3L) {
    stop("At least three distinct dose values are required.", call. = FALSE)
  }
  if (require_positive_dose && any(filtered[[dose_col]] <= 0)) {
    stop("The log-logistic model requires every dose to be greater than zero.", call. = FALSE)
  }
  if (stats::sd(filtered[[response_col]]) == 0) {
    stop("The response has zero variance and cannot be modelled.", call. = FALSE)
  }

  result <- data.frame(
    .row = which(complete),
    .dose = filtered[[dose_col]],
    .response = filtered[[response_col]],
    stringsAsFactors = FALSE
  )
  if (!is.null(group_col)) {
    result$.group <- as.character(filtered[[group_col]])
  }
  result
}

# Use one Gaussian residual likelihood (including residual variance) for AIC,
# AICc, and BIC across candidates fitted to the same response observations.
.gaussian_fit_statistics <- function(observed, fitted, mean_parameters) {
  residuals <- observed - fitted
  n <- length(observed)
  rss <- sum(residuals^2)
  # Floor only the likelihood RSS so a perfect numerical fit remains defined.
  rss_for_likelihood <- max(
    rss,
    .Machine$double.eps * max(1, sum(observed^2))
  )
  sigma2_mle <- rss_for_likelihood / n
  log_likelihood <- -0.5 * n * (log(2 * pi) + 1 + log(sigma2_mle))
  total_parameters <- mean_parameters + 1 # Includes residual variance.
  aic <- -2 * log_likelihood + 2 * total_parameters
  aicc <- if (n > total_parameters + 1) {
    aic + 2 * total_parameters * (total_parameters + 1) /
      (n - total_parameters - 1)
  } else {
    NA_real_
  }

  list(
    n = n,
    rss = rss,
    rmse = sqrt(mean(residuals^2)),
    log_likelihood = log_likelihood,
    parameters = total_parameters,
    aic = aic,
    aicc = aicc,
    bic = -2 * log_likelihood + log(n) * total_parameters
  )
}

# Prediction points are equally spaced on the raw dose scale; a log x-axis is
# a plotting choice and does not change model fitting or this grid.
.prediction_grid <- function(dose, prediction_points, prediction_range = NULL) {
  if (!is.numeric(prediction_points) || length(prediction_points) != 1L ||
      prediction_points < 2L || prediction_points != as.integer(prediction_points)) {
    stop("`prediction_points` must be an integer of at least 2.", call. = FALSE)
  }
  if (is.null(prediction_range)) {
    prediction_range <- range(dose)
  }
  if (!is.numeric(prediction_range) || length(prediction_range) != 2L ||
      any(!is.finite(prediction_range)) || prediction_range[1] >= prediction_range[2]) {
    stop("`prediction_range` must contain two increasing finite values.", call. = FALSE)
  }
  seq(prediction_range[1], prediction_range[2], length.out = prediction_points)
}

#' Normalize dose and response columns to the control-group means.
#'
#' Normalization is optional and explicit. This reproduces the historical
#' analysis convention x / mean(control x) and y / mean(control y) without
#' hiding the operation inside a model-fitting function.
prepare_dosage_data <- function(
    data,
    dose_col,
    response_cols,
    group_col = NULL,
    control_label = NULL,
    normalize_dose = FALSE,
    normalize_response = FALSE) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  .assert_scalar_string(dose_col, "dose_col")
  if (!is.character(response_cols) || length(response_cols) < 1L || any(!nzchar(response_cols))) {
    stop("`response_cols` must contain at least one column name.", call. = FALSE)
  }
  response_cols <- unique(response_cols)
  .validate_columns(data, c(dose_col, response_cols, group_col))

  numeric_columns <- c(dose_col, response_cols)
  if (any(!vapply(data[numeric_columns], is.numeric, logical(1)))) {
    stop("Dose and response columns must be numeric.", call. = FALSE)
  }
  if (any(!vapply(data[numeric_columns], function(x) all(is.finite(x)), logical(1)))) {
    stop("Dose and response columns must not contain missing or non-finite values.", call. = FALSE)
  }

  normalize_any <- isTRUE(normalize_dose) || isTRUE(normalize_response)
  if (normalize_any) {
    if (is.null(group_col) || is.null(control_label)) {
      stop("Control normalization requires both `group_col` and `control_label`.", call. = FALSE)
    }
    .assert_scalar_string(group_col, "group_col")
    .assert_scalar_string(control_label, "control_label")
  }

  result <- data
  factors <- data.frame(
    column = character(),
    operation = character(),
    control_label = character(),
    divisor = numeric(),
    stringsAsFactors = FALSE
  )

  if (normalize_any) {
    control_rows <- !is.na(result[[group_col]]) & result[[group_col]] == control_label
    if (!any(control_rows)) {
      stop(sprintf("No rows match control label `%s`.", control_label), call. = FALSE)
    }

    columns_to_normalize <- c(
      if (isTRUE(normalize_dose)) dose_col,
      if (isTRUE(normalize_response)) response_cols
    )
    for (column in columns_to_normalize) {
      divisor <- mean(result[[column]][control_rows])
      if (!is.finite(divisor) || divisor == 0) {
        stop(
          sprintf("Control mean for `%s` is zero or non-finite; normalization is undefined.", column),
          call. = FALSE
        )
      }
      result[[column]] <- result[[column]] / divisor
      factors <- rbind(
        factors,
        data.frame(
          column = column,
          operation = "divide_by_control_mean",
          control_label = control_label,
          divisor = divisor,
          stringsAsFactors = FALSE
        )
      )
    }
  }

  structure(
    list(data = result, normalization = factors),
    class = "prepared_dosage_data"
  )
}

#' Fit linear and quadratic models and compare them using Gaussian AIC.
runmodelsaic_poly <- function(
    data,
    dose_col,
    response_col,
    degrees = 1:2,
    group_col = NULL,
    prediction_points = 201L,
    prediction_range = NULL,
    na_action = c("fail", "omit")) {
  na_action <- match.arg(na_action)
  if (!is.numeric(degrees) || length(degrees) < 1L || any(!is.finite(degrees)) ||
      any(degrees < 1L) || any(degrees != as.integer(degrees))) {
    stop("`degrees` must contain positive integers.", call. = FALSE)
  }
  degrees <- sort(unique(as.integer(degrees)))
  if (any(degrees > 2L)) {
    stop("Only degrees 1 (linear) and 2 (quadratic) are supported.", call. = FALSE)
  }
  model_data <- .model_data(
    data,
    dose_col,
    response_col,
    group_col,
    min_rows = max(degrees) + 2L,
    na_action = na_action
  )
  if (length(unique(model_data$.dose)) <= max(degrees)) {
    stop("Polynomial degree must be smaller than the number of distinct doses.", call. = FALSE)
  }
  prediction_dose <- .prediction_grid(
    model_data$.dose,
    prediction_points,
    prediction_range
  )

  models <- list()
  tables <- list()
  coefficients <- list()
  predictions <- list()

  # Raw powers retain interpretable dose and dose-squared coefficients rather
  # than the orthogonal-polynomial basis used by R's default `poly()`.
  for (degree in degrees) {
    model_name <- sprintf("polynomial_degree_%d", degree)
    fit <- stats::lm(
      .response ~ stats::poly(.dose, degree = degree, raw = TRUE),
      data = model_data
    )
    models[[model_name]] <- fit
    fitted_values <- stats::fitted(fit)
    fit_stats <- .gaussian_fit_statistics(
      model_data$.response,
      fitted_values,
      mean_parameters = degree + 1
    )
    fit_summary <- summary(fit)

    tables[[model_name]] <- data.frame(
      method = "polynomial",
      model = model_name,
      converged = TRUE,
      n = fit_stats$n,
      mean_parameters = degree + 1,
      total_parameters = fit_stats$parameters,
      degree = degree,
      effective_df = NA_real_,
      rss = fit_stats$rss,
      rmse = fit_stats$rmse,
      r_squared = fit_summary$r.squared,
      adjusted_r_squared = fit_summary$adj.r.squared,
      log_likelihood = fit_stats$log_likelihood,
      aic = fit_stats$aic,
      aicc = fit_stats$aicc,
      bic = fit_stats$bic,
      convergence_message = "",
      stringsAsFactors = FALSE
    )

    coefficient_matrix <- fit_summary$coefficients
    coefficients[[model_name]] <- data.frame(
      method = "polynomial",
      model = model_name,
      term = c("intercept", paste0("dose^", seq_len(degree))),
      estimate = coefficient_matrix[, "Estimate"],
      standard_error = coefficient_matrix[, "Std. Error"],
      statistic = coefficient_matrix[, "t value"],
      p_value = coefficient_matrix[, "Pr(>|t|)"],
      fixed = FALSE,
      stringsAsFactors = FALSE,
      row.names = NULL
    )

    predictions[[model_name]] <- data.frame(
      method = "polynomial",
      model = model_name,
      dose = prediction_dose,
      predicted = as.numeric(stats::predict(
        fit,
        newdata = data.frame(.dose = prediction_dose)
      )),
      stringsAsFactors = FALSE
    )
  }

  nested_tests <- list()
  if (length(degrees) > 1L) {
    for (index in 2:length(degrees)) {
      lower_name <- sprintf("polynomial_degree_%d", degrees[index - 1L])
      higher_name <- sprintf("polynomial_degree_%d", degrees[index])
      # This F-test is valid only for these nested polynomial models; AIC,
      # not this test, compares the polynomial family with LL.4 and splines.
      comparison <- stats::anova(models[[lower_name]], models[[higher_name]])
      nested_tests[[higher_name]] <- data.frame(
        smaller_model = lower_name,
        larger_model = higher_name,
        numerator_df = comparison$Df[2L],
        f_statistic = comparison$F[2L],
        p_value = comparison$`Pr(>F)`[2L],
        stringsAsFactors = FALSE
      )
    }
  }

  model_table <- do.call(rbind, tables)
  rownames(model_table) <- NULL
  best_model <- model_table$model[which.min(model_table$aic)]

  structure(
    list(
      method = "polynomial",
      dose_col = dose_col,
      response_col = response_col,
      group_col = group_col,
      input = model_data,
      models = models,
      model_table = model_table,
      coefficients = do.call(rbind, coefficients),
      nested_tests = .rbind_fill(nested_tests),
      predictions = do.call(rbind, predictions),
      best_model = best_model,
      settings = list(
        degrees = degrees,
        prediction_points = prediction_points,
        prediction_range = range(prediction_dose),
        na_action = na_action
      )
    ),
    class = "dosage_model_fit"
  )
}

#' Fit a four-parameter log-logistic model with drc::LL.4().
#'
#' The model is y = c + (d - c) / (1 + exp(b * (log(x) - log(e)))),
#' where b is the slope, c and d are the lower and upper asymptotes, and e is
#' the ED50. All four parameters are estimated with `drc::drm()`.
runmodelsaic_drc <- function(
    data,
    dose_col,
    response_col,
    group_col = NULL,
    prediction_points = 201L,
    prediction_range = NULL,
    na_action = c("fail", "omit")) {
  na_action <- match.arg(na_action)
  if (!requireNamespace("drc", quietly = TRUE)) {
    stop(
      "`runmodelsaic_drc()` requires the `drc` package. Install it before fitting LL.4 models.",
      call. = FALSE
    )
  }
  model_data <- .model_data(
    data,
    dose_col,
    response_col,
    group_col,
    min_rows = 6L,
    require_positive_dose = TRUE,
    na_action = na_action
  )
  prediction_dose <- .prediction_grid(
    model_data$.dose,
    prediction_points,
    prediction_range
  )
  if (any(prediction_dose <= 0)) {
    stop("The LL.4 prediction range must remain greater than zero.", call. = FALSE)
  }

  warning_messages <- character()
  captured_messages <- character()
  fit_capture <- tryCatch((function() {
    # drc can print raw optimizer diagnostics while raising an error. Capture
    # both streams, while retaining warning text for successful fits. Register
    # cleanup before setup so partial sink/device/temp-file setup cannot leak.
    sink_file <- NULL
    sink_connection <- NULL
    output_sink_depth <- sink.number(type = "output")
    message_sink_depth <- sink.number(type = "message")
    graphics_devices_before <- grDevices::dev.list()
    on.exit({
      graphics_devices_after <- grDevices::dev.list()
      opened_graphics_devices <- setdiff(
        as.integer(graphics_devices_after),
        as.integer(graphics_devices_before)
      )
      for (device in opened_graphics_devices) {
        try(grDevices::dev.off(which = device), silent = TRUE)
      }
      try({
        while (sink.number(type = "message") > message_sink_depth) {
          sink(type = "message")
        }
      }, silent = TRUE)
      try({
        while (sink.number(type = "output") > output_sink_depth) {
          sink(type = "output")
        }
      }, silent = TRUE)
      if (inherits(sink_connection, "connection")) {
        try(close(sink_connection), silent = TRUE)
      }
      if (!is.null(sink_file)) {
        try(unlink(sink_file), silent = TRUE)
      }
    }, add = TRUE)

    sink_file <- tempfile(fileext = ".ll4-output")
    sink_connection <- file(sink_file, open = "wt")
    sink(sink_connection, type = "output")
    sink(sink_connection, type = "message")
    fit <- tryCatch(
      withCallingHandlers(
        drc::drm(
          .response ~ .dose,
          data = model_data,
          fct = drc::LL.4(),
          control = drc::drmc(noMessage = TRUE)
        ),
        warning = function(warning) {
          warning_messages <<- c(warning_messages, conditionMessage(warning))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(error) error
    )
    # Unsink before reading so all output is flushed to disk. Cleanup remains
    # registered for errors in this finalization step.
    while (sink.number(type = "message") > message_sink_depth) sink(type = "message")
    while (sink.number(type = "output") > output_sink_depth) sink(type = "output")
    if (inherits(sink_connection, "connection")) close(sink_connection)
    sink_connection <- NULL
    captured <- if (file.exists(sink_file)) {
      readLines(sink_file, warn = FALSE)
    } else character()
    list(fit = fit, captured = captured)
  })(), error = function(error) list(fit = error, captured = character()))
  fit <- fit_capture$fit
  captured_messages <- fit_capture$captured

  model_name <- "ll4"
  if (inherits(fit, "error")) {
    model_table <- data.frame(
      method = "drc_ll4",
      model = model_name,
      converged = FALSE,
      n = NROW(model_data),
      mean_parameters = 4L,
      total_parameters = 5L,
      degree = NA_integer_,
      effective_df = NA_real_,
      rss = NA_real_,
      rmse = NA_real_,
      r_squared = NA_real_,
      adjusted_r_squared = NA_real_,
      log_likelihood = NA_real_,
      aic = NA_real_,
      aicc = NA_real_,
      bic = NA_real_,
      convergence_message = paste(unique(c(
        warning_messages,
        conditionMessage(fit),
        captured_messages[nzchar(trimws(captured_messages))]
      )), collapse = " | "),
      stringsAsFactors = FALSE
    )
    coefficients <- data.frame()
    predictions <- data.frame(
      method = "drc_ll4",
      model = model_name,
      dose = prediction_dose,
      predicted = NA_real_,
      stringsAsFactors = FALSE
    )
    models <- list(ll4 = NULL)
    best_model <- NA_character_
  } else {
    fitted_values <- as.numeric(stats::fitted(fit))
    fit_stats <- .gaussian_fit_statistics(
      model_data$.response,
      fitted_values,
      mean_parameters = 4L
    )
    # Use drc's likelihood and parameter count for its nonlinear LL.4 fit;
    # failed fits remain explicit NA rows rather than aborting other candidates.
    fit_log_likelihood <- stats::logLik(fit)
    fit_stats$log_likelihood <- as.numeric(fit_log_likelihood)
    fit_stats$parameters <- attr(fit_log_likelihood, "df")
    fit_stats$aic <- stats::AIC(fit)
    fit_stats$aicc <- if (fit_stats$n > fit_stats$parameters + 1L) {
      fit_stats$aic + 2 * fit_stats$parameters * (fit_stats$parameters + 1L) /
        (fit_stats$n - fit_stats$parameters - 1L)
    } else {
      NA_real_
    }
    fit_stats$bic <- stats::BIC(fit)
    tss <- sum((model_data$.response - mean(model_data$.response))^2)
    fit_summary <- summary(fit)$coefficients

    model_table <- data.frame(
      method = "drc_ll4",
      model = model_name,
      converged = TRUE,
      n = fit_stats$n,
      mean_parameters = 4L,
      total_parameters = fit_stats$parameters,
      degree = NA_integer_,
      effective_df = NA_real_,
      rss = fit_stats$rss,
      rmse = fit_stats$rmse,
      r_squared = 1 - fit_stats$rss / tss,
      adjusted_r_squared = NA_real_,
      log_likelihood = fit_stats$log_likelihood,
      aic = fit_stats$aic,
      aicc = fit_stats$aicc,
      bic = fit_stats$bic,
      convergence_message = paste(unique(c(
        warning_messages,
        captured_messages[nzchar(trimws(captured_messages))]
      )), collapse = " | "),
      stringsAsFactors = FALSE
    )
    coefficients <- data.frame(
      method = "drc_ll4",
      model = model_name,
      term = rownames(fit_summary),
      estimate = fit_summary[, "Estimate"],
      standard_error = fit_summary[, "Std. Error"],
      statistic = fit_summary[, "t-value"],
      p_value = fit_summary[, "p-value"],
      fixed = FALSE,
      stringsAsFactors = FALSE,
      row.names = NULL
    )
    predictions <- data.frame(
      method = "drc_ll4",
      model = model_name,
      dose = prediction_dose,
      predicted = as.numeric(stats::predict(
        fit,
        newdata = data.frame(.dose = prediction_dose)
      )),
      stringsAsFactors = FALSE
    )
    models <- list(ll4 = fit)
    best_model <- model_name
  }

  structure(
    list(
      method = "drc_ll4",
      dose_col = dose_col,
      response_col = response_col,
      group_col = group_col,
      input = model_data,
      models = models,
      model_table = model_table,
      coefficients = coefficients,
      nested_tests = data.frame(),
      predictions = predictions,
      best_model = best_model,
      settings = list(
        model_function = "drc::LL.4",
        prediction_points = prediction_points,
        prediction_range = range(prediction_dose),
        na_action = na_action,
        package_version = as.character(utils::packageVersion("drc"))
      )
    ),
    class = "dosage_model_fit"
  )
}

#' Fit fixed-effective-df smoothing splines and calculate approximate AIC.
runmodelsaic_spline <- function(
    data,
    dose_col,
    response_col,
    df = 3,
    group_col = NULL,
    prediction_points = 201L,
    prediction_range = NULL,
    na_action = c("fail", "omit")) {
  na_action <- match.arg(na_action)
  if (!is.numeric(df) || length(df) < 1L || any(!is.finite(df)) || any(df <= 1)) {
    stop("`df` must contain numeric values greater than 1.", call. = FALSE)
  }
  df <- sort(unique(df))
  model_data <- .model_data(
    data,
    dose_col,
    response_col,
    group_col,
    min_rows = max(6L, ceiling(max(df)) + 2L),
    na_action = na_action
  )
  prediction_dose <- .prediction_grid(
    model_data$.dose,
    prediction_points,
    prediction_range
  )

  models <- list()
  tables <- list()
  predictions <- list()
  warnings <- list()

  # `df` is a smoothing target; the achieved effective df is used as the
  # mean-parameter count for the approximate spline information criteria.
  for (requested_df in df) {
    model_name <- paste0("spline_df_", format(requested_df, trim = TRUE))
    warning_message <- ""
    fit <- withCallingHandlers(
      stats::smooth.spline(
        x = model_data$.dose,
        y = model_data$.response,
        df = requested_df,
        control.spar = list(high = 2.5)
      ),
      warning = function(warning) {
        warning_message <<- conditionMessage(warning)
        invokeRestart("muffleWarning")
      }
    )
    models[[model_name]] <- fit
    fitted_values <- stats::predict(fit, x = model_data$.dose)$y
    fit_stats <- .gaussian_fit_statistics(
      model_data$.response,
      fitted_values,
      mean_parameters = fit$df
    )
    tss <- sum((model_data$.response - mean(model_data$.response))^2)

    tables[[model_name]] <- data.frame(
      method = "smoothing_spline",
      model = model_name,
      converged = TRUE,
      n = fit_stats$n,
      mean_parameters = fit$df,
      total_parameters = fit_stats$parameters,
      degree = NA_integer_,
      requested_df = requested_df,
      effective_df = fit$df,
      rss = fit_stats$rss,
      rmse = fit_stats$rmse,
      r_squared = 1 - fit_stats$rss / tss,
      adjusted_r_squared = NA_real_,
      log_likelihood = fit_stats$log_likelihood,
      aic = fit_stats$aic,
      aicc = fit_stats$aicc,
      bic = fit_stats$bic,
      convergence_message = warning_message,
      stringsAsFactors = FALSE
    )
    predictions[[model_name]] <- data.frame(
      method = "smoothing_spline",
      model = model_name,
      dose = prediction_dose,
      predicted = stats::predict(fit, x = prediction_dose)$y,
      stringsAsFactors = FALSE
    )
    warnings[[model_name]] <- warning_message
  }

  model_table <- .rbind_fill(tables)
  best_model <- model_table$model[which.min(model_table$aic)]

  structure(
    list(
      method = "smoothing_spline",
      dose_col = dose_col,
      response_col = response_col,
      group_col = group_col,
      input = model_data,
      models = models,
      model_table = model_table,
      coefficients = data.frame(),
      nested_tests = data.frame(),
      predictions = do.call(rbind, predictions),
      best_model = best_model,
      warnings = warnings,
      settings = list(
        requested_df = df,
        prediction_points = prediction_points,
        prediction_range = range(prediction_dose),
        na_action = na_action,
        aic_caution = paste(
          "Spline AIC uses effective degrees of freedom plus one residual-variance parameter.",
          "It is an approximation and should be interpreted with diagnostics, not alone."
        )
      )
    ),
    class = "dosage_model_fit"
  )
}

#' Combine fitted-model tables and calculate delta AIC and Akaike weights.
compare_dosage_models <- function(...) {
  fits <- list(...)
  if (length(fits) == 1L && is.list(fits[[1L]]) &&
      !inherits(fits[[1L]], "dosage_model_fit")) {
    fits <- fits[[1L]]
  }
  if (length(fits) < 1L || any(!vapply(fits, inherits, logical(1), "dosage_model_fit"))) {
    stop("Supply one or more `dosage_model_fit` objects.", call. = FALSE)
  }

  tables <- lapply(fits, `[[`, "model_table")
  result <- .rbind_fill(tables)
  # Weights are conditional on converged finite-AIC fits for this one response,
  # and are not probabilities that a biological mechanism is true.
  result <- result[result$converged & is.finite(result$aic), , drop = FALSE]
  if (NROW(result) == 0L) {
    stop("No converged model with a finite AIC is available.", call. = FALSE)
  }
  result <- result[order(result$aic, result$model), , drop = FALSE]
  result$delta_aic <- result$aic - min(result$aic)
  relative_likelihood <- exp(-0.5 * result$delta_aic)
  result$akaike_weight <- relative_likelihood / sum(relative_likelihood)
  result$aic_rank <- seq_len(NROW(result))
  rownames(result) <- NULL
  result
}

#' Plot observations with the AIC-best linear/quadratic and LL.4 curves.
plot_dosage_models <- function(
    data,
    dose_col,
    response_col,
    fits,
    group_col = NULL,
    file = NULL,
    width = 7,
    height = 5,
    log_dose = FALSE,
    main = NULL) {
  model_data <- .model_data(
    data,
    dose_col,
    response_col,
    group_col,
    min_rows = 5L,
    require_positive_dose = isTRUE(log_dose),
    na_action = "fail"
  )
  if (inherits(fits, "dosage_model_fit")) {
    fits <- list(fits)
  }
  if (!is.list(fits) || length(fits) < 1L ||
      any(!vapply(fits, inherits, logical(1), "dosage_model_fit"))) {
    stop("`fits` must contain one or more `dosage_model_fit` objects.", call. = FALSE)
  }

  opened_device <- FALSE
  if (!is.null(file)) {
    .assert_scalar_string(file, "file")
    directory <- dirname(file)
    if (!dir.exists(directory)) {
      dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    }
    extension <- tolower(tools::file_ext(file))
    if (identical(extension, "pdf")) {
      grDevices::pdf(file, width = width, height = height, useDingbats = FALSE)
    } else if (identical(extension, "png")) {
      grDevices::png(file, width = width, height = height, units = "in", res = 300)
    } else {
      stop("Plot file must end in `.pdf` or `.png`.", call. = FALSE)
    }
    opened_device <- TRUE
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  # Show only each family's AIC-best polynomial/LL.4 curve; splines remain in
  # tabular comparison by design because this graphic is descriptive.
  display_fits <- Filter(function(fit) !identical(fit$method, "smoothing_spline"), fits)
  if (length(display_fits) == 0L) {
    stop("At least one polynomial or drc::LL.4 fit is required for plotting.", call. = FALSE)
  }
  curve_frames <- lapply(display_fits, function(fit) {
    if (is.na(fit$best_model)) {
      return(NULL)
    }
    curve <- fit$predictions[fit$predictions$model == fit$best_model, , drop = FALSE]
    curve$family <- fit$method
    curve
  })
  curves <- .rbind_fill(curve_frames)
  y_range <- range(c(model_data$.response, curves$predicted), finite = TRUE)
  x_range <- range(model_data$.dose)

  if (is.null(main)) {
    main <- sprintf("%s versus %s", response_col, dose_col)
  }

  if (is.null(group_col)) {
    point_colors <- rep("grey30", NROW(model_data))
    groups <- NULL
    group_colors <- NULL
  } else {
    groups <- unique(model_data$.group)
    group_colors <- setNames(grDevices::hcl.colors(length(groups), "Dark 3"), groups)
    point_colors <- unname(group_colors[model_data$.group])
  }

  graphics::plot(
    model_data$.dose,
    model_data$.response,
    pch = 19,
    col = point_colors,
    xlab = dose_col,
    ylab = response_col,
    main = main,
    xlim = x_range,
    ylim = y_range,
    log = if (isTRUE(log_dose)) "x" else ""
  )

  family_colors <- c(
    polynomial = "#009E73",
    drc_ll4 = "#D55E00"
  )
  family_labels <- c(
    polynomial = "Polynomial",
    drc_ll4 = "drc::LL.4"
  )
  plotted_families <- character()
  for (fit in display_fits) {
    if (is.na(fit$best_model)) {
      next
    }
    curve <- fit$predictions[fit$predictions$model == fit$best_model, , drop = FALSE]
    valid <- is.finite(curve$dose) & is.finite(curve$predicted)
    if (!any(valid)) {
      next
    }
    graphics::lines(
      curve$dose[valid],
      curve$predicted[valid],
      col = family_colors[[fit$method]],
      lwd = 2
    )
    plotted_families <- c(plotted_families, fit$method)
  }

  plotted_families <- unique(plotted_families)
  if (length(plotted_families) > 0L) {
    if (!is.null(groups)) {
      graphics::legend(
        "topleft",
        legend = groups,
        col = group_colors,
        pch = 19,
        bty = "n",
        cex = 0.8
      )
    }
    graphics::legend(
      "bottomright",
      legend = unname(family_labels[plotted_families]),
      col = unname(family_colors[plotted_families]),
      lwd = 2,
      bty = "n",
      cex = 0.8
    )
  } else {
    graphics::text(
      mean(x_range),
      mean(y_range),
      labels = "No converged curve available",
      col = "grey30",
      cex = 1.1
    )
  }

  invisible(list(
    file = file,
    device_opened = opened_device,
    plotted_methods = unique(plotted_families)
  ))
}

#' Return observation-level diagnostics for every candidate in a dosage fit list.
#'
#' The `.row` field retained by the fitting functions is the row number in the
#' original data frame.  Failed candidates are represented by one row per input
#' observation, with unavailable quantities set to NA.  Influence quantities
#' are intentionally calculated only for stats::lm fits.
diagnose_dosage_models <- function(fits) {
  if (inherits(fits, "dosage_model_fit")) fits <- list(fits)
  if (!is.list(fits) || length(fits) < 1L ||
      any(!vapply(fits, inherits, logical(1), "dosage_model_fit"))) {
    stop("`fits` must contain one or more `dosage_model_fit` objects.", call. = FALSE)
  }

  rows <- list()
  for (family_fit in fits) {
    input <- family_fit$input
    if (!is.data.frame(input) || !all(c(".row", ".dose", ".response") %in% names(input))) {
      stop("Each dosage fit must retain `.row`, `.dose`, and `.response` input fields.", call. = FALSE)
    }
    candidate_names <- names(family_fit$models)
    if (length(candidate_names) == 0L) {
      candidate_names <- as.character(family_fit$model_table$model)
    }
    for (model_name in candidate_names) {
      table_row <- family_fit$model_table[
        as.character(family_fit$model_table$model) == model_name, , drop = FALSE
      ]
      if (NROW(table_row) == 0L) {
        table_row <- data.frame(
          converged = FALSE,
          convergence_message = "No model summary was available.",
          stringsAsFactors = FALSE
        )
      } else {
        table_row <- table_row[1L, , drop = FALSE]
      }
      converged <- isTRUE(table_row$converged[[1L]])
      message <- if ("convergence_message" %in% names(table_row)) {
        as.character(table_row$convergence_message[[1L]])
      } else ""
      model_fit <- family_fit$models[[model_name]]
      fitted_values <- rep(NA_real_, NROW(input))
      if (converged && !is.null(model_fit)) {
        fitted_values <- tryCatch({
          if (inherits(model_fit, "lm")) {
            as.numeric(stats::fitted(model_fit))
          } else if (inherits(model_fit, "smooth.spline")) {
            as.numeric(stats::predict(model_fit, x = input$.dose)$y)
          } else {
            as.numeric(stats::fitted(model_fit))
          }
        }, error = function(error) rep(NA_real_, NROW(input)))
        if (length(fitted_values) != NROW(input)) {
          fitted_values <- rep(NA_real_, NROW(input))
        }
      }
      residual_values <- input$.response - fitted_values
      standardized <- leverage <- cooks <- rep(NA_real_, NROW(input))
      if (inherits(model_fit, "lm") && converged) {
        standardized <- tryCatch(as.numeric(stats::rstandard(model_fit)),
          error = function(error) standardized)
        leverage <- tryCatch(as.numeric(stats::hatvalues(model_fit)),
          error = function(error) leverage)
        cooks <- tryCatch(as.numeric(stats::cooks.distance(model_fit)),
          error = function(error) cooks)
      }
      rows[[length(rows) + 1L]] <- data.frame(
        response = family_fit$response_col,
        model = model_name,
        family = family_fit$method,
        row_id = as.integer(input$.row),
        original_row_id = as.integer(input$.row),
        dose = as.numeric(input$.dose),
        observed = as.numeric(input$.response),
        fitted = fitted_values,
        residual = as.numeric(residual_values),
        standardized_residual = standardized,
        leverage = leverage,
        cooks_distance = cooks,
        converged = converged,
        convergence_status = if (converged) "converged" else "failed",
        convergence_message = message,
        stringsAsFactors = FALSE
      )
    }
  }
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

#' Draw observed-versus-fitted and residual-versus-fitted dosage diagnostics.
#'
#' All panels are descriptive checks.  No observation is removed or flagged for
#' automatic deletion; failed candidates receive an explicit annotation.
plot_dosage_model_diagnostics <- function(
    fits,
    file = NULL,
    width = 11,
    height = 7,
    main = NULL) {
  if (inherits(fits, "dosage_model_fit")) fits <- list(fits)
  diagnostics <- diagnose_dosage_models(fits)
  models <- unique(as.character(diagnostics$model))
  if (length(models) < 1L) stop("No candidate models are available.", call. = FALSE)

  device_opened <- FALSE
  if (!is.null(file)) {
    .assert_scalar_string(file, "file")
    directory <- dirname(file)
    if (!dir.exists(directory)) dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    extension <- tolower(tools::file_ext(file))
    if (extension == "pdf") {
      grDevices::pdf(file, width = width, height = height, useDingbats = FALSE)
    } else if (extension == "png") {
      grDevices::png(file, width = width, height = height, units = "in", res = 300)
    } else {
      stop("Diagnostic plot file must end in `.pdf` or `.png`.", call. = FALSE)
    }
    device_opened <- TRUE
  }
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    if (device_opened && grDevices::dev.cur() > 1L) grDevices::dev.off()
  }, add = TRUE)
  n_columns <- max(1L, min(3L, length(models)))
  graphics::par(mfrow = c(2L, n_columns), mar = c(4, 4, 2.8, 1), oma = c(0, 0, 2.2, 0))
  colors <- grDevices::hcl.colors(length(models), "Dark 3")
  names(colors) <- models
  for (model_name in models) {
    rows <- diagnostics[diagnostics$model == model_name, , drop = FALSE]
    valid <- is.finite(rows$observed) & is.finite(rows$fitted)
    panel_title <- sprintf("%s: observed vs fitted", model_name)
    if (any(valid)) {
      limits <- range(c(rows$observed[valid], rows$fitted[valid]), finite = TRUE)
      if (diff(limits) == 0) limits <- limits + c(-0.5, 0.5)
      graphics::plot(rows$fitted[valid], rows$observed[valid], pch = 19,
        col = colors[[model_name]], xlab = "Fitted", ylab = "Observed",
        main = panel_title, xlim = limits, ylim = limits)
      graphics::abline(0, 1, lty = 2, col = "grey45")
    } else {
      graphics::plot.new(); graphics::title(main = panel_title)
      graphics::text(0.5, 0.55, "Fit unavailable", col = "grey35")
      graphics::text(0.5, 0.42, unique(rows$convergence_message), col = "grey35", cex = 0.7)
    }
    valid_residual <- is.finite(rows$fitted) & is.finite(rows$residual)
    panel_title <- sprintf("%s: residual vs fitted", model_name)
    if (any(valid_residual)) {
      graphics::plot(rows$fitted[valid_residual], rows$residual[valid_residual], pch = 19,
        col = colors[[model_name]], xlab = "Fitted", ylab = "Residual", main = panel_title)
      graphics::abline(h = 0, lty = 2, col = "grey45")
    } else {
      graphics::plot.new(); graphics::title(main = panel_title)
      graphics::text(0.5, 0.5, "Fit unavailable", col = "grey35")
    }
  }
  graphics::mtext(if (is.null(main)) "Dosage model diagnostics" else main,
    side = 3, outer = TRUE, line = 0.6, cex = 1.1)
  graphics::mtext("Descriptive diagnostics; observations are not automatically deleted or treated as outliers.",
    side = 3, outer = TRUE, line = -0.6, cex = 0.8, col = "grey35")
  invisible(list(file = file, device_opened = device_opened, panels = 2L * length(models)))
}

# Short aliases for callers using the noun-first API.
dosage_model_diagnostics <- diagnose_dosage_models
plot_dosage_diagnostics <- plot_dosage_model_diagnostics

print.dosage_model_fit <- function(x, ...) {
  cat(sprintf("<dosage_model_fit: %s>\n", x$method))
  cat(sprintf("Response: %s; dose: %s\n", x$response_col, x$dose_col))
  cat(sprintf("Best model by AIC: %s\n", x$best_model))
  print(x$model_table, row.names = FALSE)
  invisible(x)
}
