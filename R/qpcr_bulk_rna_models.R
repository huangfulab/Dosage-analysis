# Standalone qPCR--bulk RNA-seq concordance functions.
#
# The analysis unit is a gRNA species. Four sample-level observations are
# averaged for each species before fitting the three prespecified models:
# linear, quadratic, and drc::LL.4() with bulk RNA-seq asymptotes fixed to the
# observed gRNA-mean minimum and maximum.

.qpcr_bulk_assert_string = function(value, name) {
  if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
    stop(sprintf("`%s` must be one non-empty character string.", name), call. = F)
  }
}

.qpcr_bulk_validate_columns = function(data, columns) {
  missing_columns = setdiff(columns, names(data))
  if (length(missing_columns) > 0L) {
    stop(
      sprintf("Missing required column(s): %s", paste(missing_columns, collapse = ", ")),
      call. = F
    )
  }
}

# AICc is undefined when the gRNA-level sample cannot support the fitted
# parameters plus its small-sample correction.
.qpcr_bulk_aicc = function(aic, n, parameters) {
  if (n <= parameters + 1L) {
    return(NA_real_)
  }
  aic + 2 * parameters * (parameters + 1L) / (n - parameters - 1L)
}

# Derive in-sample likelihood metrics from each fitted curve; RMSE and R² are
# descriptive summaries of the limited gRNA-level data, not validation scores.
.qpcr_bulk_fit_metrics = function(fit, model, method) {
  fitted_values = as.numeric(stats::fitted(fit))
  residual_values = as.numeric(stats::residuals(fit))
  observed_values = fitted_values + residual_values
  log_likelihood = stats::logLik(fit)
  parameters = attr(log_likelihood, "df")
  rss = sum(residual_values^2)
  tss = sum((observed_values - mean(observed_values))^2)
  aic = stats::AIC(fit)

  data.frame(
    model = model,
    method = method,
    converged = T,
    n = length(observed_values),
    parameters = parameters,
    rss = rss,
    rmse = sqrt(mean(residual_values^2)),
    r_squared = 1 - rss / tss,
    log_likelihood = as.numeric(log_likelihood),
    aic = aic,
    aicc = .qpcr_bulk_aicc(aic, length(observed_values), parameters),
    bic = stats::BIC(fit),
    convergence_message = "",
    stringsAsFactors = F
  )
}

.qpcr_bulk_coefficients = function(fit, model, method) {
  coefficient_table = summary(fit)$coefficients
  statistic_column = if ("t value" %in% colnames(coefficient_table)) {
    "t value"
  } else {
    "t-value"
  }
  p_value_column = if ("Pr(>|t|)" %in% colnames(coefficient_table)) {
    "Pr(>|t|)"
  } else {
    "p-value"
  }

  data.frame(
    model = model,
    method = method,
    term = rownames(coefficient_table),
    estimate = coefficient_table[, "Estimate"],
    standard_error = coefficient_table[, "Std. Error"],
    statistic = coefficient_table[, statistic_column],
    p_value = coefficient_table[, p_value_column],
    fixed = FALSE,
    stringsAsFactors = F,
    row.names = NULL
  )
}

.qpcr_bulk_anova_frame = function(anova_result, analysis) {
  result = as.data.frame(anova_result, check.names = F)
  result$term = rownames(result)
  result$analysis = analysis
  result = result[c("analysis", "term", setdiff(names(result), c("analysis", "term")))]
  rownames(result) = NULL
  result
}

.qpcr_bulk_empty_diagnostic_rows = function(point_data, model, message = "") {
  data.frame(
    model = model,
    grna = point_data$grna,
    qpcr_percent = point_data$qpcr_percent,
    observed_bulk_rna_percent = point_data$bulk_rna_percent,
    fitted_bulk_rna_percent = NA_real_,
    residual = NA_real_,
    standardized_residual = NA_real_,
    leverage = NA_real_,
    cooks_distance = NA_real_,
    converged = FALSE,
    convergence_message = message,
    stringsAsFactors = FALSE
  )
}

#' Derive one observed/fitted diagnostic row per gRNA and candidate model.
#'
#' Standardized residual, leverage, and Cook's distance are supplied only for
#' stats::lm models.  They are not defined for LL.4 here and remain NA; a
#' failed candidate model is represented by rows of NAs rather than an error.
qpcr_bulk_rna_model_diagnostics = function(fit) {
  if (!inherits(fit, "qpcr_bulk_rna_model_fit")) {
    stop("`fit` must be returned by `fit_qpcr_bulk_rna_models()`.", call. = FALSE)
  }
  point_data = fit$point_data
  result = lapply(c("linear", "quadratic", "ll4"), function(model) {
    model_fit = fit$models[[model]]
    if (is.null(model_fit)) {
      return(.qpcr_bulk_empty_diagnostic_rows(
        point_data, model, "Model fit unavailable."
      ))
    }
    fitted_values = tryCatch(
      as.numeric(stats::fitted(model_fit)),
      error = function(error) rep(NA_real_, NROW(point_data))
    )
    if (length(fitted_values) != NROW(point_data)) {
      fitted_values = rep(NA_real_, NROW(point_data))
    }
    residual_values = point_data$bulk_rna_percent - fitted_values
    standardized = leverage = cooks = rep(NA_real_, NROW(point_data))
    if (inherits(model_fit, "lm")) {
      standardized = tryCatch(as.numeric(stats::rstandard(model_fit)), error = function(error) standardized)
      leverage = tryCatch(as.numeric(stats::hatvalues(model_fit)), error = function(error) leverage)
      cooks = tryCatch(as.numeric(stats::cooks.distance(model_fit)), error = function(error) cooks)
    }
    data.frame(
      model = model,
      grna = point_data$grna,
      qpcr_percent = point_data$qpcr_percent,
      observed_bulk_rna_percent = point_data$bulk_rna_percent,
      fitted_bulk_rna_percent = fitted_values,
      residual = residual_values,
      standardized_residual = standardized,
      leverage = leverage,
      cooks_distance = cooks,
      converged = TRUE,
      convergence_message = "",
      stringsAsFactors = FALSE
    )
  })
  result = do.call(rbind, result)
  rownames(result) = NULL
  result
}

# Backward/forward-friendly alias for callers that use the verb diagnose.
diagnose_qpcr_bulk_rna_models = qpcr_bulk_rna_model_diagnostics

#' Aggregate sample-level qPCR and bulk RNA-seq measurements by gRNA species.
#'
#' The expected four samples per gRNA and eleven resulting gRNA points are
#' validated explicitly because they define the planned analysis unit. The
#' point-level output also retains descriptive within-gRNA mean, SD, SE, minimum,
#' and maximum values for each assay.
aggregate_grna_expression = function(
    data,
    grna_col,
    qpcr_col,
    bulk_rna_col,
    expected_replicates = 4L,
    expected_points = 11L,
    na_action = c("fail", "omit")) {
  na_action = match.arg(na_action)
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = F)
  }
  .qpcr_bulk_assert_string(grna_col, "grna_col")
  .qpcr_bulk_assert_string(qpcr_col, "qpcr_col")
  .qpcr_bulk_assert_string(bulk_rna_col, "bulk_rna_col")
  .qpcr_bulk_validate_columns(data, c(grna_col, qpcr_col, bulk_rna_col))

  if (!is.numeric(data[[qpcr_col]]) || !is.numeric(data[[bulk_rna_col]])) {
    stop("qPCR and bulk RNA-seq columns must be numeric.", call. = F)
  }
  if (!is.null(expected_replicates) &&
      (!is.numeric(expected_replicates) || length(expected_replicates) != 1L ||
        expected_replicates < 1L || expected_replicates != as.integer(expected_replicates))) {
    stop("`expected_replicates` must be NULL or one positive integer.", call. = F)
  }
  if (!is.null(expected_points) &&
      (!is.numeric(expected_points) || length(expected_points) != 1L ||
        expected_points < 1L || expected_points != as.integer(expected_points))) {
    stop("`expected_points` must be NULL or one positive integer.", call. = F)
  }

  complete = !is.na(data[[grna_col]]) & nzchar(as.character(data[[grna_col]])) &
    is.finite(data[[qpcr_col]]) & is.finite(data[[bulk_rna_col]])
  if (any(!complete) && identical(na_action, "fail")) {
    stop(
      sprintf(
        "Found %d row(s) with missing or non-finite gRNA, qPCR, or bulk RNA-seq values.",
        sum(!complete)
      ),
      call. = F
    )
  }

  # Models use unweighted gRNA means; replicate SD/SE are retained only as
  # descriptive summaries and are not used as model weights or uncertainty.
  sample_data = data[complete, c(grna_col, qpcr_col, bulk_rna_col), drop = F]
  names(sample_data) = c("grna", "qpcr_percent", "bulk_rna_percent")
  sample_data$grna = as.character(sample_data$grna)
  if (NROW(sample_data) == 0L) {
    stop("No complete sample-level observations remain.", call. = F)
  }
  if (any(sample_data$qpcr_percent <= 0)) {
    stop(
      "Every qPCR percentage must be greater than zero because LL.4 is defined on log(x).",
      call. = F
    )
  }

  grouped_data = split(sample_data, sample_data$grna, drop = T)
  point_data = do.call(rbind, lapply(grouped_data, function(group) {
    n_samples = NROW(group)
    qpcr_mean = mean(group$qpcr_percent)
    bulk_rna_mean = mean(group$bulk_rna_percent)
    # SD/SE are deliberately descriptive within-gRNA summaries.  A single
    # observation has no estimable within-group variability.
    qpcr_sd = if (n_samples >= 2L) stats::sd(group$qpcr_percent) else NA_real_
    bulk_rna_sd = if (n_samples >= 2L) stats::sd(group$bulk_rna_percent) else NA_real_
    data.frame(
      grna = group$grna[1L],
      n_samples = n_samples,
      qpcr_percent = qpcr_mean,
      bulk_rna_percent = bulk_rna_mean,
      qpcr_mean = qpcr_mean,
      qpcr_sd = qpcr_sd,
      qpcr_se = if (n_samples >= 2L) qpcr_sd / sqrt(n_samples) else NA_real_,
      qpcr_min = min(group$qpcr_percent),
      qpcr_max = max(group$qpcr_percent),
      bulk_rna_mean = bulk_rna_mean,
      bulk_rna_sd = bulk_rna_sd,
      bulk_rna_se = if (n_samples >= 2L) bulk_rna_sd / sqrt(n_samples) else NA_real_,
      bulk_rna_min = min(group$bulk_rna_percent),
      bulk_rna_max = max(group$bulk_rna_percent),
      stringsAsFactors = F
    )
  }))
  rownames(point_data) = NULL
  point_data = point_data[order(point_data$grna), , drop = F]

  if (!is.null(expected_replicates) && any(point_data$n_samples != expected_replicates)) {
    incorrect_counts = point_data[point_data$n_samples != expected_replicates, , drop = F]
    stop(
      sprintf(
        "Each gRNA must have %d complete samples. Mismatches: %s",
        expected_replicates,
        paste(
          sprintf("%s=%d", incorrect_counts$grna, incorrect_counts$n_samples),
          collapse = ", "
        )
      ),
      call. = F
    )
  }
  if (!is.null(expected_points) && NROW(point_data) != expected_points) {
    stop(
      sprintf("Expected %d gRNA-level points but found %d.", expected_points, NROW(point_data)),
      call. = F
    )
  }
  if (NROW(point_data) < 6L || length(unique(point_data$qpcr_percent)) < 4L) {
    stop("At least six gRNA points with four distinct qPCR values are required.", call. = F)
  }
  if (stats::sd(point_data$bulk_rna_percent) == 0) {
    stop("The gRNA-level bulk RNA-seq response has zero variance.", call. = F)
  }

  structure(
    list(
      sample_data = sample_data,
      point_data = point_data,
      settings = list(
        grna_col = grna_col,
        qpcr_col = qpcr_col,
        bulk_rna_col = bulk_rna_col,
        expected_replicates = expected_replicates,
        expected_points = expected_points,
        na_action = na_action
      )
    ),
    class = "grna_expression_aggregate"
  )
}

#' Fit prespecified linear, quadratic, and fixed-asymptote LL.4 models.
#'
#' The LL.4 lower and upper asymptotes are fixed to the observed bulk RNA-seq
#' minimum and maximum on the fitted gRNA-level means. Only the slope and ED50
#' remain free.
fit_qpcr_bulk_rna_models = function(aggregate, prediction_points = 201L) {
  if (!inherits(aggregate, "grna_expression_aggregate")) {
    stop("`aggregate` must be returned by `aggregate_grna_expression()`.", call. = F)
  }
  if (!requireNamespace("drc", quietly = T)) {
    stop("`fit_qpcr_bulk_rna_models()` requires the `drc` package.", call. = F)
  }
  if (!is.numeric(prediction_points) || length(prediction_points) != 1L ||
      prediction_points < 2L || prediction_points != as.integer(prediction_points)) {
    stop("`prediction_points` must be an integer of at least 2.", call. = F)
  }

  # Fit bulk RNA-seq response against qPCR predictor on the gRNA-level means,
  # matching the prespecified analysis unit rather than the 44 raw samples.
  point_data = aggregate$point_data
  response_min = min(point_data$bulk_rna_percent)
  response_max = max(point_data$bulk_rna_percent)
  if (!is.finite(response_min) || !is.finite(response_max)) {
    stop("LL.4 asymptotes require finite bulk RNA-seq minimum and maximum values.", call. = F)
  }
  if (!(response_min < response_max)) {
    stop(
      "LL.4 asymptotes require a strictly smaller bulk RNA-seq minimum than maximum so c and d can be fixed.",
      call. = F
    )
  }
  prediction_data = data.frame(
    qpcr_percent = seq(
      min(point_data$qpcr_percent),
      max(point_data$qpcr_percent),
      length.out = prediction_points
    )
  )
  # Candidate curves are linear, raw-scale quadratic, and a positive-x LL.4
  # with response asymptotes fixed to the observed bulk RNA-seq range.
  linear = stats::lm(bulk_rna_percent ~ qpcr_percent, data = point_data)
  quadratic = stats::lm(
    bulk_rna_percent ~ qpcr_percent + I(qpcr_percent^2),
    data = point_data
  )

  ll4_warnings = character()
  ll4 = tryCatch((function() {
    # drc can print raw optimizer diagnostics while raising an error. Sink
    # both output streams, while retaining warning text for successful fits.
    # Register cleanup before either sink is configured: tempfile(), file(),
    # and sink() can each fail, and a partially configured capture must not
    # leak a sink, connection, or temporary file.
    sink_file = NULL
    sink_connection = NULL
    output_sink_depth = sink.number(type = "output")
    message_sink_connection = getConnection(sink.number(type="message"))
    graphics_devices_before = grDevices::dev.list()
    on.exit({
      # drc may open a default graphics device while reporting an optimizer
      # failure.  Close only devices absent from the caller's device set; in
      # particular, never close a device that was already open on entry.
      graphics_devices_after = grDevices::dev.list()
      opened_graphics_devices = setdiff(
        as.integer(graphics_devices_after),
        as.integer(graphics_devices_before)
      )
      for (device in opened_graphics_devices) {
        try(grDevices::dev.off(which = device), silent = TRUE)
      }
      try({
        sink(message_sink_connection, type="message")
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
    sink_file = tempfile(fileext = ".ll4-output")
    sink_connection = file(sink_file, open = "wt")
    sink(sink_connection, type = "output")
    sink(sink_connection, type = "message")
    withCallingHandlers(
      drc::drm(
        bulk_rna_percent ~ qpcr_percent,
        data = point_data,
        # Fix c/d to the fitted bulk RNA-seq range; estimate only b and e.
        fct = drc::LL.4(fixed = c(NA, response_min, response_max, NA)),
        control = drc::drmc(noMessage = T)
      ),
      warning = function(warning) {
        ll4_warnings <<- c(ll4_warnings, conditionMessage(warning))
        invokeRestart("muffleWarning")
      }
    )
  })(), error = function(error) error)

  model_table = rbind(
    .qpcr_bulk_fit_metrics(linear, "linear", "stats::lm"),
    .qpcr_bulk_fit_metrics(quadratic, "quadratic", "stats::lm")
  )
  coefficient_table = rbind(
    .qpcr_bulk_coefficients(linear, "linear", "stats::lm"),
    .qpcr_bulk_coefficients(quadratic, "quadratic", "stats::lm")
  )
  predictions = rbind(
    data.frame(
      model = "linear",
      qpcr_percent = prediction_data$qpcr_percent,
      predicted_bulk_rna_percent = as.numeric(stats::predict(linear, prediction_data)),
      stringsAsFactors = F
    ),
    data.frame(
      model = "quadratic",
      qpcr_percent = prediction_data$qpcr_percent,
      predicted_bulk_rna_percent = as.numeric(stats::predict(quadratic, prediction_data)),
      stringsAsFactors = F
    )
  )

  # Free mean parameters are slope b and ED50 e; c and d are fixed asymptotes.
  free_mean_parameters = 2L
  expected_total_parameters = free_mean_parameters + 1L
  fixed_ll4_coefficients = data.frame(
    model = "ll4",
    method = "drc::LL.4",
    term = c("c:(Intercept)", "d:(Intercept)"),
    estimate = c(response_min, response_max),
    standard_error = NA_real_,
    statistic = NA_real_,
    p_value = NA_real_,
    fixed = TRUE,
    stringsAsFactors = F
  )
  if (inherits(ll4, "error")) {
    ll4_metrics = data.frame(
      model = "ll4",
      method = "drc::LL.4",
      converged = F,
      n = NROW(point_data),
      parameters = expected_total_parameters,
      rss = NA_real_,
      rmse = NA_real_,
      r_squared = NA_real_,
      log_likelihood = NA_real_,
      aic = NA_real_,
      aicc = NA_real_,
      bic = NA_real_,
      convergence_message = conditionMessage(ll4),
      stringsAsFactors = F
    )
    ll4_coefficients = fixed_ll4_coefficients
    ll4_predictions = data.frame(
      model = "ll4",
      qpcr_percent = prediction_data$qpcr_percent,
      predicted_bulk_rna_percent = NA_real_,
      stringsAsFactors = F
    )
    ll4_model = NULL
  } else {
    ll4_metrics = .qpcr_bulk_fit_metrics(ll4, "ll4", "drc::LL.4")
    ll4_metrics$convergence_message = paste(unique(ll4_warnings), collapse = " | ")
    ll4_coefficients = rbind(
      .qpcr_bulk_coefficients(ll4, "ll4", "drc::LL.4"),
      fixed_ll4_coefficients
    )
    ll4_predictions = data.frame(
      model = "ll4",
      qpcr_percent = prediction_data$qpcr_percent,
      predicted_bulk_rna_percent = as.numeric(stats::predict(ll4, prediction_data)),
      stringsAsFactors = F
    )
    ll4_model = ll4
  }

  model_table = rbind(model_table, ll4_metrics)
  coefficient_table = rbind(coefficient_table, ll4_coefficients)
  predictions = rbind(predictions, ll4_predictions)
  # Rank only converged candidates fitted to these same gRNA means. Akaike
  # weights express relative in-sample support, not assay concordance or mechanism.
  valid_models = model_table$converged & is.finite(model_table$aic)
  model_table$delta_aic = NA_real_
  model_table$akaike_weight = NA_real_
  model_table$aic_rank = NA_integer_
  if (any(valid_models)) {
    model_table$delta_aic[valid_models] =
      model_table$aic[valid_models] - min(model_table$aic[valid_models])
    relative_likelihood = exp(-0.5 * model_table$delta_aic[valid_models])
    model_table$akaike_weight[valid_models] =
      relative_likelihood / sum(relative_likelihood)
    model_table$aic_rank[valid_models] =
      rank(model_table$aic[valid_models], ties.method = "first")
  }

  fit_result = structure(
    list(
      point_data = point_data,
      models = list(linear = linear, quadratic = quadratic, ll4 = ll4_model),
      model_table = model_table,
      coefficients = coefficient_table,
      predictions = predictions,
      # One-model ANOVAs are sequential lm term tests; the nested F-test below
      # compares only linear with quadratic, never LL.4 or assay concordance.
      linear_anova = .qpcr_bulk_anova_frame(stats::anova(linear), "linear"),
      quadratic_anova = .qpcr_bulk_anova_frame(stats::anova(quadratic), "quadratic"),
      linear_vs_quadratic_anova = .qpcr_bulk_anova_frame(
        stats::anova(linear, quadratic),
        "linear_vs_quadratic"
      ),
      settings = list(
        model_set = c("linear", "quadratic", "drc::LL.4"),
        ll4_fixed_parameters = c(c = response_min, d = response_max),
        ll4_free_parameters = c("b", "e"),
        prediction_points = prediction_points,
        drc_version = as.character(utils::packageVersion("drc"))
      )
    ),
    class = "qpcr_bulk_rna_model_fit"
  )
  fit_result$diagnostics = qpcr_bulk_rna_model_diagnostics(fit_result)
  fit_result
}

#' Draw the qPCR--bulk RNA-seq scatterplot on the original percentage scale.
plot_qpcr_bulk_rna_models = function(
    fit,
    file = NULL,
    width = 7,
    height = 5.5,
    label_points = T,
    main = "qPCR and bulk RNA-seq expression concordance") {
  if (!inherits(fit, "qpcr_bulk_rna_model_fit")) {
    stop("`fit` must be returned by `fit_qpcr_bulk_rna_models()`.", call. = F)
  }
  if (!is.logical(label_points) || length(label_points) != 1L || is.na(label_points)) {
    stop("`label_points` must be TRUE or FALSE.", call. = F)
  }

  device_opened = F
  if (!is.null(file)) {
    .qpcr_bulk_assert_string(file, "file")
    output_directory = dirname(file)
    if (!dir.exists(output_directory)) {
      dir.create(output_directory, recursive = T, showWarnings = F)
    }
    extension = tolower(tools::file_ext(file))
    if (identical(extension, "pdf")) {
      grDevices::pdf(file, width = width, height = height, useDingbats = F)
    } else if (identical(extension, "png")) {
      grDevices::png(file, width = width, height = height, units = "in", res = 300)
    } else {
      stop("Plot file must end in `.pdf` or `.png`.", call. = F)
    }
    device_opened = T
    on.exit(grDevices::dev.off(), add = T)
  }

  point_data = fit$point_data
  prediction_values = fit$predictions$predicted_bulk_rna_percent
  qpcr_error_low = point_data$qpcr_percent - point_data$qpcr_se
  qpcr_error_high = point_data$qpcr_percent + point_data$qpcr_se
  bulk_error_low = point_data$bulk_rna_percent - point_data$bulk_rna_se
  bulk_error_high = point_data$bulk_rna_percent + point_data$bulk_rna_se
  # Include the control-normalized 100% reference in both axes so departures
  # from control remain visible even if observations lie wholly on one side.
  x_limits = range(c(point_data$qpcr_percent, qpcr_error_low, qpcr_error_high, 100), finite = T)
  y_limits = range(c(point_data$bulk_rna_percent, bulk_error_low, bulk_error_high, prediction_values, 100), finite = T)
  x_span = diff(x_limits)
  y_span = diff(y_limits)
  x_limits = x_limits + c(-0.03, 0.10) * if (x_span > 0) x_span else 1
  y_limits = y_limits + c(-0.03, 0.06) * if (y_span > 0) y_span else 1
  graphics::plot(
    point_data$qpcr_percent,
    point_data$bulk_rna_percent,
    pch = 19,
    col = "#252525",
    xlab = "Normalized qPCR expression (% of controls)",
    ylab = "Normalized bulk RNA-seq expression (% of controls)",
    main = main,
    xlim = x_limits,
    ylim = y_limits
  )
  graphics::abline(v = 100, h = 100, lty = 2, col = "grey55")
  # These bars summarize replicate variation at each gRNA; they are
  # descriptive SEs and do not assert inferential independence of gRNAs.
  valid_x_se = is.finite(point_data$qpcr_se)
  valid_y_se = is.finite(point_data$bulk_rna_se)
  if (any(valid_x_se)) {
    graphics::arrows(
      qpcr_error_low[valid_x_se], point_data$bulk_rna_percent[valid_x_se],
      qpcr_error_high[valid_x_se], point_data$bulk_rna_percent[valid_x_se],
      angle = 90, code = 3, length = 0.04, col = "grey35", lwd = 1
    )
  }
  if (any(valid_y_se)) {
    graphics::arrows(
      point_data$qpcr_percent[valid_y_se], bulk_error_low[valid_y_se],
      point_data$qpcr_percent[valid_y_se], bulk_error_high[valid_y_se],
      angle = 90, code = 3, length = 0.04, col = "grey35", lwd = 1
    )
  }
  graphics::mtext("SE bars are descriptive replicate summaries; not inferential independence claims.",
    side = 1, line = 3.1, cex = 0.65, col = "grey35")
  if (isTRUE(label_points)) {
    graphics::text(
      point_data$qpcr_percent,
      point_data$bulk_rna_percent,
      labels = point_data$grna,
      pos = 3,
      cex = 0.75
    )
  }

  model_colors = c(linear = "#0072B2", quadratic = "#009E73", ll4 = "#D55E00")
  model_labels = c(
    linear = "Linear",
    quadratic = "Quadratic",
    ll4 = "LL.4"
  )
  plotted_models = character()
  for (model in names(model_colors)) {
    curve = fit$predictions[fit$predictions$model == model, , drop = F]
    valid = is.finite(curve$qpcr_percent) & is.finite(curve$predicted_bulk_rna_percent)
    if (any(valid)) {
      graphics::lines(
        curve$qpcr_percent[valid],
        curve$predicted_bulk_rna_percent[valid],
        col = model_colors[[model]],
        lwd = 2
      )
      plotted_models = c(plotted_models, model)
    }
  }
  if (length(plotted_models) > 0L) {
    graphics::legend(
      "bottomright",
      legend = unname(model_labels[plotted_models]),
      col = unname(model_colors[plotted_models]),
      lwd = 2,
      bty = "n"
    )
  }

  invisible(list(file = file, device_opened = device_opened))
}

#' Write observed-versus-fitted and residual-versus-fitted diagnostic panels.
#'
#' Panels use the complete candidate set, including a safe message for a failed
#' LL.4 fit. These plots are descriptive diagnostics, not validation or
#' inferential evidence.
plot_qpcr_bulk_rna_diagnostics = function(
    fit,
    file = NULL,
    width = 10,
    height = 6.5,
    main = "qPCR and bulk RNA-seq model diagnostics") {
  if (!inherits(fit, "qpcr_bulk_rna_model_fit")) {
    stop("`fit` must be returned by `fit_qpcr_bulk_rna_models()`.", call. = FALSE)
  }
  device_opened = FALSE
  # Capture caller parameters only when a caller device exists; otherwise
  # par() would create an unintended default Rplots.pdf.
  caller_has_device = !is.null(grDevices::dev.list())
  old_par = if (caller_has_device) graphics::par(no.readonly = TRUE) else NULL
  if (!is.null(file)) {
    .qpcr_bulk_assert_string(file, "file")
    output_directory = dirname(file)
    if (!dir.exists(output_directory)) dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
    extension = tolower(tools::file_ext(file))
    if (extension == "pdf") {
      grDevices::pdf(file, width = width, height = height, useDingbats = FALSE)
    } else if (extension == "png") {
      grDevices::png(file, width = width, height = height, units = "in", res = 300)
    } else {
      stop("Diagnostic plot file must end in `.pdf` or `.png`.", call. = FALSE)
    }
    device_opened = TRUE
  }
  # Restore the plotting device while it is active, then close only an
  # output device owned by this function. This also preserves file=NULL callers.
  on.exit({
    if (!is.null(old_par)) {
      graphics::par(old_par)
      # Setting a layout can reset cex; restore it after layout-dependent values.
      graphics::par(cex = old_par$cex)
    }
    if (device_opened) grDevices::dev.off()
  }, add = TRUE)
  diagnostics = fit$diagnostics
  graphics::par(mfrow = c(2, 3), mar = c(4, 4, 2.5, 1), oma = c(0, 0, 2, 0))
  colors = c(linear = "#0072B2", quadratic = "#009E73", ll4 = "#D55E00")
  for (model in c("linear", "quadratic", "ll4")) {
    rows = diagnostics[diagnostics$model == model, , drop = FALSE]
    valid = is.finite(rows$observed_bulk_rna_percent) & is.finite(rows$fitted_bulk_rna_percent)
    if (any(valid)) {
      limits = range(c(rows$observed_bulk_rna_percent[valid], rows$fitted_bulk_rna_percent[valid]), finite = TRUE)
      graphics::plot(rows$fitted_bulk_rna_percent[valid], rows$observed_bulk_rna_percent[valid],
        pch = 19, col = colors[[model]], xlab = "Fitted", ylab = "Observed",
        main = sprintf("%s: observed vs fitted", toupper(model)), xlim = limits, ylim = limits)
      graphics::abline(0, 1, lty = 2, col = "grey55")
    } else {
      graphics::plot.new()
      graphics::title(main = sprintf("%s: observed vs fitted", toupper(model)))
      graphics::text(0.5, 0.5, "Fit unavailable", col = "grey40")
    }
    valid_resid = is.finite(rows$fitted_bulk_rna_percent) & is.finite(rows$residual)
    if (any(valid_resid)) {
      graphics::plot(rows$fitted_bulk_rna_percent[valid_resid], rows$residual[valid_resid],
        pch = 19, col = colors[[model]], xlab = "Fitted", ylab = "Residual",
        main = sprintf("%s: residual vs fitted", toupper(model)))
      graphics::abline(h = 0, lty = 2, col = "grey55")
    } else {
      graphics::plot.new()
      graphics::title(main = sprintf("%s: residual vs fitted", toupper(model)))
      graphics::text(0.5, 0.5, "Fit unavailable", col = "grey40")
    }
  }
  # Draw the page title here, while this function owns the nested plotting
  # layout. Callers must not add an outer title after this function returns.
  graphics::mtext(main, side = 3, line = 0.7, outer = TRUE, cex = 1)
  graphics::mtext("Descriptive model diagnostics; not external validation", outer = TRUE, line = -0.2, cex = 0.9)
  invisible(list(file = file, device_opened = device_opened, panels = 6L))
}

plot_qpcr_bulk_rna_model_diagnostics = plot_qpcr_bulk_rna_diagnostics

#' Refit all candidate models while omitting each gRNA in turn.
#'
#' This is a descriptive sensitivity analysis; it is not external validation.
leave_one_grna_out_qpcr_bulk_rna_models = function(fit_or_aggregate, prediction_points = NULL) {
  if (inherits(fit_or_aggregate, "qpcr_bulk_rna_model_fit")) {
    aggregate = structure(list(point_data = fit_or_aggregate$point_data), class = "grna_expression_aggregate")
    if (is.null(prediction_points)) prediction_points = fit_or_aggregate$settings$prediction_points
  } else if (inherits(fit_or_aggregate, "grna_expression_aggregate")) {
    aggregate = fit_or_aggregate
    if (is.null(prediction_points)) prediction_points = 201L
  } else {
    stop("Input must be a qPCR--bulk RNA-seq model fit or aggregate.", call. = FALSE)
  }
  if (!is.numeric(prediction_points) || length(prediction_points) != 1L ||
      prediction_points < 2L || prediction_points != as.integer(prediction_points)) {
    stop("`prediction_points` must be an integer of at least 2.", call. = FALSE)
  }
  grnas = as.character(aggregate$point_data$grna)
  model_names = c("linear", "quadratic", "ll4")
  rows = list()
  # Each refit drops one gRNA point; its metrics quantify influence sensitivity,
  # not held-out prediction performance or external validation.
  for (index in seq_along(grnas)) {
    omitted = grnas[[index]]
    subset_points = aggregate$point_data[aggregate$point_data$grna != omitted, , drop = FALSE]
    subset_aggregate = structure(list(point_data = subset_points), class = "grna_expression_aggregate")
    refit = tryCatch(
      fit_qpcr_bulk_rna_models(subset_aggregate, prediction_points = prediction_points),
      error = function(error) error
    )
    for (model in model_names) {
      if (inherits(refit, "error")) {
        rows[[length(rows) + 1L]] = data.frame(
          omission_index = index, omitted_grna = omitted,
          n_remaining_grna = NROW(subset_points), model = model,
          aic = NA_real_, aic_rank = NA_integer_, rmse = NA_real_,
          converged = FALSE, convergence_message = conditionMessage(refit),
          stringsAsFactors = FALSE
        )
      } else {
        metric = refit$model_table[refit$model_table$model == model, , drop = FALSE]
        rows[[length(rows) + 1L]] = data.frame(
          omission_index = index, omitted_grna = omitted,
          n_remaining_grna = NROW(subset_points), model = model,
          aic = metric$aic[[1L]], aic_rank = metric$aic_rank[[1L]],
          rmse = metric$rmse[[1L]], converged = metric$converged[[1L]],
          convergence_message = metric$convergence_message[[1L]],
          stringsAsFactors = FALSE
        )
      }
    }
  }
  result = do.call(rbind, rows)
  rownames(result) = NULL
  structure(list(
    results = result,
    # Named aliases make the tabular result easy to discover while retaining
    # omission metadata and the prediction setting in the structured object.
    table = result,
    model_table = result,
    omission_metadata = data.frame(
      omission_index = seq_along(grnas), omitted_grna = grnas,
      n_remaining_grna = length(grnas) - 1L, stringsAsFactors = FALSE
    ),
    settings = list(prediction_points = prediction_points, note = "Descriptive sensitivity; not external validation.")
  ), class = "qpcr_bulk_rna_loo_sensitivity")
}

run_qpcr_bulk_rna_leave_one_grna_out = leave_one_grna_out_qpcr_bulk_rna_models

#' Fit qPCR--bulk RNA-seq concordance models for every selected DEG.
#'
#' `data` must be in long form: one row per sample and gene. Each gene is
#' independently checked for four samples per gRNA and eleven gRNA-level
#' points before the three prespecified models are fitted.
fit_qpcr_bulk_rna_deg_batch = function(
    data,
    sample_id_col,
    grna_col,
    qpcr_col,
    gene_col,
    bulk_rna_col,
    genes = NULL,
    expected_replicates = 4L,
    expected_points = 11L,
    prediction_points = 201L) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = F)
  }
  for (column_name in c(sample_id_col, grna_col, qpcr_col, gene_col, bulk_rna_col)) {
    .qpcr_bulk_assert_string(column_name, "column name")
  }
  .qpcr_bulk_validate_columns(
    data,
    c(sample_id_col, grna_col, qpcr_col, gene_col, bulk_rna_col)
  )
  if (!is.numeric(data[[qpcr_col]]) || !is.numeric(data[[bulk_rna_col]])) {
    stop("qPCR and bulk RNA-seq columns must be numeric.", call. = F)
  }
  if (any(is.na(data[[sample_id_col]]) | !nzchar(as.character(data[[sample_id_col]])))) {
    stop("Every row must have a non-missing sample identifier.", call. = F)
  }
  if (any(is.na(data[[gene_col]]) | !nzchar(as.character(data[[gene_col]])))) {
    stop("Every row must have a non-missing gene identifier.", call. = F)
  }

  data[[sample_id_col]] = as.character(data[[sample_id_col]])
  data[[gene_col]] = as.character(data[[gene_col]])
  available_genes = sort(unique(data[[gene_col]]))
  if (is.null(genes)) {
    genes = available_genes
  } else {
    if (!is.character(genes) || length(genes) == 0L || any(is.na(genes)) || any(!nzchar(genes))) {
      stop("`genes` must be NULL or a non-empty character vector.", call. = F)
    }
    genes = unique(genes)
    missing_genes = setdiff(genes, available_genes)
    if (length(missing_genes) > 0L) {
      stop(
        sprintf("Selected gene(s) not found in input: %s", paste(missing_genes, collapse = ", ")),
        call. = F
      )
    }
  }

  fits = list()
  aggregates = list()
  reference_samples = NULL
  # Enforce a shared sample-to-gRNA/qPCR design across genes so differences
  # reflect gene-specific bulk RNA-seq values, not mismatched sample composition.
  for (gene in genes) {
    gene_data = data[data[[gene_col]] == gene, , drop = F]
    if (anyDuplicated(gene_data[[sample_id_col]])) {
      duplicated_samples = unique(gene_data[[sample_id_col]][duplicated(gene_data[[sample_id_col]])])
      stop(
        sprintf(
          "Gene `%s` has duplicate sample identifiers: %s",
          gene,
          paste(duplicated_samples, collapse = ", ")
        ),
        call. = F
      )
    }
    aggregate = aggregate_grna_expression(
      gene_data,
      grna_col = grna_col,
      qpcr_col = qpcr_col,
      bulk_rna_col = bulk_rna_col,
      expected_replicates = expected_replicates,
      expected_points = expected_points
    )
    sample_metadata = gene_data[, c(sample_id_col, grna_col, qpcr_col), drop = F]
    names(sample_metadata) = c("sample_id", "grna", "qpcr_percent")
    if (is.null(reference_samples)) {
      reference_samples = sample_metadata
    } else {
      if (!setequal(reference_samples$sample_id, sample_metadata$sample_id)) {
        stop(
          sprintf("Gene `%s` does not contain the same sample identifiers as the first gene.", gene),
          call. = F
        )
      }
      sample_metadata = sample_metadata[match(reference_samples$sample_id, sample_metadata$sample_id), , drop = F]
      if (any(reference_samples$grna != sample_metadata$grna) ||
          !isTRUE(all.equal(reference_samples$qpcr_percent, sample_metadata$qpcr_percent))) {
        stop(
          sprintf("Gene `%s` has inconsistent gRNA or qPCR values for shared samples.", gene),
          call. = F
        )
      }
    }
    aggregates[[gene]] = aggregate
    fits[[gene]] = fit_qpcr_bulk_rna_models(
      aggregate,
      prediction_points = prediction_points
    )
    fits[[gene]]$leave_one_grna_out = leave_one_grna_out_qpcr_bulk_rna_models(
      fits[[gene]], prediction_points = prediction_points
    )
  }

  structure(
    list(
      genes = genes,
      aggregates = aggregates,
      fits = fits,
      settings = list(
        sample_id_col = sample_id_col,
        grna_col = grna_col,
        qpcr_col = qpcr_col,
        gene_col = gene_col,
        bulk_rna_col = bulk_rna_col,
        expected_replicates = expected_replicates,
        expected_points = expected_points,
        prediction_points = prediction_points
      )
    ),
    class = "qpcr_bulk_rna_deg_batch"
  )
}

#' Extract one gene-labelled result table from a DEG batch fit.
collect_qpcr_bulk_rna_deg_results = function(batch_fit, component) {
  if (!inherits(batch_fit, "qpcr_bulk_rna_deg_batch")) {
    stop("`batch_fit` must be returned by `fit_qpcr_bulk_rna_deg_batch()`.", call. = F)
  }
  .qpcr_bulk_assert_string(component, "component")
  tables = lapply(batch_fit$genes, function(gene) {
    table = batch_fit$fits[[gene]][[component]]
    if (identical(component, "leave_one_grna_out")) {
      table = table$results
    }
    if (!is.data.frame(table)) {
      stop(sprintf("Component `%s` is not a data frame.", component), call. = F)
    }
    table$gene = gene
    table = table[c("gene", setdiff(names(table), "gene"))]
    table
  })
  result = do.call(rbind, tables)
  rownames(result) = NULL
  result
}

#' Write one qPCR--bulk RNA-seq scatterplot page per DEG to a PDF.
plot_qpcr_bulk_rna_deg_batch = function(
    batch_fit,
    file,
    width = 7,
    height = 5.5,
    label_points = T) {
  if (!inherits(batch_fit, "qpcr_bulk_rna_deg_batch")) {
    stop("`batch_fit` must be returned by `fit_qpcr_bulk_rna_deg_batch()`.", call. = F)
  }
  .qpcr_bulk_assert_string(file, "file")
  if (!identical(tolower(tools::file_ext(file)), "pdf")) {
    stop("Batch plots must be written to a `.pdf` file.", call. = F)
  }
  output_directory = dirname(file)
  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = T, showWarnings = F)
  }
  grDevices::pdf(file, width = width, height = height, useDingbats = F)
  on.exit(grDevices::dev.off(), add = T)
  for (gene in batch_fit$genes) {
    plot_qpcr_bulk_rna_models(
      batch_fit$fits[[gene]],
      file = NULL,
      label_points = label_points,
      main = sprintf("%s: qPCR and bulk RNA-seq expression concordance", gene)
    )
  }
  invisible(list(file = file, pages = length(batch_fit$genes)))
}

#' Write a separate multi-page diagnostic PDF, one page per gene.
plot_qpcr_bulk_rna_deg_batch_diagnostics = function(
    batch_fit, file, width = 10, height = 6.5) {
  if (!inherits(batch_fit, "qpcr_bulk_rna_deg_batch")) {
    stop("`batch_fit` must be returned by `fit_qpcr_bulk_rna_deg_batch()`.", call. = FALSE)
  }
  .qpcr_bulk_assert_string(file, "file")
  if (!identical(tolower(tools::file_ext(file)), "pdf")) {
    stop("Batch diagnostic plots must be written to a `.pdf` file.", call. = FALSE)
  }
  output_directory = dirname(file)
  if (!dir.exists(output_directory)) dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(file, width = width, height = height, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (gene in batch_fit$genes) {
    plot_qpcr_bulk_rna_diagnostics(
      batch_fit$fits[[gene]],
      file = NULL,
      main = gene
    )
  }
  invisible(list(file = file, pages = length(batch_fit$genes)))
}

#' Calculate bulk-RNA minus qPCR percentages for every gene by gRNA.
#'
#' The result is a descriptive disagreement table in percentage points.  Gene
#' and gRNA order follows the batch object and is therefore reproducible.
qpcr_bulk_rna_deg_concordance = function(batch_fit) {
  if (!inherits(batch_fit, "qpcr_bulk_rna_deg_batch")) {
    stop("`batch_fit` must be returned by `fit_qpcr_bulk_rna_deg_batch()`.", call. = FALSE)
  }
  genes <- as.character(batch_fit$genes)
  grnas <- sort(unique(as.character(batch_fit$aggregates[[genes[[1L]]]]$point_data$grna)))
  matrix_values <- matrix(NA_real_, nrow = length(genes), ncol = length(grnas),
    dimnames = list(genes, grnas))
  rows <- list()
  for (gene in genes) {
    point_data <- batch_fit$aggregates[[gene]]$point_data
    point_data <- point_data[match(grnas, point_data$grna), , drop = FALSE]
    differences <- point_data$bulk_rna_percent - point_data$qpcr_percent
    matrix_values[gene, ] <- differences
    rows[[length(rows) + 1L]] <- data.frame(
      gene = gene,
      grna = grnas,
      qpcr_percent = as.numeric(point_data$qpcr_percent),
      bulk_rna_percent = as.numeric(point_data$bulk_rna_percent),
      difference_percentage_points = as.numeric(differences),
      stringsAsFactors = FALSE
    )
  }
  long <- do.call(rbind, rows)
  rownames(long) <- NULL
  # Symmetric zero-centred limits give equal visual weight to bulk-RNA excess
  # and deficit across all selected genes (an extreme gene can compress contrast).
  finite_values <- as.numeric(matrix_values)[is.finite(as.numeric(matrix_values))]
  limit <- if (length(finite_values) == 0L) 1 else max(abs(finite_values))
  if (!is.finite(limit) || limit == 0) limit <- 1
  structure(list(
    matrix = matrix_values,
    difference_matrix = matrix_values,
    long = long,
    long_table = long,
    limits = c(-limit, limit),
    symmetric_limits = c(-limit, limit),
    note = "Descriptive disagreement: bulk RNA-seq minus qPCR, in percentage points.",
    genes = genes,
    grnas = grnas
  ), class = "qpcr_bulk_rna_deg_concordance")
}

# Explicit aliases for callers using matrix/heatmap terminology.
compute_qpcr_bulk_rna_deg_concordance = qpcr_bulk_rna_deg_concordance
compute_qpcr_bulk_rna_deg_batch_concordance = qpcr_bulk_rna_deg_concordance
qpcr_bulk_rna_deg_batch_concordance = qpcr_bulk_rna_deg_concordance

#' Plot the gene-by-gRNA percentage-point disagreement matrix.
plot_qpcr_bulk_rna_deg_concordance_heatmap = function(
    batch_fit, file, width = 9, height = 7) {
  concordance <- qpcr_bulk_rna_deg_concordance(batch_fit)
  .qpcr_bulk_assert_string(file, "file")
  if (!identical(tolower(tools::file_ext(file)), "pdf") &&
      !identical(tolower(tools::file_ext(file)), "png")) {
    stop("Heatmap file must end in `.pdf` or `.png`.", call. = FALSE)
  }
  output_directory <- dirname(file)
  if (!dir.exists(output_directory)) dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  extension <- tolower(tools::file_ext(file))
  caller_has_device <- !is.null(grDevices::dev.list())
  old_par <- if (caller_has_device) graphics::par(no.readonly = TRUE) else NULL
  if (extension == "pdf") {
    grDevices::pdf(file, width = width, height = height, useDingbats = FALSE)
  } else {
    grDevices::png(file, width = width, height = height, units = "in", res = 300)
  }
  on.exit({
    if (!is.null(old_par)) {
      graphics::par(old_par)
      graphics::par(cex = old_par$cex)
    }
    grDevices::dev.off()
  }, add = TRUE)
  graphics::layout(matrix(c(1L, 2L), nrow = 1L), widths = c(5.5, 1.2))
  values <- concordance$matrix
  color_count <- 101L
  palette <- grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(color_count)
  graphics::par(mar = c(8, 8, 4, 1))
  # image() maps matrix columns to x and rows bottom-to-top; transpose and
  # reverse genes to align cells with the axis labels and text overlay below.
  graphics::image(seq_len(ncol(values)), seq_len(nrow(values)), t(values[nrow(values):1L, , drop = FALSE]),
    col = palette, zlim = concordance$limits, axes = FALSE,
    xlab = "gRNA", ylab = "Gene",
    main = "Batch concordance disagreement")
  graphics::axis(1L, at = seq_len(ncol(values)), labels = colnames(values), las = 2L, cex.axis = 0.75)
  graphics::axis(2L, at = seq_len(nrow(values)), labels = rev(rownames(values)), las = 2L, cex.axis = 0.8)
  for (row in seq_len(nrow(values))) for (column in seq_len(ncol(values))) {
    value <- values[row, column]
    if (is.finite(value)) graphics::text(column, nrow(values) - row + 1L, sprintf("%.1f", value), cex = 0.55)
  }
  graphics::mtext("Bulk RNA-seq minus qPCR (percentage points); descriptive disagreement, not mechanistic truth.",
    side = 1L, line = 6.5, cex = 0.72, col = "grey30")
  graphics::par(mar = c(8, 1, 4, 4))
  legend_values <- seq(concordance$limits[1L], concordance$limits[2L], length.out = color_count)
  graphics::image(1, legend_values, matrix(legend_values, nrow = 1L), col = palette,
    zlim = concordance$limits, axes = FALSE, xlab = "", ylab = "")
  graphics::axis(4L, at = pretty(concordance$limits), las = 2L)
  graphics::mtext("percentage\npoints", side = 3L, line = 0.5, cex = 0.8)
  invisible(list(file = file, matrix = values, limits = concordance$limits))
}

plot_qpcr_bulk_rna_deg_batch_concordance_heatmap = plot_qpcr_bulk_rna_deg_concordance_heatmap

#' Build a per-gene model support and winner table, retaining failed models.
summarize_qpcr_bulk_rna_deg_models = function(batch_fit) {
  if (!inherits(batch_fit, "qpcr_bulk_rna_deg_batch")) {
    stop("`batch_fit` must be returned by `fit_qpcr_bulk_rna_deg_batch()`.", call. = FALSE)
  }
  rows <- lapply(batch_fit$genes, function(gene) {
    table <- batch_fit$fits[[gene]]$model_table
    valid <- isTRUE(table$converged) & is.finite(table$aic)
    # isTRUE() is scalar, so use vectorized status explicitly for each row.
    valid <- !is.na(table$converged) & table$converged & is.finite(table$aic)
    # Treat numerical near-equality as a tie to avoid arbitrary winners caused
    # solely by floating-point noise; failed models are never eligible.
    minimum <- if (any(valid)) min(table$aic[valid]) else NA_real_
    tolerance <- if (is.finite(minimum)) max(1e-10, abs(minimum) * 1e-10) else NA_real_
    tied <- valid & is.finite(tolerance) & abs(table$aic - minimum) <= tolerance
    winner_models <- if (any(tied)) paste(as.character(table$model)[tied], collapse = "|") else ""
    data.frame(
      gene = gene,
      model = as.character(table$model),
      method = as.character(table$method),
      converged = as.logical(table$converged),
      convergence_status = ifelse(as.logical(table$converged), "converged", "failed"),
      convergence_message = as.character(table$convergence_message),
      aic = as.numeric(table$aic),
      delta_aic = as.numeric(table$delta_aic),
      akaike_weight = as.numeric(table$akaike_weight),
      relative_support = as.numeric(table$akaike_weight),
      rmse = as.numeric(table$rmse),
      aic_rank = as.integer(table$aic_rank),
      aic_winner = tied,
      winner = tied,
      winner_models = winner_models,
      aic_tie = tied & sum(tied) > 1L,
      winner_label = ifelse(!as.logical(table$converged), "failed",
        ifelse(tied & sum(tied) > 1L, "tied_winner",
          ifelse(tied, "winner", "not_winner"))),
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}
collect_qpcr_bulk_rna_deg_model_summary = summarize_qpcr_bulk_rna_deg_models
summarize_qpcr_bulk_rna_deg_batch_models = summarize_qpcr_bulk_rna_deg_models

#' Draw a descriptive per-gene overview of AIC weights and RMSE.
plot_qpcr_bulk_rna_deg_model_overview = function(
    batch_fit, file, width = 8, height = 5.5) {
  summary_table <- summarize_qpcr_bulk_rna_deg_models(batch_fit)
  .qpcr_bulk_assert_string(file, "file")
  if (!identical(tolower(tools::file_ext(file)), "pdf")) {
    stop("Model overview must be written to a `.pdf` file.", call. = FALSE)
  }
  output_directory <- dirname(file)
  if (!dir.exists(output_directory)) dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  caller_has_device <- !is.null(grDevices::dev.list())
  old_par <- if (caller_has_device) graphics::par(no.readonly = TRUE) else NULL
  grDevices::pdf(file, width = width, height = height, useDingbats = FALSE)
  on.exit({
    if (!is.null(old_par)) {
      graphics::par(old_par)
      graphics::par(cex = old_par$cex)
    }
    grDevices::dev.off()
  }, add = TRUE)
  models <- c("linear", "quadratic", "ll4")
  colors <- c(linear = "#0072B2", quadratic = "#009E73", ll4 = "#D55E00")
  # These pages are descriptive within-gene model summaries, not held-out
  # predictive validation or evidence of a biological mechanism.
  for (gene in unique(summary_table$gene)) {
    rows <- summary_table[summary_table$gene == gene, , drop = FALSE]
    rows <- rows[match(models, rows$model), , drop = FALSE]
    weights <- rows$relative_support
    weights[!is.finite(weights)] <- 0
    heights <- max(1, max(weights, na.rm = TRUE) * 1.25)
    graphics::par(mar = c(5, 4.5, 4, 4.5))
    x <- graphics::barplot(weights, names.arg = c("Linear", "Quadratic", "LL.4"),
      col = colors[models], ylim = c(0, heights), ylab = "Akaike weight",
      main = sprintf("%s: model overview", gene))
    graphics::mtext("Relative support is descriptive; models do not establish mechanistic truth.",
      side = 3L, line = 0.4, cex = 0.72, col = "grey35")
    # Label bars while the active user coordinates still represent weights.
    for (index in seq_along(models)) {
      label <- if (!isTRUE(rows$converged[index])) "failed" else if (isTRUE(rows$aic_tie[index])) "tie" else if (isTRUE(rows$aic_winner[index])) "winner" else ""
      if (nzchar(label)) graphics::text(x[index], weights[index] + heights * 0.04, label, cex = 0.75)
    }
    valid_rmse <- is.finite(rows$rmse)
    rmse_max <- if (any(valid_rmse)) max(rows$rmse[valid_rmse]) else 1
    if (!is.finite(rmse_max) || rmse_max <= 0) rmse_max <- 1
    weight_xlim = graphics::par("usr")[1:2]
    graphics::par(new = TRUE)
    graphics::plot(x, rows$rmse, type = "p", pch = 19, col = colors[models],
      axes = FALSE, xlab = "", ylab = "", xlim = weight_xlim, xaxs="i", ylim = c(0, rmse_max * 1.25))
    graphics::axis(4L); graphics::mtext("RMSE", side = 4L, line = 2.5)
    graphics::legend("topright", legend = c("AIC weight", "RMSE"), pch = c(15, 19),
      col = c("grey35", "grey35"), bty = "n", cex = 0.8)
  }
  invisible(list(file = file, pages = length(unique(summary_table$gene))))
}
plot_qpcr_bulk_rna_deg_batch_model_overview = plot_qpcr_bulk_rna_deg_model_overview

print.qpcr_bulk_rna_model_fit = function(x, ...) {
  cat("<qpcr_bulk_rna_model_fit>\n")
  cat(sprintf("gRNA-level points: %d\n", NROW(x$point_data)))
  print(x$model_table, row.names = F)
  invisible(x)
}
