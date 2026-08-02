#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument), mustWork = TRUE)
repository_directory <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(repository_directory, "R", "dosage_models.R"))

# Construct deterministic positive-dose LL.4-like and quadratic endpoints with
# controls near dose 1, giving normalization and candidate fits known targets.
set.seed(20260721)
dose <- c(seq(0.2, 0.95, length.out = 32L), rep(1, 8L))
group <- c(rep("Targeted", 32L), rep("Controls", 8L))
hill_truth <- 0.4 + 1.4 / (1 + exp(-5 * (log(dose) - log(0.6))))
quadratic_truth <- 0.7 + 0.2 * dose + 1.1 * dose^2
test_data <- data.frame(
  group = group,
  dose = dose,
  hill = hill_truth + stats::rnorm(length(dose), sd = 0.02),
  quadratic = quadratic_truth + stats::rnorm(length(dose), sd = 0.02),
  stringsAsFactors = FALSE
)

# Verify the public fitting/normalization contract before testing each family.
prepared <- prepare_dosage_data(
  test_data,
  dose_col = "dose",
  response_cols = c("hill", "quadratic"),
  group_col = "group",
  control_label = "Controls",
  normalize_dose = TRUE,
  normalize_response = TRUE
)
stopifnot(
  abs(mean(prepared$data$dose[prepared$data$group == "Controls"]) - 1) < 1e-12,
  abs(mean(prepared$data$hill[prepared$data$group == "Controls"]) - 1) < 1e-12,
  NROW(prepared$normalization) == 3L
)

polynomial <- runmodelsaic_poly(
  prepared$data,
  dose_col = "dose",
  response_col = "quadratic",
  degrees = 1:2,
  group_col = "group"
)
stopifnot(
  inherits(polynomial, "dosage_model_fit"),
  NROW(polynomial$model_table) == 2L,
  all(is.finite(polynomial$model_table$aic)),
  polynomial$best_model %in% polynomial$model_table$model,
  NROW(polynomial$nested_tests) == 1L
)

high_degree_error <- tryCatch(
  {
    runmodelsaic_poly(
      prepared$data,
      dose_col = "dose",
      response_col = "quadratic",
      degrees = 3,
      group_col = "group"
    )
    FALSE
  },
  error = function(error) grepl("Only degrees 1", conditionMessage(error), fixed = TRUE)
)
stopifnot(high_degree_error)

hill <- runmodelsaic_drc(
  prepared$data,
  dose_col = "dose",
  response_col = "hill",
  group_col = "group"
)
stopifnot(
  inherits(hill, "dosage_model_fit"),
  NROW(hill$model_table) == 1L,
  hill$model_table$model == "ll4",
  hill$model_table$converged,
  is.character(hill$best_model),
  any(grepl("e:", hill$coefficients$term, fixed = TRUE))
)

spline <- runmodelsaic_spline(
  prepared$data,
  dose_col = "dose",
  response_col = "hill",
  df = c(3, 4),
  group_col = "group"
)
stopifnot(
  inherits(spline, "dosage_model_fit"),
  NROW(spline$model_table) == 2L,
  all(is.finite(spline$model_table$aic))
)

comparison <- compare_dosage_models(polynomial, hill, spline)
stopifnot(
  all(diff(comparison$aic) >= 0),
  abs(sum(comparison$akaike_weight) - 1) < 1e-12,
  identical(comparison$aic_rank, seq_len(NROW(comparison)))
)

diagnostics <- diagnose_dosage_models(list(polynomial, hill, spline))
expected_diagnostic_columns <- c(
  "response", "model", "family", "row_id", "original_row_id", "dose", "observed",
  "fitted", "residual", "standardized_residual", "leverage", "cooks_distance",
  "converged", "convergence_status", "convergence_message"
)
stopifnot(
  identical(names(diagnostics), expected_diagnostic_columns),
  NROW(diagnostics) == NROW(prepared$data) *
    sum(vapply(list(polynomial, hill, spline), function(x) NROW(x$models), integer(1))),
  all(diagnostics$original_row_id == diagnostics$row_id),
  all.equal(diagnostics$residual,
    diagnostics$observed - diagnostics$fitted, check.attributes = FALSE) == TRUE,
  all(is.finite(diagnostics$standardized_residual[diagnostics$family == "polynomial"])),
  all(is.finite(diagnostics$leverage[diagnostics$family == "polynomial"])),
  all(is.na(diagnostics$standardized_residual[diagnostics$family != "polynomial"])),
  all(is.na(diagnostics$cooks_distance[diagnostics$family != "polynomial"]))
)
diagnostic_plot_file <- tempfile(fileext = ".pdf")
diagnostic_plot <- plot_dosage_model_diagnostics(
  list(polynomial, hill, spline), file = diagnostic_plot_file
)
stopifnot(
  diagnostic_plot$panels == 2L * NROW(unique(diagnostics$model)),
  file.exists(diagnostic_plot_file), file.info(diagnostic_plot_file)$size > 0
)
unlink(diagnostic_plot_file)

plot_file <- tempfile(fileext = ".pdf")
plot_result <- plot_dosage_models(
  prepared$data,
  dose_col = "dose",
  response_col = "hill",
  fits = list(polynomial, hill, spline),
  group_col = "group",
  file = plot_file
)
stopifnot(
  file.exists(plot_file),
  file.info(plot_file)$size > 0,
  identical(sort(plot_result$plotted_methods), c("drc_ll4", "polynomial"))
)
unlink(plot_file)

# Deliberately extreme values force the LL.4 failure path and test that its
# sink, connection, and graphics-device cleanup preserves caller resources.
failed_ll4_data <- data.frame(
  dose = exp(seq(log(1e-9), log(1e9), length.out = 6L)),
  response = c(-1e300, 1e300, -1e300, 1e300, -1e300, 1e300)
)
sink_output_before <- sink.number(type = "output")
sink_message_before <- sink.number(type = "message")
connections_before <- rownames(showConnections(all = TRUE))
caller_device_file <- tempfile(fileext = ".pdf")
grDevices::pdf(caller_device_file)
devices_before_failed_fit <- grDevices::dev.list()
device_before_failed_fit <- grDevices::dev.cur()
failed_ll4 <- suppressWarnings(runmodelsaic_drc(
  failed_ll4_data,
  dose_col = "dose",
  response_col = "response"
))
device_after_failed_fit <- grDevices::dev.cur()
devices_after_failed_fit <- grDevices::dev.list()
connections_after <- rownames(showConnections(all = TRUE))
stopifnot(
  inherits(failed_ll4, "dosage_model_fit"),
  identical(failed_ll4$method, "drc_ll4"),
  identical(failed_ll4$best_model, NA_character_),
  !failed_ll4$model_table$converged,
  nzchar(failed_ll4$model_table$convergence_message),
  identical(sink.number(type = "output"), sink_output_before),
  identical(sink.number(type = "message"), sink_message_before),
  identical(connections_after, connections_before),
  identical(devices_after_failed_fit, devices_before_failed_fit),
  identical(device_after_failed_fit, device_before_failed_fit),
  !file.exists("Rplots.pdf")
)
grDevices::dev.off()
unlink(caller_device_file)
failed_plot_file <- tempfile(fileext = ".pdf")
device_before_failed_plot <- grDevices::dev.cur()
failed_plot <- plot_dosage_models(
  failed_ll4_data,
  dose_col = "dose",
  response_col = "response",
  fits = failed_ll4,
  file = failed_plot_file
)
device_after_failed_plot <- grDevices::dev.cur()
stopifnot(
  file.exists(failed_plot_file),
  file.info(failed_plot_file)$size > 0,
  length(failed_plot$plotted_methods) == 0L,
  identical(device_before_failed_plot, device_after_failed_plot),
  !file.exists("Rplots.pdf")
)
unlink(failed_plot_file)

nonpositive_data <- prepared$data
nonpositive_data$dose[1L] <- 0
positive_dose_error <- tryCatch(
  {
    runmodelsaic_drc(nonpositive_data, "dose", "hill")
    FALSE
  },
  error = function(error) grepl("greater than zero", conditionMessage(error), fixed = TRUE)
)
stopifnot(positive_dose_error)

cat("All dosage-model tests passed.\n")
