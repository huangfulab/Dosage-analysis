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

# Returned polynomial objects must predict the same curves after the loop and
# after serialization, rather than reevaluating the final loop degree.
polynomial_rds = tempfile(fileext=".rds")
saveRDS(polynomial, polynomial_rds)
for (retained in list(polynomial, readRDS(polynomial_rds))) {
  for (model in names(retained$models)) {
    expected = retained$predictions[retained$predictions$model == model, ]
    actual = as.numeric(stats::predict(retained$models[[model]],
      newdata=data.frame(.dose=expected$dose)))
    stopifnot(isTRUE(all.equal(actual, expected$predicted)))
  }
}
unlink(polynomial_rds)
# A response offset must not erase a real difference in residual likelihood.
for (offset in c(0, 1e9)) {
  offset_fit = suppressWarnings(runmodelsaic_poly(
    data.frame(dose=1:10, response=offset + (1:10)^2), "dose", "response"))
  stopifnot(identical(offset_fit$best_model, "polynomial_degree_2"))
}
for (offset in c(1e6, 1e9)) {
  singular = tryCatch(runmodelsaic_poly(
    data.frame(dose=offset + 1:10, response=(1:10)^2), "dose", "response"), error=identity)
  stopifnot(inherits(singular, "error"), grepl("rank deficient", conditionMessage(singular), fixed=T))
}
exact = .gaussian_fit_statistics(c(0, 1e-160), c(0, 1e-160), 1)
imperfect = .gaussian_fit_statistics(c(0, 1e-160), c(1e-160, 1e-160), 1)
stopifnot(is.finite(exact$aic), exact$aic < imperfect$aic)
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
  hill$model_table$mean_parameters == 2L,
  is.character(hill$best_model),
  any(grepl("e:", hill$coefficients$term, fixed = TRUE)),
  any(hill$coefficients$term == "c:(Intercept)" & hill$coefficients$fixed),
  any(hill$coefficients$term == "d:(Intercept)" & hill$coefficients$fixed),
  isTRUE(all.equal(
    hill$coefficients$estimate[hill$coefficients$term == "c:(Intercept)"],
    min(prepared$data$hill)
  )),
  isTRUE(all.equal(
    hill$coefficients$estimate[hill$coefficients$term == "d:(Intercept)"],
    max(prepared$data$hill)
  )),
  isTRUE(all.equal(unname(hill$settings$fixed_parameters[["c"]]), min(prepared$data$hill))),
  isTRUE(all.equal(unname(hill$settings$fixed_parameters[["d"]]), max(prepared$data$hill)))
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

different_response_error = tryCatch(compare_dosage_models(polynomial, hill, spline), error=identity)
stopifnot(inherits(different_response_error, "error"),
  grepl("same response and observations", conditionMessage(different_response_error), fixed=T),
  inherits(tryCatch(plot_dosage_model_diagnostics(list(polynomial, hill)), error=identity), "error"))
# Positive comparisons/plots must use the same endpoint, unlike the old fixture.
polynomial = runmodelsaic_poly(prepared$data, "dose", "hill", group_col="group")
omitted_data = prepared$data
omitted_data$hill[[1L]] = NA_real_
omitted_fit = runmodelsaic_poly(omitted_data, "dose", "hill", na_action="omit")
omitted_error = tryCatch(compare_dosage_models(polynomial, omitted_fit), error=identity)
stopifnot(inherits(omitted_error, "error"))
omitted_pdf = tempfile(fileext=".pdf")
plot_dosage_models(omitted_data, "dose", "hill", omitted_fit, file=omitted_pdf, na_action="omit")
stopifnot(file.exists(omitted_pdf))
unlink(omitted_pdf)
for (degree in 1:2) {
  small_data = data.frame(dose=seq_len(degree + 2L), response=seq_len(degree + 2L)^2 + rep(c(0.1, -0.2), length.out=degree + 2L))
  small_fit = runmodelsaic_poly(small_data, "dose", "response", degrees=degree)
  small_pdf = tempfile(fileext=".pdf")
  plot_dosage_models(small_data, "dose", "response", small_fit, file=small_pdf)
  stopifnot(file.exists(small_pdf))
  unlink(small_pdf)
}
stopifnot(inherits(tryCatch(plot_dosage_models(prepared$data, "dose", "hill",
  list(polynomial, omitted_fit)), error=identity), "error"))
changed_data = prepared$data
changed_data$hill[[1L]] = changed_data$hill[[1L]] + 0.1
changed_fit = runmodelsaic_poly(changed_data, "dose", "hill")
stopifnot(inherits(tryCatch(compare_dosage_models(polynomial, changed_fit), error=identity), "error"))
close_df = c(3.00000001, 3.00000002)
close_spline = runmodelsaic_spline(prepared$data, "dose", "hill", df=close_df)
stopifnot(length(close_spline$models) == 2L, NROW(close_spline$model_table) == 2L,
  identical(close_spline$model_table$requested_df, close_df),
  length(unique(close_spline$predictions$model)) == 2L)
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
# Restoring layout can reset cex: verify both owned-device and caller-device paths.
caller_diagnostics_file = tempfile(fileext=".pdf")
grDevices::pdf(caller_diagnostics_file)
graphics::par(cex=1.07)
par_names = c("mfrow", "mar", "oma", "cex")
caller_par = graphics::par(par_names)
plot_dosage_model_diagnostics(polynomial, file=NULL)
stopifnot(isTRUE(all.equal(caller_par, graphics::par(par_names))))
plot_dosage_model_diagnostics(polynomial, file=diagnostic_plot_file)
stopifnot(isTRUE(all.equal(caller_par, graphics::par(par_names))))
grDevices::dev.off()
unlink(c(caller_diagnostics_file, diagnostic_plot_file))

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
