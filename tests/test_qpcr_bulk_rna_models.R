#!/usr/bin/env Rscript

script_argument = grep("^--file=", commandArgs(trailingOnly = F), value = T)
script_path = normalizePath(sub("^--file=", "", script_argument), mustWork = T)
repository_directory = normalizePath(file.path(dirname(script_path), ".."), mustWork = T)
source(file.path(repository_directory, "R", "qpcr_bulk_rna_models.R"))
# drc may open R's default graphics device while reporting an optimizer error;
# keep test-created devices outside the repository working directory.
setwd(tempdir())

# Build the planned 11-gRNA × 4-sample noisy LL.4-like fixture to test
# aggregation, candidate fitting, diagnostics, and sensitivity outputs.
set.seed(20260727)
grna = rep(sprintf("gRNA_%02d", seq_len(11L)), each = 4L)
qpcr_center = rep(c(28, 38, 52, 66, 82, 100, 125, 158, 200, 252, 318), each = 4L)
qpcr_percent = qpcr_center + stats::rnorm(44L, mean = 0, sd = 2)
bulk_rna_percent = 32 + 176 / (1 + exp(-3.7 * (log(qpcr_percent) - log(118)))) +
  stats::rnorm(44L, mean = 0, sd = 4)
sample_data = data.frame(
  sample_id = sprintf("sample_%02d", seq_len(44L)),
  grna = grna,
  qpcr_percent = qpcr_percent,
  bulk_rna_percent = bulk_rna_percent,
  stringsAsFactors = F
)

aggregate = aggregate_grna_expression(
  sample_data,
  grna_col = "grna",
  qpcr_col = "qpcr_percent",
  bulk_rna_col = "bulk_rna_percent"
)
stopifnot(
  inherits(aggregate, "grna_expression_aggregate"),
  NROW(aggregate$sample_data) == 44L,
  NROW(aggregate$point_data) == 11L,
  all(aggregate$point_data$n_samples == 4L),
  all(c("qpcr_mean", "qpcr_sd", "qpcr_se", "qpcr_min", "qpcr_max",
        "bulk_rna_mean", "bulk_rna_sd", "bulk_rna_se", "bulk_rna_min", "bulk_rna_max") %in%
        names(aggregate$point_data)),
  abs(aggregate$point_data$qpcr_mean[1L] - mean(qpcr_percent[1:4])) < 1e-12,
  abs(aggregate$point_data$qpcr_sd[1L] - stats::sd(qpcr_percent[1:4])) < 1e-12,
  abs(aggregate$point_data$bulk_rna_se[1L] - stats::sd(bulk_rna_percent[1:4]) / 2) < 1e-12
)

# Keep a caller-owned device open while exercising both successful fitting and
# leave-one-gRNA-out refits.  drc's transient default devices must not alter the
# caller's device set, and its default output must not leave a PDF behind.
normal_fit_working_directory = tempfile("qpcr-normal-fit-")
dir.create(normal_fit_working_directory, recursive = TRUE)
normal_fit_working_directory_before = getwd()
setwd(normal_fit_working_directory)
normal_caller_device_file = file.path(normal_fit_working_directory, "caller.pdf")
grDevices::pdf(normal_caller_device_file)
normal_devices_before = grDevices::dev.list()
normal_caller_device = as.integer(grDevices::dev.cur())
normal_default_device_before = getOption("device")
options(device = function(...) {
  grDevices::pdf(file.path(normal_fit_working_directory, "drc-default-%d.pdf"))
})
normal_devices_after_fit = NULL
normal_devices_after_loo = NULL
fit = tryCatch(
  {
    result = fit_qpcr_bulk_rna_models(aggregate)
    normal_devices_after_fit <<- grDevices::dev.list()
    loo_result <<- leave_one_grna_out_qpcr_bulk_rna_models(result)
    normal_devices_after_loo <<- grDevices::dev.list()
    result
  },
  finally = {
    # Capture assertions above before disposing of any leaked test devices.
    current_devices = grDevices::dev.list()
    for (device in setdiff(as.integer(current_devices), as.integer(normal_devices_before))) {
      try(grDevices::dev.off(which = device), silent = TRUE)
    }
    if (normal_caller_device %in% as.integer(grDevices::dev.list())) {
      try(grDevices::dev.off(which = normal_caller_device), silent = TRUE)
    }
    options(device = normal_default_device_before)
    setwd(normal_fit_working_directory_before)
  }
)
normal_fit_devices_unchanged = identical(
  as.integer(normal_devices_after_fit), as.integer(normal_devices_before)
)
normal_loo_devices_unchanged = identical(
  as.integer(normal_devices_after_loo), as.integer(normal_devices_before)
)
normal_default_pdfs = list.files(normal_fit_working_directory, pattern = "^drc-default-.*\\.pdf$")
normal_fit_has_no_default_pdf = length(normal_default_pdfs) == 0L
unlink(normal_fit_working_directory, recursive = TRUE)
stopifnot(
  normal_fit_devices_unchanged,
  normal_loo_devices_unchanged,
  normal_fit_has_no_default_pdf
)

stopifnot(
  inherits(fit, "qpcr_bulk_rna_model_fit"),
  identical(fit$model_table$model, c("linear", "quadratic", "ll4")),

  all(fit$model_table$converged),
  all(is.finite(fit$model_table$aic)),
  abs(sum(fit$model_table$akaike_weight) - 1) < 1e-12,
  NROW(fit$coefficients) == 9L,
  NROW(fit$linear_anova) == 2L,
  NROW(fit$quadratic_anova) == 3L,
  NROW(fit$linear_vs_quadratic_anova) == 2L,
  NROW(fit$diagnostics) == 33L,
  identical(names(fit$diagnostics), c(
    "model", "grna", "qpcr_percent", "observed_bulk_rna_percent",
    "fitted_bulk_rna_percent", "residual", "standardized_residual",
    "leverage", "cooks_distance", "converged", "convergence_message"
  )),
  all(is.finite(fit$diagnostics$standardized_residual[fit$diagnostics$model != "ll4"])),
  all(is.na(fit$diagnostics$standardized_residual[fit$diagnostics$model == "ll4"])),
  NROW(loo_result$results) == 33L
)

plot_file = tempfile(fileext = ".pdf")
plot_qpcr_bulk_rna_models(fit, file = plot_file)
stopifnot(file.exists(plot_file), file.info(plot_file)$size > 0)
unlink(plot_file)

# Diagnostic plotting must preserve caller parameters and device sets both
# when it owns a file device and when it draws on the caller device.
diagnostic_caller_file = tempfile(fileext = ".pdf")
grDevices::pdf(diagnostic_caller_file)
diagnostic_par_values = c("mfrow", "mar", "oma", "cex")
graphics::par(mfrow = c(1L, 1L), mar = c(5.1, 4.1, 4.1, 2.1), oma = c(0, 0, 0, 0), cex = 1.07)
diagnostic_par_before_file = graphics::par(diagnostic_par_values)
diagnostic_devices_before_file = as.integer(grDevices::dev.list())
diagnostic_output_file = tempfile(fileext = ".png")
plot_qpcr_bulk_rna_diagnostics(fit, file = diagnostic_output_file)
stopifnot(
  file.exists(diagnostic_output_file), file.info(diagnostic_output_file)$size > 0,
  identical(as.integer(grDevices::dev.list()), diagnostic_devices_before_file),
  isTRUE(all.equal(graphics::par(diagnostic_par_values), diagnostic_par_before_file))
)
unlink(diagnostic_output_file)
diagnostic_par_before_no_file = graphics::par(diagnostic_par_values)
diagnostic_devices_before_no_file = as.integer(grDevices::dev.list())
plot_qpcr_bulk_rna_diagnostics(fit, file = NULL)
stopifnot(
  identical(as.integer(grDevices::dev.list()), diagnostic_devices_before_no_file),
  isTRUE(all.equal(graphics::par(diagnostic_par_values), diagnostic_par_before_no_file))
)
grDevices::dev.off()
unlink(diagnostic_caller_file)

# Extreme, deterministic values exercise the LL.4 optimizer-error path. The
# model-level result must remain inspectable without emitting optimizer output.
ll4_failure_data = data.frame(
  grna = letters[1:11],
  qpcr_percent = c(2.273726e3, 5.704183e-1, 4.024335e-4, 3.563452e-2,
                   8.526349e2, 1.152769e-2, 5.701349, 8.410452e2,
                   1.884045e-4, 5.771833e1, 5.332956e1),
  bulk_rna_percent = c(4.020118e-1, -7.317482e-1, -9.916963e1,
                       -1.000001e6, -4.798441e-2, -9.855884e1,
                       9.898415e1, 1.004120e2, -3.810761e-1,
                       -9.959060e1, 2.688873),
  stringsAsFactors = F
)
# drc opens its default graphics device while reporting this optimizer error.
# Keep that known failure fixture out of the repository working directory.
ll4_failure_working_directory = tempfile("qpcr-ll4-failure-")
dir.create(ll4_failure_working_directory, recursive = TRUE)
ll4_working_directory_before = getwd()
setwd(ll4_failure_working_directory)
ll4_failure_plot_file = file.path(ll4_failure_working_directory, "drc-failure.pdf")
grDevices::pdf(ll4_failure_plot_file)
ll4_default_device_before = getOption("device")
ll4_devices_before = grDevices::dev.list()
ll4_caller_device = as.integer(grDevices::dev.cur())
options(device = function(...) {
  grDevices::pdf(file.path(ll4_failure_working_directory, "drc-default-%d.pdf"))
})
ll4_output_sink_before = sink.number(type = "output")
ll4_message_sink_before = sink.number(type = "message")
ll4_connections_before = showConnections(all = TRUE)
ll4_devices_after = NULL
ll4_default_pdfs = character()
ll4_failure_fit = tryCatch(
  fit_qpcr_bulk_rna_models(
    aggregate_grna_expression(
      ll4_failure_data, "grna", "qpcr_percent", "bulk_rna_percent",
      expected_replicates = NULL, expected_points = NULL
    )
  ),
  finally = {
    ll4_devices_after <<- grDevices::dev.list()
    ll4_default_pdfs <<- list.files(
      ll4_failure_working_directory,
      pattern = "^drc-default-.*\\.pdf$"
    )
    # Dispose of test devices only after recording the exact post-call set.
    current_devices = grDevices::dev.list()
    for (device in setdiff(as.integer(current_devices), as.integer(ll4_devices_before))) {
      try(grDevices::dev.off(which = device), silent = TRUE)
    }
    if (ll4_caller_device %in% as.integer(grDevices::dev.list())) {
      try(grDevices::dev.off(which = ll4_caller_device), silent = TRUE)
    }
    options(device = ll4_default_device_before)
    setwd(ll4_working_directory_before)
  }
)
ll4_devices_unchanged = identical(
  as.integer(ll4_devices_after), as.integer(ll4_devices_before)
)
unlink(ll4_failure_working_directory, recursive = TRUE)
stopifnot(
  ll4_devices_unchanged,
  length(ll4_default_pdfs) == 0L,
  sink.number(type = "output") == ll4_output_sink_before,
  sink.number(type = "message") == ll4_message_sink_before,
  identical(showConnections(all = TRUE), ll4_connections_before),
  !ll4_failure_fit$model_table$converged[ll4_failure_fit$model_table$model == "ll4"],
  nzchar(ll4_failure_fit$model_table$convergence_message[ll4_failure_fit$model_table$model == "ll4"]),
  is.null(ll4_failure_fit$models$ll4)
)

# A singleton gRNA has no estimable within-group SD/SE, which must remain
# descriptive NA values rather than fabricated zero variability.
singleton_data = sample_data
singleton_data = singleton_data[singleton_data$grna != "gRNA_01", , drop = FALSE]
singleton_data = rbind(singleton_data, sample_data[1L, , drop = FALSE])
singleton_aggregate = aggregate_grna_expression(
  singleton_data, "grna", "qpcr_percent", "bulk_rna_percent",
  expected_replicates = NULL, expected_points = NULL
)
stopifnot(
  is.na(singleton_aggregate$point_data$qpcr_sd[singleton_aggregate$point_data$grna == "gRNA_01"]),
  is.na(singleton_aggregate$point_data$qpcr_se[singleton_aggregate$point_data$grna == "gRNA_01"]),
  is.na(singleton_aggregate$point_data$bulk_rna_sd[singleton_aggregate$point_data$grna == "gRNA_01"])
)

# Removing one raw sample must fail the default four-replicate design check.
incorrect_replicates = sample_data[-1L, , drop = F]
replicate_error = tryCatch(
  {
    aggregate_grna_expression(
      incorrect_replicates,
      grna_col = "grna",
      qpcr_col = "qpcr_percent",
      bulk_rna_col = "bulk_rna_percent"
    )
    F
  },
  error = function(error) grepl("must have 4 complete samples", conditionMessage(error), fixed = T)
)
stopifnot(replicate_error)

cat("All qPCR--bulk RNA-seq model tests passed.\n")
