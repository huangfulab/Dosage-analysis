#!/usr/bin/env Rscript

usage <- function() {
  cat(paste(
    "Usage:",
    "  Rscript analysis/run_dosage_analysis.R \\",
    "    --input data/example_dosage.csv \\",
    "    --dose-col dose \\",
    "    --response-cols response_hill,response_quadratic \\",
    "    --group-col group \\",
    "    --control-label Controls \\",
    "    --normalize-dose true \\",
    "    --normalize-response true \\",
    "    --output-dir results/example",
    "",
    "Options:",
    "  --poly-degrees       Must be 1,2 (default: 1,2)",
    "  --spline-df           Comma-separated smoothing-spline df (default: 3)",
    "  --prediction-points   Number of curve points (default: 201)",
    "  --log-dose-plot       true or false (default: false)",
    "  --overwrite           true or false; allow replacing non-empty output (default: false)",
    "  --help                Show this message",
    sep = "\n"
  ))
}


safe_file_stem <- function(value) {
  result <- gsub("[^A-Za-z0-9._-]+", "_", value)
  result <- gsub("^_+|_+$", "", result)
  if (!nzchar(result)) "response" else result
}

# Resolve and source repository-local code, rather than relying on an installed
# package, so the manifest can identify the exact implementation that was run.
script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_argument) != 1L) {
  stop("Could not determine the analysis script path.", call. = FALSE)
}
script_path <- normalizePath(sub("^--file=", "", script_argument), mustWork = TRUE)
repository_directory <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(repository_directory, "R", "cli_helpers.R"))
source(file.path(repository_directory, "R", "dosage_models.R"))

options <- parse_cli_options(
  commandArgs(trailingOnly = TRUE),
  usage,
  c(
    "input", "dose-col", "response-cols", "group-col", "control-label",
    "normalize-dose", "normalize-response", "output-dir", "poly-degrees",
    "spline-df", "prediction-points", "log-dose-plot", "overwrite"
  )
)
required_options <- c("input", "dose-col", "response-cols", "output-dir")
missing_options <- setdiff(required_options, names(options))
if (length(missing_options) > 0L) {
  stop(
    sprintf("Missing required option(s): --%s", paste(missing_options, collapse = ", --")),
    call. = FALSE
  )
}

input_file <- options[["input"]]
dose_col <- options[["dose-col"]]
response_cols <- strsplit(options[["response-cols"]], ",", fixed = TRUE)[[1L]]
response_cols <- response_cols[nzchar(response_cols)]
output_directory <- options[["output-dir"]]
group_col <- option_or(options, "group-col")
control_label <- option_or(options, "control-label")
normalize_dose <- parse_cli_boolean(option_or(options, "normalize-dose", "false"), "normalize-dose")
normalize_response <- parse_cli_boolean(
  option_or(options, "normalize-response", "false"),
  "normalize-response"
)
poly_degrees <- parse_cli_numeric_list(
  option_or(options, "poly-degrees", "1,2"),
  "poly-degrees",
  integer = TRUE
)
if (!identical(poly_degrees, c(1L, 2L))) {
  stop("`--poly-degrees` must be exactly `1,2` for this publication workflow.", call. = FALSE)
}
spline_df <- parse_cli_numeric_list(option_or(options, "spline-df", "3"), "spline-df")
prediction_points <- parse_cli_numeric_list(
  option_or(options, "prediction-points", "201"),
  "prediction-points",
  integer = TRUE
)
if (length(prediction_points) != 1L || prediction_points < 2L) {
  stop("`--prediction-points` must be one integer of at least 2.", call. = FALSE)
}
log_dose_plot <- parse_cli_boolean(option_or(options, "log-dose-plot", "false"), "log-dose-plot")
overwrite <- parse_cli_boolean(option_or(options, "overwrite", "false"), "overwrite")

if (!file.exists(input_file)) {
  stop(sprintf("Input file does not exist: %s", input_file), call. = FALSE)
}
ensure_output_directory(output_directory, overwrite = overwrite)

raw_data <- utils::read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
prepared <- prepare_dosage_data(
  raw_data,
  dose_col = dose_col,
  response_cols = response_cols,
  group_col = group_col,
  control_label = control_label,
  normalize_dose = normalize_dose,
  normalize_response = normalize_response
)
analysis_data <- prepared$data

utils::write.csv(
  analysis_data,
  file.path(output_directory, "prepared_data.csv"),
  row.names = FALSE
)
utils::write.csv(
  prepared$normalization,
  file.path(output_directory, "normalization_factors.csv"),
  row.names = FALSE
)

all_results <- list()
all_comparisons <- list()
all_metrics <- list()
all_coefficients <- list()
all_predictions <- list()
all_nested_tests <- list()
all_diagnostics <- list()

# For each endpoint, compare prespecified simple, saturating, and flexible
# descriptive candidates on identical normalized observations; inspect
# diagnostics alongside AIC rather than treating selection as mechanism proof.
for (response_col in response_cols) {
  message(sprintf("Fitting models for `%s`...", response_col))
  polynomial <- runmodelsaic_poly(
    analysis_data,
    dose_col = dose_col,
    response_col = response_col,
    degrees = poly_degrees,
    group_col = group_col,
    prediction_points = prediction_points
  )
  hill <- runmodelsaic_drc(
    analysis_data,
    dose_col = dose_col,
    response_col = response_col,
    group_col = group_col,
    prediction_points = prediction_points
  )
  spline <- runmodelsaic_spline(
    analysis_data,
    dose_col = dose_col,
    response_col = response_col,
    df = spline_df,
    group_col = group_col,
    prediction_points = prediction_points
  )
  fits <- list(polynomial = polynomial, hill = hill, spline = spline)
  comparison <- compare_dosage_models(fits)
  comparison$response <- response_col

  metrics <- .rbind_fill(lapply(fits, `[[`, "model_table"))
  metrics$response <- response_col
  coefficients <- .rbind_fill(lapply(fits, `[[`, "coefficients"))
  if (NROW(coefficients) > 0L) coefficients$response <- response_col
  predictions <- .rbind_fill(lapply(fits, `[[`, "predictions"))
  predictions$response <- response_col
  nested_tests <- polynomial$nested_tests
  if (NROW(nested_tests) > 0L) nested_tests$response <- response_col
  diagnostics <- diagnose_dosage_models(fits)
  # Keep one explicit, stable column order in both per-response and combined
  # diagnostic CSVs, including rows for failed candidates.
  diagnostic_columns <- c(
    "response", "model", "family", "row_id", "original_row_id", "dose", "observed",
    "fitted", "residual", "standardized_residual", "leverage", "cooks_distance",
    "converged", "convergence_status", "convergence_message"
  )
  diagnostics <- diagnostics[diagnostic_columns]
  diagnostic_csv <- file.path(
    output_directory,
    paste0("model_diagnostics_", safe_file_stem(response_col), ".csv")
  )
  utils::write.csv(diagnostics, diagnostic_csv, row.names = FALSE)

  plot_file <- file.path(
    output_directory,
    paste0("model_comparison_", safe_file_stem(response_col), ".pdf")
  )
  plot_dosage_models(
    analysis_data,
    dose_col = dose_col,
    response_col = response_col,
    fits = fits,
    group_col = group_col,
    file = plot_file,
    log_dose = log_dose_plot
  )
  diagnostic_plot_file <- file.path(
    output_directory,
    paste0("model_diagnostics_", safe_file_stem(response_col), ".pdf")
  )
  plot_dosage_model_diagnostics(
    fits,
    file = diagnostic_plot_file,
    main = sprintf("%s: dosage model diagnostics", response_col)
  )

  all_results[[response_col]] <- fits
  all_comparisons[[response_col]] <- comparison
  all_metrics[[response_col]] <- metrics
  all_coefficients[[response_col]] <- coefficients
  all_predictions[[response_col]] <- predictions
  all_nested_tests[[response_col]] <- nested_tests
  all_diagnostics[[response_col]] <- diagnostics
}

# Write combined machine-readable outputs after retaining per-response files
# so downstream manuscript tables can trace every endpoint-specific fit.
utils::write.csv(
  .rbind_fill(all_comparisons),
  file.path(output_directory, "model_comparison.csv"),
  row.names = FALSE
)
utils::write.csv(
  .rbind_fill(all_metrics),
  file.path(output_directory, "model_metrics.csv"),
  row.names = FALSE
)
utils::write.csv(
  .rbind_fill(all_coefficients),
  file.path(output_directory, "model_coefficients.csv"),
  row.names = FALSE
)
utils::write.csv(
  .rbind_fill(all_predictions),
  file.path(output_directory, "model_predictions.csv"),
  row.names = FALSE
)
utils::write.csv(
  .rbind_fill(all_nested_tests),
  file.path(output_directory, "polynomial_nested_tests.csv"),
  row.names = FALSE
)
utils::write.csv(
  .rbind_fill(all_diagnostics),
  file.path(output_directory, "model_diagnostics.csv"),
  row.names = FALSE
)
saveRDS(all_results, file.path(output_directory, "model_objects.rds"), version = 3)

# Git metadata is best-effort so archived source trees still run; hashes,
# options, UTC time, and session info preserve the reproducibility record.
git_commit <- tryCatch(
  system2(
    "git",
    c("-C", shQuote(repository_directory), "rev-parse", "HEAD"),
    stdout = TRUE,
    stderr = FALSE
  ),
  error = function(error) NA_character_
)
git_status <- tryCatch(
  system2(
    "git",
    c("-C", shQuote(repository_directory), "status", "--porcelain"),
    stdout = TRUE,
    stderr = FALSE
  ),
  error = function(error) NA_character_
)
runner_source <- file.path(repository_directory, "analysis", basename(script_path))
model_source <- file.path(repository_directory, "R", "dosage_models.R")
cli_source <- file.path(repository_directory, "R", "cli_helpers.R")
manifest <- data.frame(
  key = c(
    "run_time_utc", "input", "input_md5", "dose_col", "response_cols", "group_col",
    "control_label", "normalize_dose", "normalize_response", "poly_degrees",
    "spline_df", "prediction_points", "log_dose_plot", "overwrite",
    "runner_source", "runner_source_md5", "model_source", "model_source_md5",
    "cli_helpers_source", "cli_helpers_source_md5", "git_commit", "git_worktree_dirty"
  ),
  value = c(
    format(Sys.time(), tz = "UTC", usetz = TRUE),
    input_file,
    file_md5(input_file),
    dose_col,
    paste(response_cols, collapse = ","),
    if (is.null(group_col)) "" else group_col,
    if (is.null(control_label)) "" else control_label,
    normalize_dose,
    normalize_response,
    paste(poly_degrees, collapse = ","),
    paste(spline_df, collapse = ","),
    prediction_points,
    log_dose_plot,
    overwrite,
    runner_source,
    file_md5(runner_source),
    model_source,
    file_md5(model_source),
    cli_source,
    file_md5(cli_source),
    paste(git_commit, collapse = ""),
    length(git_status) > 0L && any(nzchar(git_status))
  ),
  stringsAsFactors = FALSE
)
utils::write.table(
  manifest,
  file.path(output_directory, "run_manifest.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
session_report <- write_base_session_info(file.path(output_directory, "session_info.txt"))
# Preserve nested inputs, settings, source fingerprints, and an inventory of
# every generated artifact; TSV remains for quick inspection, JSON for reuse.
write_analysis_provenance(
  file.path(output_directory, "analysis_provenance.json"),
  workflow = "dosage_analysis",
  inputs = list(input = list(path = input_file, md5 = file_md5(input_file))),
  settings = list(
    dose_col = dose_col, response_cols = as.list(response_cols), group_col = group_col,
    control_label = control_label, normalize_dose = normalize_dose,
    normalize_response = normalize_response, poly_degrees = as.list(poly_degrees),
    spline_df = as.list(spline_df), prediction_points = prediction_points,
    log_dose_plot = log_dose_plot
  ),
  sources = list(
    runner = list(path = runner_source, md5 = file_md5(runner_source)),
    models = list(path = model_source, md5 = file_md5(model_source)),
    cli_helpers = list(path = cli_source, md5 = file_md5(cli_source)),
    git_commit = paste(git_commit, collapse = ""),
    git_worktree_dirty = length(git_status) > 0L && any(nzchar(git_status))
  ),
  summary = list(response_count = length(response_cols), input_rows = NROW(raw_data)),
  output_directory = output_directory
)

cat(sprintf(
  "Completed %d response(s). Results are in %s\n",
  length(response_cols),
  output_directory
))
cat(paste(session_report, collapse = "\n"), "\n")
