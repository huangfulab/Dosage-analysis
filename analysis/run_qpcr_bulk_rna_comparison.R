#!/usr/bin/env Rscript

usage = function() {
  cat(paste(
    "Usage:",
    "  Rscript analysis/run_qpcr_bulk_rna_comparison.R \\",
    "    --input data/example_qpcr_bulk_rna.csv \\",
    "    --grna-col grna \\",
    "    --qpcr-col qpcr_percent \\",
    "    --bulk-rna-col bulk_rna_percent \\",
    "    --output-dir results/qpcr_bulk_rna_example",
    "",
    "Options:",
    "  --expected-replicates  Samples per gRNA species (default: 4)",
    "  --expected-points      Number of gRNA species (default: 11)",
    "  --prediction-points    Points per fitted curve (default: 201)",
    "  --label-points         true or false (default: true)",
    "  --overwrite            true or false; allow replacing non-empty output (default: false)",
    "  --help                 Show this message",
    sep = "\n"
  ))
}


# Source repository-local code so this run uses a visible implementation whose
# hash is recorded below, rather than an implicit installed-package version.
script_argument = grep("^--file=", commandArgs(trailingOnly = F), value = T)
if (length(script_argument) != 1L) {
  stop("Could not determine the analysis script path.", call. = F)
}
script_path = normalizePath(sub("^--file=", "", script_argument), mustWork = T)
repository_directory = normalizePath(file.path(dirname(script_path), ".."), mustWork = T)
source(file.path(repository_directory, "R", "cli_helpers.R"))
source(file.path(repository_directory, "R", "qpcr_bulk_rna_models.R"))

options = parse_cli_options(
  commandArgs(trailingOnly = T),
  usage,
  c(
    "input", "grna-col", "qpcr-col", "bulk-rna-col", "output-dir",
    "expected-replicates", "expected-points", "prediction-points", "label-points",
    "overwrite"
  )
)
required_options = c("input", "grna-col", "qpcr-col", "bulk-rna-col", "output-dir")
missing_options = setdiff(required_options, names(options))
if (length(missing_options) > 0L) {
  stop(
    sprintf("Missing required option(s): --%s", paste(missing_options, collapse = ", --")),
    call. = F
  )
}

input_file = options[["input"]]
grna_col = options[["grna-col"]]
qpcr_col = options[["qpcr-col"]]
bulk_rna_col = options[["bulk-rna-col"]]
output_directory = options[["output-dir"]]
expected_replicates = parse_cli_positive_integer(
  option_or(options, "expected-replicates", "4"),
  "expected-replicates"
)
expected_points = parse_cli_positive_integer(option_or(options, "expected-points", "11"), "expected-points")
prediction_points = parse_cli_positive_integer(
  option_or(options, "prediction-points", "201"),
  "prediction-points"
)
label_points = parse_cli_boolean(option_or(options, "label-points", "true"), "label-points")
overwrite = parse_cli_boolean(option_or(options, "overwrite", "false"), "overwrite")

if (!file.exists(input_file)) {
  stop(sprintf("Input file does not exist: %s", input_file), call. = F)
}
ensure_output_directory(output_directory, overwrite=overwrite,
  protected_paths=c(input_file, script_path,
    file.path(repository_directory, "R", c("cli_helpers.R", "qpcr_bulk_rna_models.R"))))

message("Reading sample-level qPCR and bulk RNA-seq values...")
sample_data = utils::read.csv(input_file, stringsAsFactors = F, check.names = F)
# Aggregate the planned four sample-level observations to each gRNA before
# fitting, then use leave-one-gRNA-out refits as descriptive sensitivity only.
aggregate = aggregate_grna_expression(
  sample_data,
  grna_col = grna_col,
  qpcr_col = qpcr_col,
  bulk_rna_col = bulk_rna_col,
  expected_replicates = expected_replicates,
  expected_points = expected_points
)
message(sprintf("Fitting three models to %d gRNA-level points...", NROW(aggregate$point_data)))
fit = fit_qpcr_bulk_rna_models(aggregate, prediction_points = prediction_points)
loo = leave_one_grna_out_qpcr_bulk_rna_models(fit, prediction_points = prediction_points)

# Export input/aggregation, candidate evidence, diagnostics/sensitivity, and
# graphics separately so each published summary can be traced to its source.
utils::write.csv(
  aggregate$sample_data,
  file.path(output_directory, "complete_sample_data.csv"),
  row.names = F
)
utils::write.csv(
  aggregate$point_data,
  file.path(output_directory, "grna_level_means.csv"),
  row.names = F
)
utils::write.csv(
  aggregate$point_data,
  file.path(output_directory, "grna_replicate_summaries.csv"),
  row.names = F
)
utils::write.csv(
  fit$diagnostics,
  file.path(output_directory, "model_diagnostics.csv"),
  row.names = F
)
utils::write.csv(
  loo$results,
  file.path(output_directory, "leave_one_grna_out.csv"),
  row.names = F
)
utils::write.csv(
  fit$model_table,
  file.path(output_directory, "model_aic.csv"),
  row.names = F
)
utils::write.csv(
  fit$coefficients,
  file.path(output_directory, "model_coefficients.csv"),
  row.names = F
)
utils::write.csv(
  fit$predictions,
  file.path(output_directory, "model_predictions.csv"),
  row.names = F
)
utils::write.csv(
  fit$linear_anova,
  file.path(output_directory, "linear_anova.csv"),
  row.names = F
)
utils::write.csv(
  fit$quadratic_anova,
  file.path(output_directory, "quadratic_anova.csv"),
  row.names = F
)
utils::write.csv(
  fit$linear_vs_quadratic_anova,
  file.path(output_directory, "linear_vs_quadratic_anova.csv"),
  row.names = F
)
plot_qpcr_bulk_rna_models(
  fit,
  file = file.path(output_directory, "qpcr_bulk_rna_scatterplot.pdf"),
  label_points = label_points
)
plot_qpcr_bulk_rna_diagnostics(
  fit,
  file = file.path(output_directory, "qpcr_bulk_rna_diagnostics.pdf")
)
fit$leave_one_grna_out = loo
saveRDS(fit, file.path(output_directory, "model_objects.rds"), version = 3)

# Record hashes, options, drc version, Git state, and session details; Git
# calls are best-effort so the runner remains usable outside a checkout.
git_commit = tryCatch(
  system2(
    "git",
    c("-C", shQuote(repository_directory), "rev-parse", "HEAD"),
    stdout = T,
    stderr = F
  ),
  error = function(error) NA_character_
)
git_status = tryCatch(
  system2(
    "git",
    c("-C", shQuote(repository_directory), "status", "--porcelain"),
    stdout = T,
    stderr = F
  ),
  error = function(error) NA_character_
)
runner_source = file.path(repository_directory, "analysis", basename(script_path))
model_source = file.path(repository_directory, "R", "qpcr_bulk_rna_models.R")
cli_source = file.path(repository_directory, "R", "cli_helpers.R")
manifest = data.frame(
  key = c(
    "run_time_utc", "input", "input_md5", "grna_col", "qpcr_col", "bulk_rna_col",
    "expected_replicates", "expected_points", "prediction_points", "label_points", "overwrite",
    "drc_version", "runner_source", "runner_source_md5", "model_source", "model_source_md5",
    "cli_helpers_source", "cli_helpers_source_md5", "git_commit", "git_worktree_dirty"
  ),
  value = c(
    format(Sys.time(), tz = "UTC", usetz = T),
    input_file,
    file_md5(input_file),
    grna_col,
    qpcr_col,
    bulk_rna_col,
    expected_replicates,
    expected_points,
    prediction_points,
    label_points,
    overwrite,
    fit$settings$drc_version,
    runner_source,
    file_md5(runner_source),
    model_source,
    file_md5(model_source),
    cli_source,
    file_md5(cli_source),
    paste(git_commit, collapse = ""),
    length(git_status) > 0L && any(nzchar(git_status))
  ),
  stringsAsFactors = F
)
utils::write.table(
  manifest,
  file.path(output_directory, "run_manifest.tsv"),
  sep = "\t",
  row.names = F,
  quote = T,
  qmethod = "double"
)
session_report = write_base_session_info(file.path(output_directory, "session_info.txt"))
# Preserve nested inputs, settings, source fingerprints, and an inventory of
# every generated artifact; TSV remains for quick inspection, JSON for reuse.
write_analysis_provenance(
  file.path(output_directory, "analysis_provenance.json"),
  workflow = "qpcr_bulk_rna_comparison",
  inputs = list(input = list(path = input_file, md5 = file_md5(input_file))),
  settings = list(
    grna_col = grna_col, qpcr_col = qpcr_col, bulk_rna_col = bulk_rna_col,
    expected_replicates = expected_replicates, expected_points = expected_points,
    prediction_points = prediction_points, label_points = label_points
  ),
  sources = list(
    runner = list(path = runner_source, md5 = file_md5(runner_source)),
    models = list(path = model_source, md5 = file_md5(model_source)),
    cli_helpers = list(path = cli_source, md5 = file_md5(cli_source)),
    git_commit = paste(git_commit, collapse = ""),
    git_worktree_dirty = length(git_status) > 0L && any(nzchar(git_status))
  ),
  summary = list(
    input_rows = NROW(sample_data), grna_points = NROW(aggregate$point_data),
    models = as.list(fit$model_table$model), converged_models = sum(fit$model_table$converged)
  ),
  output_directory = output_directory
)

cat(sprintf(
  "Completed qPCR--bulk RNA-seq comparison for %d gRNA-level points. Results are in %s\n",
  NROW(aggregate$point_data),
  output_directory
))
cat(paste(session_report, collapse = "\n"), "\n")
