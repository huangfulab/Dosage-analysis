#!/usr/bin/env Rscript

usage = function() {
  cat(paste(
    "Usage:",
    "  Rscript analysis/run_qpcr_bulk_rna_deg_batch.R \\",
    "    --input data/deg_expression_long.csv \\",
    "    --sample-id-col sample_id \\",
    "    --grna-col grna \\",
    "    --qpcr-col qpcr_percent \\",
    "    --gene-col gene \\",
    "    --bulk-rna-col bulk_rna_percent \\",
    "    --output-dir results/qpcr_bulk_rna_degs",
    "",
    "Input must be long-form: one row per sample and gene.",
    "",
    "Options:",
    "  --deg-file             Optional CSV listing genes to analyse",
    "  --deg-gene-col         Gene column in --deg-file (default: gene)",
    "  --expected-replicates  Samples per gRNA per gene (default: 4)",
    "  --expected-points      gRNA-level points per gene (default: 11)",
    "  --prediction-points    Points per fitted curve (default: 201)",
    "  --label-points         true or false (default: true)",
    "  --overwrite            true or false; allow replacing non-empty output (default: false)",
    "  --help                 Show this message",
    "Outputs include descriptive disagreement heatmap/matrix and model-support overview.",
    sep = "\n"
  ))
}

# Flatten gene-specific gRNA summaries while retaining the batch gene ordering
# used in fitting, plots, and all downstream machine-readable output.
collect_grna_means = function(batch_fit) {
  tables = lapply(batch_fit$genes, function(gene) {
    table = batch_fit$aggregates[[gene]]$point_data
    table$gene = gene
    table[c("gene", setdiff(names(table), "gene"))]
  })
  result = do.call(rbind, tables)
  rownames(result) = NULL
  result
}

# Source the repository-local implementation so the manifest can hash exactly
# what executed rather than depending on an installed package.
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
    "input", "sample-id-col", "grna-col", "qpcr-col", "gene-col", "bulk-rna-col",
    "output-dir", "deg-file", "deg-gene-col", "expected-replicates", "expected-points",
    "prediction-points", "label-points", "overwrite"
  )
)
required_options = c(
  "input", "sample-id-col", "grna-col", "qpcr-col", "gene-col", "bulk-rna-col", "output-dir"
)
missing_options = setdiff(required_options, names(options))
if (length(missing_options) > 0L) {
  stop(
    sprintf("Missing required option(s): --%s", paste(missing_options, collapse = ", --")),
    call. = F
  )
}

input_file = options[["input"]]
sample_id_col = options[["sample-id-col"]]
grna_col = options[["grna-col"]]
qpcr_col = options[["qpcr-col"]]
gene_col = options[["gene-col"]]
bulk_rna_col = options[["bulk-rna-col"]]
output_directory = options[["output-dir"]]
deg_file = option_or(options, "deg-file")
deg_gene_col = option_or(options, "deg-gene-col", "gene")
expected_replicates = parse_cli_positive_integer(option_or(options, "expected-replicates", "4"), "expected-replicates")
expected_points = parse_cli_positive_integer(option_or(options, "expected-points", "11"), "expected-points")
prediction_points = parse_cli_positive_integer(option_or(options, "prediction-points", "201"), "prediction-points")
label_points = parse_cli_boolean(option_or(options, "label-points", "true"), "label-points")
overwrite = parse_cli_boolean(option_or(options, "overwrite", "false"), "overwrite")

if (!file.exists(input_file)) {
  stop(sprintf("Input file does not exist: %s", input_file), call. = F)
}
if (!is.null(deg_file) && !file.exists(deg_file)) {
  stop(sprintf("DEG file does not exist: %s", deg_file), call. = F)
}
ensure_output_directory(output_directory, overwrite = overwrite)

message("Reading long-form sample-by-gene qPCR and bulk RNA-seq data...")
long_data = utils::read.csv(input_file, stringsAsFactors = F, check.names = F)
selected_genes = NULL
# Restrict scope to externally selected DEG identifiers when supplied; this
# filter preserves file order and does not re-test differential expression.
if (!is.null(deg_file)) {
  deg_data = utils::read.csv(deg_file, stringsAsFactors = F, check.names = F)
  if (!deg_gene_col %in% names(deg_data)) {
    stop(sprintf("DEG file has no `%s` column.", deg_gene_col), call. = F)
  }
  selected_genes = unique(as.character(deg_data[[deg_gene_col]]))
  selected_genes = selected_genes[!is.na(selected_genes) & nzchar(selected_genes)]
  if (length(selected_genes) == 0L) {
    stop("The DEG file contains no non-empty gene identifiers.", call. = F)
  }
}

batch_fit = fit_qpcr_bulk_rna_deg_batch(
  long_data,
  sample_id_col = sample_id_col,
  grna_col = grna_col,
  qpcr_col = qpcr_col,
  gene_col = gene_col,
  bulk_rna_col = bulk_rna_col,
  genes = selected_genes,
  expected_replicates = expected_replicates,
  expected_points = expected_points,
  prediction_points = prediction_points
)
message(sprintf(
  "Fitted three models to %d gRNA-level points for each of %d gene(s).",
  expected_points,
  length(batch_fit$genes)
))

# Export the selected input, descriptive percentage-point disagreement, and
# per-gene AIC support separately; none establishes assay or mechanism truth.
selected_long_data = long_data[long_data[[gene_col]] %in% batch_fit$genes, , drop = F]
utils::write.csv(
  selected_long_data,
  file.path(output_directory, "selected_long_input.csv"),
  row.names = F
)
utils::write.csv(
  collect_grna_means(batch_fit),
  file.path(output_directory, "gene_grna_level_means.csv"),
  row.names = F
)
utils::write.csv(
  collect_grna_means(batch_fit),
  file.path(output_directory, "gene_grna_replicate_summaries.csv"),
  row.names = F
)
concordance <- qpcr_bulk_rna_deg_concordance(batch_fit)
concordance_matrix <- data.frame(
  gene = rownames(concordance$matrix), concordance$matrix,
  check.names = FALSE, stringsAsFactors = FALSE
)
utils::write.csv(
  concordance_matrix,
  file.path(output_directory, "batch_concordance_difference_matrix.csv"),
  row.names = FALSE
)
utils::write.csv(
  concordance$long,
  file.path(output_directory, "batch_concordance_difference_long.csv"),
  row.names = FALSE
)
model_overview <- summarize_qpcr_bulk_rna_deg_models(batch_fit)
utils::write.csv(
  model_overview,
  file.path(output_directory, "gene_model_summary.csv"),
  row.names = FALSE
)
for (component in c(
  "model_table", "coefficients", "predictions", "diagnostics", "leave_one_grna_out",
  "linear_anova", "quadratic_anova", "linear_vs_quadratic_anova"
)) {
  output_name = switch(
    component,
    model_table = "gene_model_aic.csv",
    coefficients = "gene_model_coefficients.csv",
    predictions = "gene_model_predictions.csv",
    linear_anova = "gene_linear_anova.csv",
    quadratic_anova = "gene_quadratic_anova.csv",
    linear_vs_quadratic_anova = "gene_linear_vs_quadratic_anova.csv",
    diagnostics = "gene_model_diagnostics.csv",
    leave_one_grna_out = "gene_leave_one_grna_out.csv"
  )
  utils::write.csv(
    collect_qpcr_bulk_rna_deg_results(batch_fit, component),
    file.path(output_directory, output_name),
    row.names = F
  )
}
plot_result = plot_qpcr_bulk_rna_deg_batch(
  batch_fit,
  file = file.path(output_directory, "qpcr_bulk_rna_all_degs.pdf"),
  label_points = label_points
)
diagnostic_plot_result = plot_qpcr_bulk_rna_deg_batch_diagnostics(
  batch_fit,
  file = file.path(output_directory, "qpcr_bulk_rna_all_degs_diagnostics.pdf")
)
heatmap_plot_result = plot_qpcr_bulk_rna_deg_concordance_heatmap(
  batch_fit,
  file = file.path(output_directory, "qpcr_bulk_rna_batch_concordance_heatmap.pdf")
)
model_overview_plot_result = plot_qpcr_bulk_rna_deg_model_overview(
  batch_fit,
  file = file.path(output_directory, "qpcr_bulk_rna_model_overview.pdf")
)
saveRDS(batch_fit, file.path(output_directory, "batch_model_objects.rds"), version = 3)

# Capture selected-input/source hashes, shared plot limits/page counts, options,
# Git state, and session information; Git calls are best-effort outside checkouts.
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
    "run_time_utc", "input", "input_md5", "deg_file", "deg_file_md5", "sample_id_col",
    "grna_col", "qpcr_col", "gene_col", "bulk_rna_col", "genes_fitted", "pdf_pages",
    "expected_replicates", "expected_points", "prediction_points", "label_points", "overwrite",
    "concordance_limits", "model_overview_pages", "drc_version", "runner_source", "runner_source_md5", "model_source", "model_source_md5",
    "cli_helpers_source", "cli_helpers_source_md5", "git_commit", "git_worktree_dirty"
  ),
  value = c(
    format(Sys.time(), tz = "UTC", usetz = T),
    input_file,
    file_md5(input_file),
    if (is.null(deg_file)) "" else deg_file,
    if (is.null(deg_file)) "" else file_md5(deg_file),
    sample_id_col,
    grna_col,
    qpcr_col,
    gene_col,
    bulk_rna_col,
    length(batch_fit$genes),
    plot_result$pages,
    expected_replicates,
    expected_points,
    prediction_points,
    label_points,
    overwrite,
    paste(concordance$limits, collapse = ","),
    model_overview_plot_result$pages,
    batch_fit$fits[[1L]]$settings$drc_version,
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
  quote = F
)
session_report = write_base_session_info(file.path(output_directory, "session_info.txt"))
# Preserve nested inputs, settings, source fingerprints, and an inventory of
# every generated artifact; TSV remains for quick inspection, JSON for reuse.
write_analysis_provenance(
  file.path(output_directory, "analysis_provenance.json"),
  workflow = "qpcr_bulk_rna_deg_batch",
  inputs = list(
    input = list(path = input_file, md5 = file_md5(input_file)),
    deg_file = if (is.null(deg_file)) NULL else list(path = deg_file, md5 = file_md5(deg_file))
  ),
  settings = list(
    sample_id_col = sample_id_col, grna_col = grna_col, qpcr_col = qpcr_col,
    gene_col = gene_col, bulk_rna_col = bulk_rna_col,
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
    input_rows = NROW(long_data), genes = as.list(batch_fit$genes),
    gene_count = length(batch_fit$genes), grna_points_per_gene = expected_points,
    plot_pages = plot_result$pages, concordance_limits = as.list(concordance$limits)
  ),
  output_directory = output_directory
)

cat(sprintf(
  "Completed %d gene(s); each PDF page contains %d gRNA-level points. Results are in %s\n",
  length(batch_fit$genes),
  expected_points,
  output_directory
))
cat(paste(session_report, collapse = "\n"), "\n")
