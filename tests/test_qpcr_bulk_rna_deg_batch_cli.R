#!/usr/bin/env Rscript

script_argument = grep("^--file=", commandArgs(trailingOnly = F), value = T)
script_path = normalizePath(sub("^--file=", "", script_argument), mustWork = T)
repository_directory = normalizePath(file.path(dirname(script_path), ".."), mustWork = T)
rscript = file.path(R.home("bin"), "Rscript")

# Verify the full batch CLI contract: long-form fixture, flattened outputs,
# disagreement sign, model labels, provenance, overwrite safety, and PDF pages.
temporary_directory = tempfile("qpcr-bulk-rna-deg-batch-cli-test-")
dir.create(temporary_directory, recursive = T)
on.exit(unlink(temporary_directory, recursive = T, force = T), add = T)
# Keep any graphics-device regression artifacts out of the repository root.
setwd(temporary_directory)

input_file = file.path(temporary_directory, "example_qpcr_bulk_rna_degs.csv")
output_directory = file.path(temporary_directory, "analysis_results")

run_rscript = function(arguments) {
  output = system2(rscript, shQuote(arguments), stdout = T, stderr = T)
  status = attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop(paste(output, collapse = "\n"), call. = F)
  }
  invisible(output)
}

run_rscript_error = function(arguments, pattern) {
  output = suppressWarnings(system2(rscript, shQuote(arguments), stdout = T, stderr = T))
  status = attr(output, "status")
  stopifnot(
    !is.null(status),
    status != 0L,
    any(grepl(pattern, output, fixed = T))
  )
}

run_rscript(c(
  file.path(repository_directory, "scripts", "generate_qpcr_bulk_rna_deg_example.R"),
  input_file
))
runner = file.path(repository_directory, "analysis", "run_qpcr_bulk_rna_deg_batch.R")
runner_arguments = c(
  runner,
  "--input", input_file,
  "--sample-id-col", "sample_id",
  "--grna-col", "grna",
  "--qpcr-col", "qpcr_percent",
  "--gene-col", "gene",
  "--bulk-rna-col", "bulk_rna_percent",
  "--output-dir", output_directory
)
run_rscript_error(c(runner, "--unknown", "value"), "Unknown command-line option(s): --unknown")
run_rscript_error(c(
  runner,
  "--input", input_file,
  "--sample-id-col", "sample_id",
  "--grna-col", "grna",
  "--qpcr-col", "qpcr_percent",
  "--gene-col", "gene",
  "--bulk-rna-col", "bulk_rna_percent",
  "--expected-points", "11.5",
  "--output-dir", output_directory
), "`--expected-points` must be one positive integer.")
run_rscript(runner_arguments)
# An optional selection file is also an input, even when the main input is elsewhere.
protected_deg = file.path(output_directory, "genes.csv")
write.csv(data.frame(gene="GENE_LINEAR"), protected_deg, row.names=F)
deg_hash = tools::md5sum(protected_deg)
run_rscript_error(c(runner_arguments, "--deg-file", protected_deg, "--overwrite", "true"),
  "must not contain input or source files")
stopifnot(identical(tools::md5sum(protected_deg), deg_hash))
unlink(protected_deg)

expected_outputs = c(
  "selected_long_input.csv",
  "gene_grna_level_means.csv",
  "gene_grna_replicate_summaries.csv",
  "gene_model_diagnostics.csv",
  "gene_leave_one_grna_out.csv",
  "gene_model_summary.csv",
  "batch_concordance_difference_matrix.csv",
  "batch_concordance_difference_long.csv",
  "gene_model_aic.csv",
  "gene_model_coefficients.csv",
  "gene_model_predictions.csv",
  "gene_linear_anova.csv",
  "gene_quadratic_anova.csv",
  "gene_linear_vs_quadratic_anova.csv",
  "qpcr_bulk_rna_all_degs.pdf",
  "qpcr_bulk_rna_all_degs_diagnostics.pdf",
  "qpcr_bulk_rna_batch_concordance_heatmap.pdf",
  "qpcr_bulk_rna_model_overview.pdf",
  "batch_model_objects.rds",
  "run_manifest.tsv",
  "session_info.txt",
  "analysis_provenance.json"
)
output_paths = file.path(output_directory, expected_outputs)
stopifnot(
  all(file.exists(output_paths)),
  all(file.info(output_paths)$size > 0),
  !file.exists(file.path(temporary_directory, "Rplots.pdf"))
)

grna_means = utils::read.csv(file.path(output_directory, "gene_grna_level_means.csv"), stringsAsFactors = F)
model_aic = utils::read.csv(file.path(output_directory, "gene_model_aic.csv"), stringsAsFactors = F)
manifest = utils::read.delim(file.path(output_directory, "run_manifest.tsv"), stringsAsFactors = F)
concordance_matrix = utils::read.csv(file.path(output_directory, "batch_concordance_difference_matrix.csv"), stringsAsFactors = F, check.names = F)
concordance_long = utils::read.csv(file.path(output_directory, "batch_concordance_difference_long.csv"), stringsAsFactors = F)
model_summary = utils::read.csv(file.path(output_directory, "gene_model_summary.csv"), stringsAsFactors = F)
stopifnot(
  NROW(grna_means) == 33L,
  all(table(grna_means$gene) == 11L),
  all(grna_means$n_samples == 4L),
  nrow(concordance_matrix) == 3L, ncol(concordance_matrix) == 12L,
  NROW(concordance_long) == 33L,
  all.equal(concordance_long$difference_percentage_points,
    concordance_long$bulk_rna_percent - concordance_long$qpcr_percent,
    check.attributes = FALSE) == TRUE,
  NROW(model_summary) == 9L,
  all(model_summary$winner_label %in% c("winner", "not_winner", "tied_winner", "failed")),
  NROW(model_aic) == 9L,
  all(table(model_aic$gene) == 3L),
  NROW(utils::read.csv(file.path(output_directory, "gene_model_diagnostics.csv"))) == 99L,
  NROW(utils::read.csv(file.path(output_directory, "gene_leave_one_grna_out.csv"))) == 99L,
  identical(manifest$value[manifest$key == "pdf_pages"], "3")
)
manifest_keys = c("input_md5", "runner_source_md5", "model_source_md5", "cli_helpers_source_md5")
stopifnot(all(manifest_keys %in% manifest$key), all(nzchar(manifest$value[manifest$key %in% manifest_keys])))
manifest_hashes = manifest$value[match(manifest_keys, manifest$key)]
expected_hashes = c(
  unname(tools::md5sum(input_file)),
  unname(tools::md5sum(runner)),
  unname(tools::md5sum(file.path(repository_directory, "R", "qpcr_bulk_rna_models.R"))),
  unname(tools::md5sum(file.path(repository_directory, "R", "cli_helpers.R")))
)
stopifnot(identical(manifest_hashes, as.character(expected_hashes)))
# JSON keeps the nested batch settings, gene list, and artifact inventory that
# a flat TSV cannot represent without lossy comma-separated values.
provenance = paste(readLines(file.path(output_directory, "analysis_provenance.json")), collapse = "")
stopifnot(
  grepl('"schema_version":"1.0"', provenance, fixed = T),
  grepl('"workflow":"qpcr_bulk_rna_deg_batch"', provenance, fixed = T),
  grepl('"artifacts"', provenance, fixed = T),
  grepl('batch_model_objects.rds', provenance, fixed = T)
)
stopifnot(any(grepl("R version", readLines(file.path(output_directory, "session_info.txt")), fixed = T)))
run_rscript_error(runner_arguments, "Output directory is not empty")
run_rscript(c(runner_arguments, "--overwrite", "true"))
stopifnot(!file.exists(file.path(temporary_directory, "Rplots.pdf")))

cat("End-to-end batch qPCR--bulk RNA-seq command-line test passed.\n")
