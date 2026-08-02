#!/usr/bin/env Rscript

script_argument = grep("^--file=", commandArgs(trailingOnly = F), value = T)
script_path = normalizePath(sub("^--file=", "", script_argument), mustWork = T)
repository_directory = normalizePath(file.path(dirname(script_path), ".."), mustWork = T)
rscript = file.path(R.home("bin"), "Rscript")

# Verify the qPCR CLI's strict typing, artifacts, provenance hashes, output
# safety, and absence of unintended graphics files end to end.
temporary_directory = tempfile("qpcr-bulk-rna-cli-test-")
dir.create(temporary_directory, recursive = T)
on.exit(unlink(temporary_directory, recursive = T, force = T), add = T)
# The CLI exercises drc failure handling during leave-one-gRNA-out fitting;
# run child R processes from the disposable directory so any default graphics
# device is cleaned up with the other test outputs rather than written to root.
setwd(temporary_directory)

input_file = file.path(temporary_directory, "example_qpcr_bulk_rna.csv")
output_directory = file.path(temporary_directory, "analysis_results")

run_rscript = function(arguments) {
  output = system2(
    rscript,
    shQuote(arguments),
    stdout = T,
    stderr = T
  )
  status = attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop(paste(output, collapse = "\n"), call. = F)
  }
  invisible(output)
}

run_rscript_error = function(arguments, pattern) {
  output = suppressWarnings(system2(
    rscript,
    shQuote(arguments),
    stdout = T,
    stderr = T
  ))
  status = attr(output, "status")
  stopifnot(
    !is.null(status),
    status != 0L,
    any(grepl(pattern, output, fixed = T))
  )
}

run_rscript(c(
  file.path(repository_directory, "scripts", "generate_qpcr_bulk_rna_example.R"),
  input_file
))
runner = file.path(repository_directory, "analysis", "run_qpcr_bulk_rna_comparison.R")
runner_arguments = c(
  runner,
  "--input", input_file,
  "--grna-col", "grna",
  "--qpcr-col", "qpcr_percent",
  "--bulk-rna-col", "bulk_rna_percent",
  "--output-dir", output_directory
)
run_rscript_error(c(runner, "--unknown", "value"), "Unknown command-line option(s): --unknown")
run_rscript_error(c(
  runner,
  "--input", input_file,
  "--grna-col", "grna",
  "--qpcr-col", "qpcr_percent",
  "--bulk-rna-col", "bulk_rna_percent",
  "--expected-points", "11.5",
  "--output-dir", output_directory
), "`--expected-points` must be one positive integer.")
run_rscript(runner_arguments)

expected_outputs = c(
  "complete_sample_data.csv",
  "grna_level_means.csv",
  "grna_replicate_summaries.csv",
  "model_diagnostics.csv",
  "leave_one_grna_out.csv",
  "model_aic.csv",
  "model_coefficients.csv",
  "model_predictions.csv",
  "linear_anova.csv",
  "quadratic_anova.csv",
  "linear_vs_quadratic_anova.csv",
  "qpcr_bulk_rna_scatterplot.pdf",
  "qpcr_bulk_rna_diagnostics.pdf",
  "model_objects.rds",
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

grna_means = utils::read.csv(
  file.path(output_directory, "grna_level_means.csv"),
  stringsAsFactors = F
)
model_aic = utils::read.csv(
  file.path(output_directory, "model_aic.csv"),
  stringsAsFactors = F
)
stopifnot(
  NROW(grna_means) == 11L,
  all(grna_means$n_samples == 4L),
  identical(model_aic$model, c("linear", "quadratic", "ll4")),
  all(model_aic$converged),
  all(is.finite(model_aic$aic)),
  NROW(utils::read.csv(file.path(output_directory, "model_diagnostics.csv"))) == 33L,
  NROW(utils::read.csv(file.path(output_directory, "leave_one_grna_out.csv"))) == 33L
)
manifest = utils::read.delim(file.path(output_directory, "run_manifest.tsv"), stringsAsFactors = F)
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
# JSON keeps the nested configuration and checksums that a flat TSV cannot.
provenance = paste(readLines(file.path(output_directory, "analysis_provenance.json")), collapse = "")
stopifnot(
  grepl('"schema_version":"1.0"', provenance, fixed = T),
  grepl('"workflow":"qpcr_bulk_rna_comparison"', provenance, fixed = T),
  grepl('"artifacts"', provenance, fixed = T),
  grepl('model_objects.rds', provenance, fixed = T)
)
stopifnot(any(grepl("R version", readLines(file.path(output_directory, "session_info.txt")), fixed = T)))
run_rscript_error(runner_arguments, "Output directory is not empty")
run_rscript(c(runner_arguments, "--overwrite", "true"))
stopifnot(!file.exists(file.path(temporary_directory, "Rplots.pdf")))

cat("End-to-end qPCR--bulk RNA-seq command-line test passed.\n")
