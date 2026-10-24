#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument), mustWork = TRUE)
repository_directory <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
rscript <- file.path(R.home("bin"), "Rscript")
source(file.path(repository_directory, "R", "cli_helpers.R"))
# Provenance must remain valid JSON for quoted names, Windows paths, and every
# representable ASCII control character, not just ordinary fixture paths.
stopifnot(
  identical(.json_encode('a"b'), '"a\\"b"'),
  identical(.json_encode(paste0("C:", intToUtf8(92), "data")), '"C:\\\\data"'),
  identical(.json_escape(intToUtf8(1:31)), paste(sprintf("\\u%04x", 1:31), collapse="")),
  identical(.json_encode("한글"), '"한글"'),
  identical(as.numeric(.json_encode(pi)), pi)
)

# Exercise the public end-to-end contract: generated input, strict CLI errors,
# outputs/provenance, and refusal versus explicit destructive overwrite.
temporary_directory <- tempfile("dosage-cli-test-")
dir.create(temporary_directory, recursive = TRUE)
on.exit(unlink(temporary_directory, recursive = TRUE, force = TRUE), add = TRUE)

input_file <- file.path(temporary_directory, if (.Platform$OS.type == "unix")
  "example\t\"data\"\n.csv" else "example data.csv")
output_directory <- file.path(temporary_directory, "analysis results")

run_rscript <- function(arguments) {
  output <- system2(
    rscript,
    shQuote(arguments),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop(paste(output, collapse = "\n"), call. = FALSE)
  }
  invisible(output)
}

run_rscript(c(
  file.path(repository_directory, "scripts", "generate_example_data.R"),
  input_file
))
runner <- file.path(repository_directory, "analysis", "run_dosage_analysis.R")
runner_arguments <- c(
  runner,
  "--input", input_file,
  "--dose-col", "dose",
  "--response-cols", "response_hill,response_quadratic",
  "--group-col", "group",
  "--control-label", "Controls",
  "--normalize-dose", "true",
  "--normalize-response", "true",
  "--output-dir", output_directory
)
run_rscript_error <- function(arguments, pattern) {
  output <- suppressWarnings(system2(rscript, shQuote(arguments), stdout = TRUE, stderr = TRUE))
  status <- attr(output, "status")
  stopifnot(!is.null(status), status != 0L, any(grepl(pattern, output, fixed = TRUE)))
}
# Distinct endpoint names that sanitize to the same filename must fail before
# creating an output directory or silently overwriting per-response artifacts.
collision_input = file.path(temporary_directory, "collisions.csv")
utils::write.csv(data.frame(dose=1:6, "RNA/A"=2:7, "RNA A"=3:8, check.names=F),
  collision_input, row.names=F)
collision_output = file.path(temporary_directory, "collision-output")
run_rscript_error(c(runner, "--input", collision_input, "--dose-col", "dose",
  "--response-cols", "RNA/A,RNA A", "--output-dir", collision_output), "unique output filename stems")
stopifnot(!dir.exists(collision_output))
protected_input_hash = file_md5(input_file)
protected_arguments = runner_arguments
protected_arguments[which(protected_arguments == "--output-dir") + 1L] = temporary_directory
run_rscript_error(c(protected_arguments, "--overwrite", "true"), "must not contain input or source files")
stopifnot(identical(file_md5(input_file), protected_input_hash))
# A missing intermediate component must not bypass the containment check.
protected_arguments[which(protected_arguments == "--output-dir") + 1L] =
  file.path(temporary_directory, "new-component", "..")
run_rscript_error(c(protected_arguments, "--overwrite", "true"), "must not contain input or source files")
stopifnot(identical(file_md5(input_file), protected_input_hash))
if (.Platform$OS.type == "unix") {
  link_output = file.path(temporary_directory, "symlink-output")
  dir.create(link_output)
  linked_input = file.path(link_output, "input.csv")
  stopifnot(file.symlink(input_file, linked_input))
  link_error = tryCatch(ensure_output_directory(link_output, overwrite=T,
    protected_paths=linked_input), error=identity)
  stopifnot(inherits(link_error, "error"), file.exists(linked_input),
    identical(file_md5(input_file), protected_input_hash))
}
run_rscript_error(c(runner, "--unknown", "value"), "Unknown command-line option(s): --unknown")
run_rscript(runner_arguments)

expected_outputs <- c(
  "prepared_data.csv",
  "normalization_factors.csv",
  "model_comparison.csv",
  "model_metrics.csv",
  "model_coefficients.csv",
  "model_predictions.csv",
  "polynomial_nested_tests.csv",
  "model_diagnostics.csv",
  "model_diagnostics_response_hill.csv",
  "model_diagnostics_response_quadratic.csv",
  "model_objects.rds",
  "model_comparison_response_hill.pdf",
  "model_diagnostics_response_hill.pdf",
  "model_diagnostics_response_quadratic.pdf",
  "model_comparison_response_quadratic.pdf",
  "run_manifest.tsv",
  "session_info.txt",
  "analysis_provenance.json"
)
output_paths <- file.path(output_directory, expected_outputs)
stopifnot(
  all(file.exists(output_paths)),
  all(file.info(output_paths)$size > 0)
)

comparison <- utils::read.csv(
  file.path(output_directory, "model_comparison.csv"),
  stringsAsFactors = FALSE
)
diagnostics <- utils::read.csv(file.path(output_directory, "model_diagnostics.csv"), stringsAsFactors = FALSE)
stopifnot(
  setequal(unique(comparison$response), c("response_hill", "response_quadratic")),
  all(comparison$converged),
  all(is.finite(comparison$aic)),
  all(comparison$aic_rank >= 1L),
  !any(comparison$model %in% c("polynomial_degree_3", "polynomial_degree_4")),
  NROW(diagnostics) == 320L,
  all.equal(diagnostics$residual, diagnostics$observed - diagnostics$fitted,
    check.attributes = FALSE) == TRUE,
  all(is.finite(diagnostics$leverage[diagnostics$family == "polynomial"])),
  all(is.na(diagnostics$leverage[diagnostics$family != "polynomial"]))
)
manifest <- utils::read.delim(file.path(output_directory, "run_manifest.tsv"), stringsAsFactors = FALSE)
manifest_keys <- c("input_md5", "runner_source_md5", "model_source_md5", "cli_helpers_source_md5")
stopifnot(identical(manifest$value[manifest$key == "input"], input_file))
stopifnot(all(manifest_keys %in% manifest$key), all(nzchar(manifest$value[manifest$key %in% manifest_keys])))
manifest_hashes <- manifest$value[match(manifest_keys, manifest$key)]
expected_hashes <- c(
  unname(tools::md5sum(input_file)),
  unname(tools::md5sum(runner)),
  unname(tools::md5sum(file.path(repository_directory, "R", "dosage_models.R"))),
  unname(tools::md5sum(file.path(repository_directory, "R", "cli_helpers.R")))
)
stopifnot(identical(manifest_hashes, as.character(expected_hashes)))
# JSON keeps the nested configuration and checksums that a flat TSV cannot.
provenance <- paste(readLines(file.path(output_directory, "analysis_provenance.json")), collapse = "")
stopifnot(
  grepl('"schema_version":"1.0"', provenance, fixed = TRUE),
  grepl('"workflow":"dosage_analysis"', provenance, fixed = TRUE),
  grepl('"artifacts"', provenance, fixed = TRUE),
  grepl('model_objects.rds', provenance, fixed = TRUE)
)
stopifnot(any(grepl("R version", readLines(file.path(output_directory, "session_info.txt")), fixed = TRUE)))
run_rscript_error(runner_arguments, "Output directory is not empty")
stale_output <- file.path(output_directory, "stale-output.txt")
writeLines("stale", stale_output)
stopifnot(file.exists(stale_output))
run_rscript(c(runner_arguments, "--overwrite", "true"))
stopifnot(!file.exists(stale_output))

cat("End-to-end command-line test passed.\n")
