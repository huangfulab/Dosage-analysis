#!/usr/bin/env Rscript

# Verify that the downloaded CSV and delivered TSV contain the same expression
# data, while allowing delimiter and sample-header differences.

options(stringsAsFactors=F)
if (!requireNamespace("devtools", quietly=T)) stop("devtools is required for session information", call.=F)
# Old sessioninfo versions cannot format an empty non-base attached-package
# list. Attach the reporting dependency explicitly before the final report.
suppressPackageStartupMessages(library(devtools))

script_dir = function() {
  file_arg = grep("^--file=", commandArgs(trailingOnly=F), value=T)
  # If not executed via Rscript analysis.R, use the current working directory.
  if (length(file_arg) == 0L) return(normalizePath(getwd(), mustWork=T))

  script_file = sub("^--file=", "", file_arg[[1L]])
  return(dirname(normalizePath(script_file, mustWork=T)))
}

args = commandArgs(trailingOnly=T)
project_dir = script_dir()
# Study-local placeholders; supply the real files here or pass explicit paths.
csv_path = if (length(args) >= 1L) args[[1L]] else file.path(project_dir, "inputs", "K2LLPS-expression-matrix.csv")
tsv_path = if (length(args) >= 2L) args[[2L]] else file.path(project_dir, "inputs", "K2LLPS-expression-matrix.tsv")
output_dir = if (length(args) >= 3L) args[[3L]] else file.path(project_dir, "outputs")
tolerance = if (length(args) >= 4L) suppressWarnings(as.numeric(args[[4L]])) else 1e-14

if (!is.finite(tolerance) || tolerance < 0) stop("Tolerance must be a non-negative number", call.=F)
csv_path = normalizePath(csv_path, mustWork=T)
tsv_path = normalizePath(tsv_path, mustWork=T)
dir.create(output_dir, recursive=T, showWarnings=F)

message("Reading CSV: ", csv_path)
csv_data = read.csv(csv_path, check.names=F)
message("Reading TSV: ", tsv_path)
tsv_data = read.delim(tsv_path, check.names=F)

not_empty = NROW(csv_data) > 0 && NROW(tsv_data) > 0 && NCOL(csv_data) > 0 && NCOL(tsv_data) > 0
same_dimensions = identical(dim(csv_data), dim(tsv_data))
annotation_columns = c("gene_id", "gene_name", "gene_biotype")
annotation_columns_present = all(annotation_columns %in% names(csv_data)) &&
  all(annotation_columns %in% names(tsv_data))

finite_values = F
if (!not_empty || !same_dimensions || !annotation_columns_present) {
  same_annotations = F
  same_expression_exact = F
  same_expression_within_tolerance = F
  differing_expression_cells = NA_integer_
  maximum_absolute_difference = NA_real_
} else {
  same_annotations = identical(
    csv_data[, annotation_columns, drop=F],
    tsv_data[, annotation_columns, drop=F]
  )

  csv_expression_columns = setdiff(names(csv_data), annotation_columns)
  tsv_expression_columns = setdiff(names(tsv_data), annotation_columns)
  same_expression_column_count = length(csv_expression_columns) > 0L &&
    length(csv_expression_columns) == length(tsv_expression_columns)

  if (!same_expression_column_count) {
    same_expression_exact = F
    same_expression_within_tolerance = F
    differing_expression_cells = NA_integer_
    maximum_absolute_difference = NA_real_
  } else {
    # Compare expression columns by position because the CSV uses
    # informative guide/replicate names while the TSV uses K2LLPS_N names.
    csv_expression = as.matrix(csv_data[, csv_expression_columns, drop=F])
    tsv_expression = as.matrix(tsv_data[, tsv_expression_columns, drop=F])
    storage.mode(csv_expression) = "double"
    storage.mode(tsv_expression) = "double"

    finite_values = all(is.finite(csv_expression)) && all(is.finite(tsv_expression))
    differences = abs(csv_expression - tsv_expression)
    differing_expression_cells = sum(differences > tolerance)
    maximum_absolute_difference = max(differences)
    same_expression_exact = identical(unname(csv_expression), unname(tsv_expression))
    same_expression_within_tolerance = finite_values && differing_expression_cells == 0L
  }
}

header_map = data.frame(
  column_index = seq_len(min(NCOL(csv_data), NCOL(tsv_data))),
  csv_header = names(csv_data)[seq_len(min(NCOL(csv_data), NCOL(tsv_data)))],
  tsv_header = names(tsv_data)[seq_len(min(NCOL(csv_data), NCOL(tsv_data)))],
  header_identical = names(csv_data)[seq_len(min(NCOL(csv_data), NCOL(tsv_data)))] ==
    names(tsv_data)[seq_len(min(NCOL(csv_data), NCOL(tsv_data)))],
  stringsAsFactors=F
)
write.csv(header_map, file.path(output_dir, "expression-header-map.csv"), row.names=F, quote=T)

different_header_count = sum(!header_map$header_identical)
practically_equal = same_dimensions && annotation_columns_present &&
  same_annotations && finite_values && same_expression_within_tolerance

report = c(
  "K2LLPS expression-file equivalence check",
  sprintf("CSV: %s", csv_path),
  sprintf("TSV: %s", tsv_path),
  sprintf("Tolerance: %.17g", tolerance),
  sprintf("CSV dimensions: %d rows x %d columns", NROW(csv_data), NCOL(csv_data)),
  sprintf("TSV dimensions: %d rows x %d columns", NROW(tsv_data), NCOL(tsv_data)),
  sprintf("Same dimensions: %s", same_dimensions),
  sprintf("Required annotation columns present: %s", annotation_columns_present),
  sprintf("Annotation rows identical: %s", same_annotations),
  sprintf("Expression values exactly identical: %s", same_expression_exact),
  sprintf("Expression values identical within tolerance: %s", same_expression_within_tolerance),
  sprintf("Expression cells exceeding tolerance: %s", differing_expression_cells),
  sprintf("Maximum absolute expression difference: %s", format(maximum_absolute_difference, scientific=T)),
  sprintf("Differently named headers: %d", different_header_count),
  sprintf("PRACTICALLY_EQUAL: %s", practically_equal),
  "Interpretation: delimiter and header renaming are allowed; annotation and expression content must match."
)

report_path = file.path(output_dir, "expression-file-equivalence.txt")
writeLines(report, report_path)
cat(paste(report, collapse="\n"), "\n")
message("Header map: ", file.path(output_dir, "expression-header-map.csv"))
message("Report: ", report_path)

if (!practically_equal) quit(save="no", status=1L)

devtools::session_info(pkgs="attached")
