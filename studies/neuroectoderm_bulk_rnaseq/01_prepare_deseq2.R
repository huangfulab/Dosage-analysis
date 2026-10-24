#!/usr/bin/env Rscript

# K2LLPS 3'-seq: blocked DESeq2 comparison of OCT4-targeting guides against
# pooled nontargeting and safe-harbor controls.

options(stringsAsFactors=F)

script_dir = function() {
  file_arg = grep("^--file=", commandArgs(trailingOnly=F), value=T)
  # If not executed via Rscript analysis.R, use the current working directory.
  if (length(file_arg) == 0L) return(normalizePath(getwd(), mustWork=T))

  script_file = sub("^--file=", "", file_arg[[1L]])
  return(dirname(normalizePath(script_file, mustWork=T)))
}
stop_with = function(...) stop(sprintf(...), call.=F)
md5_value = function(path) unname(tools::md5sum(path)[[1L]])
max_or_na = function(values) {
  values = values[is.finite(values)]
  if (length(values) == 0L) return(NA_real_)
  max(values)
}
quantile_or_na = function(values, probability) {
  values = values[is.finite(values)]
  if (length(values) == 0L) return(NA_real_)
  unname(quantile(values, probs=probability, names=F))
}

args = commandArgs(trailingOnly=T)
project_dir = script_dir()
# Study-local input placeholders are not distributed with this code.
input_path = if (length(args) >= 1L) args[[1L]] else file.path(project_dir, "inputs", "K2LLPS-expression-matrix.tsv")
metadata_path = if (length(args) >= 2L) args[[2L]] else file.path(project_dir, "sample_metadata_draft.csv")
output_dir = if (length(args) >= 3L) args[[3L]] else file.path(project_dir, "outputs")
design_workbook_path = if (length(args) >= 4L) args[[4L]] else file.path(project_dir, "inputs", "260320_OCT4_qPCR.xlsx")

# "round" is dependency-light and compatible with standard DESeq2.
# "glmGamPoi" retains fractional featureCounts weights but requires glmGamPoi.
count_strategy = Sys.getenv("K2LLPS_COUNT_STRATEGY", unset="round")
valid_strategies = c("round", "glmGamPoi")
if (!count_strategy %in% valid_strategies) {
  stop_with("K2LLPS_COUNT_STRATEGY must be one of: %s", paste(valid_strategies, collapse=", "))
}

min_total_count = suppressWarnings(as.numeric(Sys.getenv("K2LLPS_MIN_TOTAL_COUNT", unset="10")))
if (!is.finite(min_total_count) || min_total_count < 0) {
  stop_with("K2LLPS_MIN_TOTAL_COUNT must be a non-negative number")
}

alpha = suppressWarnings(as.numeric(Sys.getenv("K2LLPS_ALPHA", unset="0.05")))
if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) {
  stop_with("K2LLPS_ALPHA must be a number strictly between 0 and 1")
}

# Hard effect-size threshold for the DEG heatmap/violin gene set. Raw (unshrunken)
# LFC is used so the cutoff is not deflated by normal-prior shrinkage.
deg_lfc_threshold = suppressWarnings(as.numeric(Sys.getenv("K2LLPS_DEG_LFC", unset="1")))
if (!is.finite(deg_lfc_threshold) || deg_lfc_threshold < 0) {
  stop_with("K2LLPS_DEG_LFC must be a non-negative number")
}

if (!requireNamespace("DESeq2", quietly=T)) {
  stop_with("DESeq2 is required but is not installed")
}
if (!requireNamespace("openxlsx", quietly=T)) {
  stop_with("openxlsx is required to read the design workbook but is not installed")
}
if (!requireNamespace("devtools", quietly=T)) {
  stop_with("devtools is required to capture the session information but is not installed")
}
if (!requireNamespace("ggplot2", quietly=T)) {
  stop_with("ggplot2 is required for the PCA eigenvector arrow plot but is not installed")
}
if (!requireNamespace("ggrepel", quietly=T)) {
  stop_with("ggrepel is required for non-overlapping gene labels on the PCA eigenvector plot but is not installed")
}
if (!requireNamespace("pheatmap", quietly=T)) {
  stop_with("pheatmap is required for the DEG expression heatmap but is not installed")
}
if (count_strategy == "glmGamPoi" && !requireNamespace("glmGamPoi", quietly=T)) {
  stop_with("glmGamPoi mode requested, but the glmGamPoi package is not installed")
}

suppressPackageStartupMessages(library(DESeq2))

dir.create(output_dir, recursive=T, showWarnings=F)
input_path = normalizePath(input_path, mustWork=T)
design_workbook_path = normalizePath(design_workbook_path, mustWork=T)

message("Reading expression matrix: ", input_path)
expression = read.delim(input_path, check.names=F)

required_annotation = c("gene_id", "gene_name", "gene_biotype")
missing_annotation = setdiff(required_annotation, names(expression))
if (length(missing_annotation) > 0L) {
  stop_with("Missing annotation column(s): %s", paste(missing_annotation, collapse=", "))
}
if (anyNA(expression$gene_id) || any(expression$gene_id == "") || anyDuplicated(expression$gene_id)) {
  stop_with("gene_id must be non-missing, non-empty, and unique")
}

count_columns = grep("_count$", names(expression), value=T)
cpm_columns = grep("_cpm$", names(expression), value=T)
if (length(count_columns) == 0L) stop_with("No *_count columns were found")

sample_ids = sub("_count$", "", count_columns)
cpm_sample_ids = sub("_cpm$", "", cpm_columns)
if (!setequal(sample_ids, cpm_sample_ids)) {
  stop_with("The *_count and *_cpm sample sets do not match")
}
if (anyDuplicated(sample_ids)) stop_with("Sample identifiers derived from *_count columns are not unique")
if (anyDuplicated(cpm_sample_ids)) stop_with("Sample identifiers derived from *_cpm columns are not unique")

count_matrix = as.matrix(expression[, count_columns, drop=F])
storage.mode(count_matrix) = "double"
rownames(count_matrix) = expression$gene_id
colnames(count_matrix) = sample_ids

if (anyNA(count_matrix) || any(!is.finite(count_matrix)) || any(count_matrix < 0)) {
  stop_with("Count matrix contains missing, non-finite, or negative values")
}

fractional_mask = abs(count_matrix - round(count_matrix)) > 1e-8
fractional_cells = sum(fractional_mask)
message(sprintf(
  "Loaded %s genes x %s samples; %s count cells (%.2f%%) are fractional",
  format(NROW(count_matrix), big.mark=","),
  NCOL(count_matrix),
  format(fractional_cells, big.mark=","),
  100 * fractional_cells / length(count_matrix)
))

# The supplied CPM columns are not used. They are retained only in the source
# file and are intentionally excluded from the DESeq2 object.
if (count_strategy == "round") {
  analysis_counts = round(count_matrix)
  skip_integer_mode = F
  dispersion_fit_type = "parametric"
} else {
  analysis_counts = count_matrix
  skip_integer_mode = T
  dispersion_fit_type = "glmGamPoi"
}

if (!file.exists(metadata_path)) {
  # The TSV headers contain only K2LLPS_* identifiers, so guide and replicate
  # cannot be inferred from the TSV alone. Recover the biological sample names
  # from the CSV-to-TSV header map produced by 00_verify_expression_files.R,
  # but only after confirming that the verifier passed.
  verification_dir = output_dir
  equivalence_report_path = file.path(verification_dir, "expression-file-equivalence.txt")
  header_map_path = file.path(verification_dir, "expression-header-map.csv")

  if (!file.exists(equivalence_report_path) || !file.exists(header_map_path)) {
    stop_with("Run 00_verify_expression_files.R before generating draft metadata")
  }

  equivalence_report = readLines(equivalence_report_path, warn=F)
  if (!any(trimws(equivalence_report) == "PRACTICALLY_EQUAL: TRUE")) {
    stop_with("Cannot generate metadata because the expression-file equivalence check did not pass")
  }

  verified_csv_line = grep("^CSV: ", equivalence_report, value=T)
  verified_tsv_line = grep("^TSV: ", equivalence_report, value=T)
  if (length(verified_csv_line) != 1L || length(verified_tsv_line) != 1L) {
    stop_with("Equivalence report must identify exactly one CSV and one TSV input")
  }
  verified_csv_path = normalizePath(sub("^CSV: ", "", verified_csv_line), mustWork=T)
  verified_tsv_path = normalizePath(sub("^TSV: ", "", verified_tsv_line), mustWork=T)
  if (!identical(verified_tsv_path, input_path)) {
    stop_with("Equivalence report was generated for a different TSV input")
  }
  newest_input_time = max(file.info(c(verified_csv_path, verified_tsv_path))$mtime)
  oldest_verification_time = min(file.info(c(equivalence_report_path, header_map_path))$mtime)
  if (oldest_verification_time < newest_input_time) {
    stop_with("Expression input is newer than the verification outputs; rerun 00_verify_expression_files.R")
  }

  header_map = read.csv(header_map_path, check.names=F)
  required_header_map = c("csv_header", "tsv_header")
  missing_header_map = setdiff(required_header_map, names(header_map))
  if (length(missing_header_map) > 0L) {
    stop_with("Header map must contain column(s): %s", paste(missing_header_map, collapse=", "))
  }
  if (anyNA(header_map$csv_header) || anyNA(header_map$tsv_header)) {
    stop_with("Header map contains missing CSV or TSV headers")
  }

  count_header_map = header_map[
    grepl("_count$", header_map$csv_header) & grepl("_count$", header_map$tsv_header),
    required_header_map,
    drop=F
  ]
  if (anyDuplicated(count_header_map$csv_header) || anyDuplicated(count_header_map$tsv_header)) {
    stop_with("Count headers in the header map must be unique")
  }
  if (!setequal(count_columns, count_header_map$tsv_header)) {
    stop_with("TSV count columns do not match the verified header map")
  }

  count_header_map = count_header_map[match(count_columns, count_header_map$tsv_header), , drop=F]
  biological_sample_ids = sub("_count$", "", count_header_map$csv_header)
  if (anyDuplicated(biological_sample_ids) || !all(grepl("_rep_[^_]+$", biological_sample_ids))) {
    stop_with("CSV count headers do not contain unique '<guide>_rep_<replicate>' sample names")
  }

  inferred_metadata = data.frame(
    count_column = count_header_map$tsv_header,
    sample_id = biological_sample_ids,
    guide = sub("_rep_[^_]+$", "", biological_sample_ids),
    replicate = sub("^.*_rep_", "", biological_sample_ids),
    stringsAsFactors=F
  )

  dir.create(dirname(metadata_path), recursive=T, showWarnings=F)
  write.csv(inferred_metadata, metadata_path, row.names=F, quote=F)
  message("Created draft metadata from the verified CSV-to-TSV header map: ", metadata_path)
}
metadata_path = normalizePath(metadata_path, mustWork=T)

metadata = read.csv(metadata_path, check.names=F)
required_metadata = c("count_column", "sample_id", "guide", "replicate")
missing_metadata = setdiff(required_metadata, names(metadata))
if (length(missing_metadata) > 0L) {
  stop_with("Metadata must contain column(s): %s", paste(missing_metadata, collapse=", "))
}
if (anyNA(metadata$count_column) || any(metadata$count_column == "") || anyDuplicated(metadata$count_column)) {
  stop_with("metadata$count_column must be non-missing, non-empty, and unique")
}
if (anyNA(metadata$sample_id) || any(metadata$sample_id == "") || anyDuplicated(metadata$sample_id)) {
  stop_with("metadata$sample_id must be non-missing, non-empty, and unique")
}
if (!setequal(count_columns, metadata$count_column)) {
  missing_from_metadata = setdiff(count_columns, metadata$count_column)
  extra_in_metadata = setdiff(metadata$count_column, count_columns)
  stop_with(
    "Metadata/count-column mismatch. Missing: [%s]. Extra: [%s]",
    paste(missing_from_metadata, collapse=", "),
    paste(extra_in_metadata, collapse=", ")
  )
}

metadata = metadata[match(count_columns, metadata$count_column), , drop=F]
if (!identical(count_columns, metadata$count_column)) {
  stop_with("Internal error while aligning metadata to expression columns")
}
metadata$guide = as.character(metadata$guide)
metadata$replicate = as.character(metadata$replicate)
expected_sample_ids = paste0(metadata$guide, "_rep_", metadata$replicate)
if (!identical(metadata$sample_id, expected_sample_ids)) {
  stop_with("Every sample_id must equal '<guide>_rep_<replicate>' after count-column alignment")
}

# Link the biological guide labels and post-treatment qPCR values from the design
# workbook by explicit guide and replicate keys.
target_guides = c("guide01", "guide02", "guide03", "guide04", "guide05", "guide06", "guide14", "guide15")
nontargeting_guides = c("guide18", "guide19")
safe_harbor_guides = "guide20"
control_guides = c(nontargeting_guides, safe_harbor_guides)
expected_guides = c(target_guides, control_guides)
expected_replicates = as.character(2:5)

message("Reading design workbook: ", design_workbook_path)
design_sheet_name = "Sheet1"
design_sheet = openxlsx::read.xlsx(
  design_workbook_path,
  sheet=design_sheet_name,
  check.names=F,
  na.strings=c("NA", "")
)
required_design_columns = c(
  "HZ_gRNA_name", "JY_gRNA_name",
  "rep2_OCT4", "rep3_OCT4", "rep4_OCT4", "rep5_OCT4"
)
missing_design_columns = setdiff(required_design_columns, names(design_sheet))
if (length(missing_design_columns) > 0L) {
  stop_with("Design workbook is missing column(s): %s", paste(missing_design_columns, collapse=", "))
}
qpcr_columns = grep("^rep[0-9]+_OCT4$", names(design_sheet), value=T)
expected_qpcr_columns = paste0("rep", expected_replicates, "_OCT4")
if (!setequal(qpcr_columns, expected_qpcr_columns)) {
  stop_with(
    "Design workbook qPCR columns must be exactly [%s]; found [%s]",
    paste(expected_qpcr_columns, collapse=", "),
    paste(qpcr_columns, collapse=", ")
  )
}
qpcr_columns = expected_qpcr_columns

hz_guide_number = suppressWarnings(as.numeric(design_sheet$HZ_gRNA_name))
if (
  anyNA(hz_guide_number) || any(!is.finite(hz_guide_number)) ||
  any(abs(hz_guide_number - round(hz_guide_number)) > 1e-8) ||
  anyDuplicated(hz_guide_number)
) {
  stop_with("HZ_gRNA_name must contain unique, finite whole-number guide IDs")
}
jy_guide_label = trimws(as.character(design_sheet$JY_gRNA_name))
if (anyNA(jy_guide_label) || any(jy_guide_label == "") || anyDuplicated(jy_guide_label)) {
  stop_with("JY_gRNA_name must contain unique, non-missing, non-empty labels")
}
design_sheet$guide = sprintf("guide%02d", as.integer(hz_guide_number))
if (!setequal(design_sheet$guide, expected_guides)) {
  stop_with(
    "Workbook guide mismatch. Expected [%s]; found [%s]",
    paste(expected_guides, collapse=", "),
    paste(design_sheet$guide, collapse=", ")
  )
}

qpcr_long_parts = vector("list", length(qpcr_columns))
for (column_index in seq_along(qpcr_columns)) {
  qpcr_column = qpcr_columns[[column_index]]
  replicate_value = sub("^rep([0-9]+)_OCT4$", "\\1", qpcr_column)
  raw_values = design_sheet[[qpcr_column]]
  numeric_values = suppressWarnings(as.numeric(raw_values))
  invalid_numeric = !is.na(raw_values) & trimws(as.character(raw_values)) != "" & is.na(numeric_values)
  if (any(invalid_numeric)) {
    stop_with("Column %s contains nonnumeric qPCR value(s)", qpcr_column)
  }

  qpcr_long_parts[[column_index]] = data.frame(
    guide = design_sheet$guide,
    replicate = replicate_value,
    hz_guide_number = as.integer(hz_guide_number),
    jy_guide_label = jy_guide_label,
    qpcr_oct4_post_treatment = numeric_values,
    stringsAsFactors=F,
    check.names=F
  )
}
qpcr_long = do.call(rbind, qpcr_long_parts)
qpcr_long$key = paste(qpcr_long$guide, qpcr_long$replicate, sep="::")
metadata$key = paste(metadata$guide, metadata$replicate, sep="::")

if (NROW(metadata) != 44L || anyDuplicated(metadata$key)) {
  stop_with("Metadata must contain exactly 44 unique guide/replicate combinations")
}
if (NROW(qpcr_long) != 44L || anyDuplicated(qpcr_long$key)) {
  stop_with("Design workbook must yield exactly 44 unique guide/replicate combinations")
}
if (!setequal(unique(metadata$guide), expected_guides)) {
  stop_with("Metadata guide set does not match the expected 8 targeting and 3 control guides")
}
if (!setequal(unique(metadata$replicate), expected_replicates)) {
  stop_with("Metadata replicate set must be exactly 2, 3, 4, and 5")
}
guide_replicate_table = table(metadata$guide, metadata$replicate)
if (
  !identical(sort(rownames(guide_replicate_table)), sort(expected_guides)) ||
  !identical(sort(colnames(guide_replicate_table)), sort(expected_replicates)) ||
  any(guide_replicate_table != 1L)
) {
  stop_with("Every guide must occur exactly once in each of replicates 2, 3, 4, and 5")
}

qpcr_match = match(metadata$key, qpcr_long$key)
if (anyNA(qpcr_match) || anyDuplicated(qpcr_match) || !setequal(metadata$key, qpcr_long$key)) {
  stop_with("Metadata and workbook guide/replicate keys do not form a one-to-one match")
}
metadata$hz_guide_number = qpcr_long$hz_guide_number[qpcr_match]
metadata$jy_guide_label = qpcr_long$jy_guide_label[qpcr_match]
metadata$qpcr_oct4_post_treatment = qpcr_long$qpcr_oct4_post_treatment[qpcr_match]
metadata$qpcr_missing = is.na(metadata$qpcr_oct4_post_treatment)
metadata$guide_class = ifelse(
  metadata$guide %in% target_guides,
  "OCT4_targeting",
  ifelse(metadata$guide %in% nontargeting_guides, "nontargeting_control", "safe_harbor_control")
)
metadata$comparison_group = ifelse(metadata$guide %in% target_guides, "OCT4_targeting", "pooled_control")

expected_missing_qpcr_sample = "guide04_rep_5"
missing_qpcr_samples = metadata$sample_id[metadata$qpcr_missing]
if (!identical(missing_qpcr_samples, expected_missing_qpcr_sample)) {
  stop_with(
    "Expected the only missing qPCR value at %s; found [%s]",
    expected_missing_qpcr_sample,
    paste(missing_qpcr_samples, collapse=", ")
  )
}

metadata$guide = factor(metadata$guide, levels=expected_guides)
metadata$replicate = factor(metadata$replicate, levels=expected_replicates)
metadata$guide_class = factor(
  metadata$guide_class,
  levels=c("OCT4_targeting", "nontargeting_control", "safe_harbor_control")
)
metadata$comparison_group = factor(
  metadata$comparison_group,
  levels=c("pooled_control", "OCT4_targeting")
)
metadata$key = NULL
rownames(metadata) = metadata$sample_id

if (sum(metadata$comparison_group == "OCT4_targeting") != 32L || sum(metadata$comparison_group == "pooled_control") != 12L) {
  stop_with("Expected 32 OCT4-targeting libraries and 12 pooled-control libraries")
}
if (length(unique(as.character(metadata$guide[metadata$comparison_group == "OCT4_targeting"]))) != 8L ||
    length(unique(as.character(metadata$guide[metadata$comparison_group == "pooled_control"]))) != 3L) {
  stop_with("Expected 8 OCT4-targeting guides and 3 pooled-control guides")
}

write.csv(
  data.frame(metadata, row.names=NULL, check.names=F),
  file.path(output_dir, "analysis-sample-metadata.csv"),
  row.names=F,
  quote=F,
  na="NA"
)

# The TSV column name is the stable input key; sample_id is the biological name
# used for the DESeq2 object's columns after the explicit mapping above.
colnames(count_matrix) = metadata$sample_id
colnames(analysis_counts) = metadata$sample_id

design_formula = ~ replicate + comparison_group
if (any(grepl("qpcr", all.vars(design_formula), ignore.case=T))) {
  stop_with("Post-treatment qPCR measurements must not be included in the DESeq2 design")
}
model_matrix = model.matrix(design_formula, data=metadata)
model_rank = qr(model_matrix)$rank
if (model_rank != NCOL(model_matrix)) {
  stop_with("The DESeq2 model matrix is not full rank")
}
residual_degrees_freedom = NROW(model_matrix) - model_rank
if (residual_degrees_freedom <= 0L) {
  stop_with("The DESeq2 model has no residual degrees of freedom")
}

if (skip_integer_mode) {
  # skipIntegerMode is an argument to DESeqDataSet(), not to the convenience
  # DESeqDataSetFromMatrix() constructor.
  se = SummarizedExperiment::SummarizedExperiment(
    assays=list(counts=analysis_counts),
    colData=S4Vectors::DataFrame(metadata)
  )
  DESeq2Obj = DESeqDataSet(se, design=design_formula, skipIntegerMode=T)
} else {
  DESeq2Obj = DESeqDataSetFromMatrix(
    countData=analysis_counts,
    colData=metadata,
    design=design_formula
  )
}

mcols(DESeq2Obj)$gene_name = expression$gene_name
mcols(DESeq2Obj)$gene_biotype = expression$gene_biotype

# This intentionally permissive filter removes only genes with almost no data;
# DESeq2's independent filtering is also applied when results are extracted.
genes_before_filter = NROW(DESeq2Obj)
keep = rowSums(counts(DESeq2Obj)) >= min_total_count
DESeq2Obj = DESeq2Obj[keep, ]
genes_after_filter = NROW(DESeq2Obj)
message(sprintf(
  "Prefilter retained %s / %s genes (total analysis count >= %s)",
  format(genes_after_filter, big.mark=","),
  format(genes_before_filter, big.mark=","),
  format(min_total_count, scientific=F)
))

# Estimate size factors and dispersions, then fit the Wald model.
DESeq2Obj = estimateSizeFactors(DESeq2Obj)
DESeq2Obj = estimateDispersions(DESeq2Obj, fitType=dispersion_fit_type)

DESeq2Obj = nbinomWaldTest(DESeq2Obj, betaPrior=F)
coefficient_name = "comparison_group_OCT4_targeting_vs_pooled_control"
if (!coefficient_name %in% resultsNames(DESeq2Obj)) {
  stop_with(
    "Expected coefficient '%s' was not found. Available coefficients: %s",
    coefficient_name,
    paste(resultsNames(DESeq2Obj), collapse=", ")
  )
}
diffExprLst = results(
  DESeq2Obj,
  name=coefficient_name,
  alpha=alpha,
  independentFiltering=T,
  cooksCutoff=T
)
diffExprLst.shrunk = lfcShrink(
  DESeq2Obj,
  coef=coefficient_name,
  res=diffExprLst,
  type="normal"
)
saveRDS(DESeq2Obj, file.path(output_dir, "DESeq2Obj-after-wald-test.rds"))

normalized_counts = counts(DESeq2Obj, normalized=T)
targeting_samples = colData(DESeq2Obj)$comparison_group == "OCT4_targeting"
control_samples = colData(DESeq2Obj)$comparison_group == "pooled_control"
mean_targeting = rowMeans(normalized_counts[, targeting_samples, drop=F])
mean_control = rowMeans(normalized_counts[, control_samples, drop=F])

diffExprLst.df = as.data.frame(diffExprLst)
diffExprLst.shrunk.df = as.data.frame(diffExprLst.shrunk)
result_gene_ids = rownames(diffExprLst.df)
result_table = data.frame(
  gene_id = result_gene_ids,
  gene_name = mcols(DESeq2Obj)[result_gene_ids, "gene_name"],
  gene_biotype = mcols(DESeq2Obj)[result_gene_ids, "gene_biotype"],
  baseMean = diffExprLst.df$baseMean,
  mean_normalized_count_oct4_targeting = mean_targeting[result_gene_ids],
  mean_normalized_count_pooled_control = mean_control[result_gene_ids],
  log2FoldChange_raw = diffExprLst.df$log2FoldChange,
  lfcSE_raw = diffExprLst.df$lfcSE,
  log2FoldChange_shrunken = diffExprLst.shrunk.df$log2FoldChange,
  lfcSE_shrunken = diffExprLst.shrunk.df$lfcSE,
  stat = diffExprLst.df$stat,
  pvalue = diffExprLst.df$pvalue,
  padj = diffExprLst.df$padj,
  stringsAsFactors=F,
  check.names=F
)
result_table$is_tested = !is.na(result_table$pvalue)
result_table$has_adjusted_pvalue = !is.na(result_table$padj)
result_table$independent_filter_flag = !is.na(result_table$pvalue) & is.na(result_table$padj)
result_table$cooks_or_other_filter_flag = is.na(result_table$pvalue) & result_table$baseMean > 0
result_table$significant_fdr = !is.na(result_table$padj) & result_table$padj < alpha
result_table$significant_up = result_table$significant_fdr & !is.na(result_table$log2FoldChange_shrunken) & result_table$log2FoldChange_shrunken > 0
result_table$significant_down = result_table$significant_fdr & !is.na(result_table$log2FoldChange_shrunken) & result_table$log2FoldChange_shrunken < 0
result_table$direction = ifelse(
  result_table$significant_up,
  "higher_in_OCT4_targeting",
  ifelse(result_table$significant_down, "lower_in_OCT4_targeting", "not_significant")
)
result_table$result_status = ifelse(
  result_table$baseMean == 0,
  "all_zero",
  ifelse(
    result_table$cooks_or_other_filter_flag,
    "not_tested_cooks_or_other",
    ifelse(
      result_table$independent_filter_flag,
      "independent_filter",
      ifelse(result_table$significant_fdr, "significant_fdr", "tested_not_significant")
    )
  )
)

padj_sort = ifelse(is.na(result_table$padj), Inf, result_table$padj)
effect_sort = ifelse(is.na(result_table$log2FoldChange_shrunken), -Inf, abs(result_table$log2FoldChange_shrunken))
result_order = order(!result_table$significant_fdr, padj_sort, -effect_sort, result_table$gene_id)
result_table = result_table[result_order, , drop=F]
rownames(result_table) = NULL
significant_table = result_table[result_table$significant_fdr, , drop=F]

complete_results_path = file.path(output_dir, "de-results-oct4-targeting-vs-pooled-controls-complete.csv")
significant_results_path = file.path(output_dir, "de-results-oct4-targeting-vs-pooled-controls-significant.csv")
write.csv(result_table, complete_results_path, row.names=F, quote=F, na="NA")
write.csv(significant_table, significant_results_path, row.names=F, quote=F, na="NA")

cooks_matrix = assays(DESeq2Obj)[["cooks"]]
cooks_max = apply(cooks_matrix, 2L, max_or_na)
cooks_q99 = apply(cooks_matrix, 2L, quantile_or_na, probability=0.99)
analysis_library_sums = colSums(analysis_counts)
median_analysis_library_sum = median(analysis_library_sums)

sample_qc = data.frame(
  count_column = metadata$count_column,
  sample_id = metadata$sample_id,
  guide = as.character(metadata$guide),
  hz_guide_number = metadata$hz_guide_number,
  jy_guide_label = metadata$jy_guide_label,
  replicate = as.character(metadata$replicate),
  guide_class = as.character(metadata$guide_class),
  comparison_group = as.character(metadata$comparison_group),
  qpcr_oct4_post_treatment = metadata$qpcr_oct4_post_treatment,
  qpcr_missing = metadata$qpcr_missing,
  original_count_sum = colSums(count_matrix),
  analysis_count_sum = analysis_library_sums,
  analysis_count_sum_over_median = analysis_library_sums / median_analysis_library_sum,
  count_sum_change = analysis_library_sums - colSums(count_matrix),
  size_factor = sizeFactors(DESeq2Obj),
  cooks_max = cooks_max[metadata$sample_id],
  cooks_q99 = cooks_q99[metadata$sample_id],
  included_in_analysis = T,
  row.names=NULL,
  check.names=F
)
write.csv(sample_qc, file.path(output_dir, "sample_qc.csv"), row.names=F, quote=F, na="NA")

# Diagnostic plots use a validated two-category palette and retain a table-view
# companion through the exported PCA coordinates and complete result table.
plot_surface = "#fcfcfb"
plot_text = "#0b0b0b"
plot_muted = "#898781"
plot_grid = "#e1e0d9"
plot_axis = "#c3c2b7"
group_colors = c(pooled_control="#2a78d6", OCT4_targeting="#1baf7a")
group_shapes = c(pooled_control=21L, OCT4_targeting=24L)

png(
  filename=file.path(output_dir, "dispersion-estimates.png"),
  width=10,
  height=8,
  units="in",
  res=300,
  bg=plot_surface
)
par(fg=plot_text, col.axis=plot_text, col.lab=plot_text, col.main=plot_text)
plotDispEsts(DESeq2Obj)
invisible(dev.off())

vst_data = varianceStabilizingTransformation(DESeq2Obj, blind=F)
pca_data = plotPCA(vst_data, intgroup=c("comparison_group", "guide_class"), returnData=T)
percent_variance = round(100 * attr(pca_data, "percentVar"), 1)
pca_data$sample_id = rownames(pca_data)
pca_data$guide = as.character(metadata[pca_data$sample_id, "guide"])
pca_data$replicate = as.character(metadata[pca_data$sample_id, "replicate"])
pca_data$jy_guide_label = metadata[pca_data$sample_id, "jy_guide_label"]
pca_data$qpcr_oct4_post_treatment = metadata[pca_data$sample_id, "qpcr_oct4_post_treatment"]
pca_data = pca_data[, c(
  "sample_id", "guide", "jy_guide_label", "replicate", "comparison_group",
  "guide_class", "qpcr_oct4_post_treatment", "PC1", "PC2"
)]
write.csv(pca_data, file.path(output_dir, "pca-coordinates.csv"), row.names=F, quote=F, na="NA")

png(
  filename=file.path(output_dir, "pca-oct4-targeting-vs-pooled-controls.png"),
  width=10,
  height=8,
  units="in",
  res=300,
  bg=plot_surface
)
par(
  mar=c(5.2, 5.2, 4.5, 2),
  bg=plot_surface,
  fg=plot_text,
  col.axis=plot_text,
  col.lab=plot_text,
  col.main=plot_text,
  family="sans"
)
plot(
  pca_data$PC1,
  pca_data$PC2,
  type="n",
  xlab=sprintf("PC1: %.1f%% variance", percent_variance[[1L]]),
  ylab=sprintf("PC2: %.1f%% variance", percent_variance[[2L]]),
  main="K2LLPS variance-stabilized PCA",
  sub="Shape and color distinguish the pooled controls from OCT4-targeting libraries",
  axes=F
)
abline(h=axTicks(2), v=axTicks(1), col=plot_grid, lwd=0.7)
axis(1, col=plot_axis, col.ticks=plot_axis, col.axis=plot_text)
axis(2, col=plot_axis, col.ticks=plot_axis, col.axis=plot_text, las=1)
box(col=plot_axis)
for (group_name in names(group_colors)) {
  group_rows = as.character(pca_data$comparison_group) == group_name
  points(
    pca_data$PC1[group_rows],
    pca_data$PC2[group_rows],
    pch=group_shapes[[group_name]],
    bg=group_colors[[group_name]],
    col=plot_surface,
    lwd=2,
    cex=1.35
  )
}
qc_label_sample = "guide15_rep_4"
qc_label_row = pca_data$sample_id == qc_label_sample
if (any(qc_label_row)) {
  text(
    pca_data$PC1[qc_label_row],
    pca_data$PC2[qc_label_row],
    labels=qc_label_sample,
    pos=4,
    offset=0.7,
    cex=0.72,
    col=plot_text
  )
}
legend(
  "topright",
  legend=c("Pooled control", "OCT4 targeting"),
  pch=unname(group_shapes),
  pt.bg=unname(group_colors),
  col=plot_surface,
  pt.cex=1.3,
  bty="n",
  text.col=plot_text
)
invisible(dev.off())

# PCA loadings for the same variance-stabilized matrix and gene set used by
# DESeq2::plotPCA (top ntop genes by row variance). Positive PC1 loadings push
# samples rightward on the companion sample PCA plot.
pca_ntop = 500L
pca_loadings_top_n = 30L
vst_matrix = assay(vst_data)
row_variances = matrixStats::rowVars(vst_matrix)
names(row_variances) = rownames(vst_matrix)
finite_variance_genes = names(row_variances)[is.finite(row_variances)]
if (length(finite_variance_genes) == 0L) {
  stop_with("No genes with finite row variance were available for PCA loadings")
}
pca_gene_count = min(pca_ntop, length(finite_variance_genes))
pca_selected_genes = names(sort(row_variances[finite_variance_genes], decreasing=T))[seq_len(pca_gene_count)]
pca_fit = prcomp(t(vst_matrix[pca_selected_genes, , drop=F]), center=T, scale.=F)
pca_percent_variance = 100 * (pca_fit$sdev ^ 2) / sum(pca_fit$sdev ^ 2)

# Guard that this reconstruction matches the sample PCA already exported above.
pca_score_check = data.frame(
  sample_id = rownames(pca_fit$x),
  PC1 = pca_fit$x[, 1L],
  PC2 = pca_fit$x[, 2L],
  stringsAsFactors=F,
  check.names=F
)
pca_score_check = pca_score_check[match(pca_data$sample_id, pca_score_check$sample_id), , drop=F]
if (anyNA(pca_score_check$sample_id)) {
  stop_with("Internal PCA loadings reconstruction lost one or more sample IDs")
}
if (
  !isTRUE(all.equal(pca_score_check$PC1, pca_data$PC1, tolerance=1e-6, check.attributes=F)) ||
  !isTRUE(all.equal(pca_score_check$PC2, pca_data$PC2, tolerance=1e-6, check.attributes=F))
) {
  stop_with("Reconstructed PCA scores do not match plotPCA coordinates")
}
if (!isTRUE(all.equal(round(pca_percent_variance[1:2], 1), percent_variance, tolerance=1e-8, check.attributes=F))) {
  stop_with("Reconstructed PCA percent variance does not match plotPCA labels")
}

pca_gene_names = as.character(mcols(DESeq2Obj)[pca_selected_genes, "gene_name"])
missing_gene_name = is.na(pca_gene_names) | pca_gene_names == ""
pca_gene_names[missing_gene_name] = pca_selected_genes[missing_gene_name]
pca_loadings_table = data.frame(
  gene_id = pca_selected_genes,
  gene_name = pca_gene_names,
  gene_biotype = as.character(mcols(DESeq2Obj)[pca_selected_genes, "gene_biotype"]),
  row_variance = unname(row_variances[pca_selected_genes]),
  PC1_loading = unname(pca_fit$rotation[, 1L]),
  PC2_loading = unname(pca_fit$rotation[, 2L]),
  stringsAsFactors=F,
  check.names=F
)
pca_loadings_table$abs_PC1_loading = abs(pca_loadings_table$PC1_loading)
pca_loadings_table$abs_PC2_loading = abs(pca_loadings_table$PC2_loading)
pca_loadings_table$PC_loading_radius = sqrt(
  pca_loadings_table$PC1_loading ^ 2 + pca_loadings_table$PC2_loading ^ 2
)
pca_loadings_table$PC1_direction = ifelse(
  pca_loadings_table$PC1_loading > 0,
  "positive_PC1",
  ifelse(pca_loadings_table$PC1_loading < 0, "negative_PC1", "zero_PC1")
)
pca_loadings_table$PC2_direction = ifelse(
  pca_loadings_table$PC2_loading > 0,
  "positive_PC2",
  ifelse(pca_loadings_table$PC2_loading < 0, "negative_PC2", "zero_PC2")
)
# Keep the full loadings table sorted by decreasing PC1 for a stable primary export.
pca_loadings_table = pca_loadings_table[order(-pca_loadings_table$PC1_loading, pca_loadings_table$gene_id), , drop=F]
rownames(pca_loadings_table) = NULL

pca_loadings_path = file.path(output_dir, "pca-loadings-oct4-targeting-vs-pooled-controls.csv")
write.csv(pca_loadings_table, pca_loadings_path, row.names=F, quote=F, na="NA")

make_loading_display_label = function(gene_name, gene_id) {
  duplicate_gene_name = duplicated(gene_name) | duplicated(gene_name, fromLast=T)
  ifelse(
    gene_name == gene_id,
    gene_id,
    ifelse(duplicate_gene_name, paste0(gene_name, " (", gene_id, ")"), gene_name)
  )
}

# For each axis, export the top/bottom contributors and draw both a 1D bar plot
# and a native PC1-PC2 arrow plot using that same gene set.
pca_axis_specs = list(
  list(
    axis="PC1",
    loading_column="PC1_loading",
    percent_index=1L,
    extremes_csv="pca-pc1-loadings-top-bottom-30-oct4-targeting-vs-pooled-controls.csv",
    extremes_png="pca-pc1-loadings-top-bottom-30-oct4-targeting-vs-pooled-controls.png",
    arrows_csv="pca-pc1-eigenvector-arrows-oct4-targeting-vs-pooled-controls.csv",
    arrows_png="pca-pc1-eigenvector-arrows-oct4-targeting-vs-pooled-controls.png",
    bar_main="K2LLPS PCA PC1 loadings",
    arrow_title="K2LLPS PCA gene eigenvectors (PC1 extremes)",
    positive_legend="Top 30 positive PC1 loading",
    negative_legend="Bottom 30 negative PC1 loading"
  ),
  list(
    axis="PC2",
    loading_column="PC2_loading",
    percent_index=2L,
    extremes_csv="pca-pc2-loadings-top-bottom-30-oct4-targeting-vs-pooled-controls.csv",
    extremes_png="pca-pc2-loadings-top-bottom-30-oct4-targeting-vs-pooled-controls.png",
    arrows_csv="pca-pc2-eigenvector-arrows-oct4-targeting-vs-pooled-controls.csv",
    arrows_png="pca-pc2-eigenvector-arrows-oct4-targeting-vs-pooled-controls.png",
    bar_main="K2LLPS PCA PC2 loadings",
    arrow_title="K2LLPS PCA gene eigenvectors (PC2 extremes)",
    positive_legend="Top 30 positive PC2 loading",
    negative_legend="Bottom 30 negative PC2 loading"
  )
)

for (axis_spec in pca_axis_specs) {
  loading_values_all = pca_loadings_table[[axis_spec$loading_column]]
  axis_order = order(-loading_values_all, pca_loadings_table$gene_id)
  axis_sorted = pca_loadings_table[axis_order, , drop=F]
  positive_count = min(pca_loadings_top_n, sum(loading_values_all > 0))
  negative_count = min(pca_loadings_top_n, sum(loading_values_all < 0))
  if (positive_count == 0L || negative_count == 0L) {
    stop_with("%s loadings do not contain both positive and negative contributors", axis_spec$axis)
  }

  pca_loadings_extremes = rbind(
    head(axis_sorted, positive_count),
    tail(axis_sorted, negative_count)
  )
  pca_loadings_extremes = pca_loadings_extremes[
    order(pca_loadings_extremes[[axis_spec$loading_column]], pca_loadings_extremes$gene_id),
    ,
    drop=F
  ]
  rownames(pca_loadings_extremes) = NULL
  pca_loadings_extremes$display_label = make_loading_display_label(
    pca_loadings_extremes$gene_name,
    pca_loadings_extremes$gene_id
  )
  write.csv(
    pca_loadings_extremes,
    file.path(output_dir, axis_spec$extremes_csv),
    row.names=F,
    quote=F,
    na="NA"
  )

  loading_values = pca_loadings_extremes[[axis_spec$loading_column]]
  loading_colors = ifelse(
    loading_values >= 0,
    group_colors[["OCT4_targeting"]],
    group_colors[["pooled_control"]]
  )
  loading_limit = max(abs(loading_values))
  if (!is.finite(loading_limit) || loading_limit <= 0) loading_limit = 1

  png(
    filename=file.path(output_dir, axis_spec$extremes_png),
    width=11,
    height=14,
    units="in",
    res=300,
    bg=plot_surface
  )
  par(
    mar=c(5.2, 12.5, 4.8, 1.8),
    bg=plot_surface,
    fg=plot_text,
    col.axis=plot_text,
    col.lab=plot_text,
    col.main=plot_text,
    family="sans"
  )
  barplot(
    loading_values,
    horiz=T,
    names.arg=pca_loadings_extremes$display_label,
    las=1,
    col=loading_colors,
    border=NA,
    xlim=c(-loading_limit * 1.08, loading_limit * 1.08),
    xlab=sprintf(
      "%s loading (%.1f%% variance; top %d variable genes)",
      axis_spec$axis,
      percent_variance[[axis_spec$percent_index]],
      pca_gene_count
    ),
    main=axis_spec$bar_main,
    sub=sprintf(
      "Bottom %d and top %d contributors among the same gene set used for the sample PCA",
      negative_count,
      positive_count
    ),
    axes=F,
    cex.names=0.72
  )
  abline(v=0, col=plot_axis, lwd=1.2)
  axis(1, col=plot_axis, col.ticks=plot_axis, col.axis=plot_text)
  box(col=plot_axis)
  legend(
    "bottomright",
    legend=c(
      sprintf("Positive %s loading", axis_spec$axis),
      sprintf("Negative %s loading", axis_spec$axis)
    ),
    fill=c(group_colors[["OCT4_targeting"]], group_colors[["pooled_control"]]),
    border=NA,
    bty="n",
    text.col=plot_text
  )
  invisible(dev.off())

  pca_arrows = pca_loadings_extremes
  pca_arrows$PC1_arrow = pca_arrows$PC1_loading
  pca_arrows$PC2_arrow = pca_arrows$PC2_loading
  pca_arrow_count = NROW(pca_arrows)
  max_loading_radius = max(pca_arrows$PC_loading_radius, na.rm=T)
  if (!is.finite(max_loading_radius) || max_loading_radius <= 0) {
    stop_with("Selected %s gene loadings have non-positive radius; cannot draw arrows", axis_spec$axis)
  }
  write.csv(pca_arrows, file.path(output_dir, axis_spec$arrows_csv), row.names=F, quote=F, na="NA")

  x_range = range(c(pca_arrows$PC1_arrow, 0), finite=T)
  y_range = range(c(pca_arrows$PC2_arrow, 0), finite=T)
  x_pad = max(diff(x_range) * 0.18, max_loading_radius * 0.12)
  y_pad = max(diff(y_range) * 0.18, max_loading_radius * 0.12)
  if (!is.finite(x_pad) || x_pad == 0) x_pad = 0.01
  if (!is.finite(y_pad) || y_pad == 0) y_pad = 0.01

  direction_values = pca_arrows[[axis_spec$loading_column]]
  pca_arrows$direction_label = factor(
    ifelse(direction_values >= 0, axis_spec$positive_legend, axis_spec$negative_legend),
    levels=c(axis_spec$positive_legend, axis_spec$negative_legend)
  )
  pca_arrow_direction_colors = c(
    stats::setNames(unname(group_colors[["OCT4_targeting"]]), axis_spec$positive_legend),
    stats::setNames(unname(group_colors[["pooled_control"]]), axis_spec$negative_legend)
  )

  pca_arrow_plot = ggplot2::ggplot(
    pca_arrows,
    ggplot2::aes(x=PC1_arrow, y=PC2_arrow, color=direction_label)
  ) +
    ggplot2::geom_hline(yintercept=0, color=plot_axis, linewidth=0.45) +
    ggplot2::geom_vline(xintercept=0, color=plot_axis, linewidth=0.45) +
    ggplot2::geom_segment(
      ggplot2::aes(x=0, y=0, xend=PC1_arrow, yend=PC2_arrow),
      arrow=ggplot2::arrow(length=grid::unit(0.18, "cm"), type="closed"),
      linewidth=0.7,
      alpha=0.9,
      lineend="round"
    ) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label=display_label),
      size=2.7,
      fontface="plain",
      family="sans",
      box.padding=0.35,
      point.padding=0.2,
      min.segment.length=0,
      segment.size=0.25,
      segment.alpha=0.55,
      segment.color=plot_muted,
      force=1.4,
      force_pull=0.6,
      max.overlaps=Inf,
      max.time=2,
      max.iter=20000,
      seed=1,
      show.legend=F
    ) +
    ggplot2::scale_color_manual(values=pca_arrow_direction_colors, name=NULL) +
    ggplot2::coord_cartesian(
      xlim=c(x_range[[1L]] - x_pad, x_range[[2L]] + x_pad),
      ylim=c(y_range[[1L]] - y_pad, y_range[[2L]] + y_pad),
      expand=F,
      clip="off"
    ) +
    ggplot2::labs(
      title=axis_spec$arrow_title,
      subtitle=sprintf(
        "Top %d and bottom %d %s loadings (same gene set as the %s bar plot); native coordinates from the VST / top-%d gene PCA",
        positive_count,
        negative_count,
        axis_spec$axis,
        axis_spec$axis,
        pca_gene_count
      ),
      x=sprintf("PC1 loading (%.1f%% variance)", percent_variance[[1L]]),
      y=sprintf("PC2 loading (%.1f%% variance)", percent_variance[[2L]])
    ) +
    ggplot2::theme_bw(base_family="sans") +
    ggplot2::theme(
      plot.background=ggplot2::element_rect(fill=plot_surface, color=NA),
      panel.background=ggplot2::element_rect(fill=plot_surface, color=NA),
      panel.border=ggplot2::element_rect(color=plot_axis, fill=NA, linewidth=0.6),
      panel.grid.major=ggplot2::element_line(color=plot_grid, linewidth=0.35),
      panel.grid.minor=ggplot2::element_blank(),
      axis.ticks=ggplot2::element_line(color=plot_axis),
      axis.text=ggplot2::element_text(color=plot_text, size=10),
      axis.title=ggplot2::element_text(color=plot_text, size=11),
      plot.title=ggplot2::element_text(color=plot_text, size=14, face="plain", hjust=0.5),
      plot.subtitle=ggplot2::element_text(color=plot_muted, size=9, hjust=0.5),
      legend.position="top",
      legend.justification="left",
      legend.background=ggplot2::element_blank(),
      legend.key=ggplot2::element_blank(),
      legend.text=ggplot2::element_text(color=plot_text, size=10),
      plot.margin=ggplot2::margin(12, 18, 12, 12)
    )

  ggplot2::ggsave(
    filename=file.path(output_dir, axis_spec$arrows_png),
    plot=pca_arrow_plot,
    width=11,
    height=9,
    units="in",
    dpi=300,
    bg=plot_surface
  )
}

# Sample PC1 scores (gene-loading projection) aggregated by gRNA species:
# 8 OCT4-targeting guides + 1 pooled-control group; x-axis sorted by median PC1.
pc1_by_species = data.frame(
  sample_id = pca_data$sample_id,
  guide = as.character(pca_data$guide),
  comparison_group = as.character(pca_data$comparison_group),
  PC1 = pca_data$PC1,
  stringsAsFactors=F,
  check.names=F
)
pc1_by_species$gRNA_species = ifelse(
  pc1_by_species$comparison_group == "pooled_control",
  "pooled_control",
  pc1_by_species$guide
)
pc1_species_levels_expected = c(
  "pooled_control", "guide01", "guide02", "guide03", "guide04",
  "guide05", "guide06", "guide14", "guide15"
)
if (!setequal(pc1_by_species$gRNA_species, pc1_species_levels_expected)) {
  stop_with("Unexpected gRNA species while aggregating PC1 scores")
}
pc1_species_medians = tapply(pc1_by_species$PC1, pc1_by_species$gRNA_species, stats::median)
pc1_species_order = names(sort(pc1_species_medians, decreasing=F))
pc1_by_species$gRNA_species = factor(pc1_by_species$gRNA_species, levels=pc1_species_order)
pc1_species_colors = c(
  pooled_control=unname(group_colors[["pooled_control"]]),
  guide01="#74c476", guide02="#41ab5d", guide03="#238b45", guide04="#006d2c",
  guide05="#00441b", guide06="#66c2a4", guide14="#2ca25f", guide15="#99d8c9"
)
pc1_species_table = data.frame(
  sample_id = pc1_by_species$sample_id,
  gRNA_species = as.character(pc1_by_species$gRNA_species),
  guide = pc1_by_species$guide,
  comparison_group = pc1_by_species$comparison_group,
  PC1 = pc1_by_species$PC1,
  stringsAsFactors=F,
  check.names=F
)
pc1_species_table = pc1_species_table[
  order(match(pc1_species_table$gRNA_species, pc1_species_order), pc1_species_table$sample_id),
  ,
  drop=F
]
rownames(pc1_species_table) = NULL
write.csv(
  pc1_species_table,
  file.path(output_dir, "pca-pc1-score-by-grna-species.csv"),
  row.names=F,
  quote=F,
  na="NA"
)
pc1_species_median_table = data.frame(
  gRNA_species = pc1_species_order,
  n_samples = as.integer(table(pc1_by_species$gRNA_species)[pc1_species_order]),
  PC1_median = as.numeric(pc1_species_medians[pc1_species_order]),
  PC1_mean = as.numeric(tapply(pc1_by_species$PC1, pc1_by_species$gRNA_species, mean)[pc1_species_order]),
  stringsAsFactors=F,
  check.names=F
)
write.csv(
  pc1_species_median_table,
  file.path(output_dir, "pca-pc1-score-by-grna-species-medians.csv"),
  row.names=F,
  quote=F,
  na="NA"
)
pc1_species_violin = ggplot2::ggplot(
  pc1_by_species,
  ggplot2::aes(x=gRNA_species, y=PC1, fill=gRNA_species)
) +
  ggplot2::geom_hline(yintercept=0, color=plot_axis, linewidth=0.45) +
  ggplot2::geom_violin(scale="width", color=plot_text, linewidth=0.3, alpha=0.9, width=0.9) +
  ggplot2::geom_boxplot(width=0.12, outlier.shape=NA, fill=plot_surface, color=plot_text, linewidth=0.3) +
  ggplot2::geom_point(
    position=ggplot2::position_jitter(width=0.08, height=0, seed=1),
    size=1.8, alpha=0.85, color=plot_text
  ) +
  ggplot2::stat_summary(fun=stats::median, geom="point", size=2.2, shape=23, fill="#f7f7f7", color=plot_text) +
  ggplot2::scale_fill_manual(values=pc1_species_colors, guide="none") +
  ggplot2::labs(
    title="PC1 scores aggregated by gRNA species",
    subtitle="Sample PC1 values from the VST PCA (gene-loading projection); x-axis sorted by median PC1",
    x=NULL,
    y="PC1 score"
  ) +
  ggplot2::theme_bw(base_family="sans") +
  ggplot2::theme(
    plot.background=ggplot2::element_rect(fill=plot_surface, color=NA),
    panel.background=ggplot2::element_rect(fill=plot_surface, color=NA),
    panel.border=ggplot2::element_rect(color=plot_axis, fill=NA, linewidth=0.6),
    panel.grid.major.x=ggplot2::element_blank(),
    panel.grid.minor=ggplot2::element_blank(),
    panel.grid.major.y=ggplot2::element_line(color=plot_grid, linewidth=0.35),
    axis.ticks=ggplot2::element_line(color=plot_axis),
    axis.text.x=ggplot2::element_text(color=plot_text, size=10, angle=45, hjust=1),
    axis.text.y=ggplot2::element_text(color=plot_text, size=10),
    axis.title=ggplot2::element_text(color=plot_text, size=11),
    plot.title=ggplot2::element_text(color=plot_text, size=14, hjust=0.5),
    plot.subtitle=ggplot2::element_text(color=plot_muted, size=9, hjust=0.5),
    plot.margin=ggplot2::margin(12, 14, 12, 12)
  )
ggplot2::ggsave(
  filename=file.path(output_dir, "pca-pc1-score-violin-by-grna-species.png"),
  plot=pc1_species_violin,
  width=10,
  height=6.5,
  units="in",
  dpi=300,
  bg=plot_surface
)

finite_ma_values = c(result_table$log2FoldChange_raw, result_table$log2FoldChange_shrunken)
finite_ma_values = finite_ma_values[is.finite(finite_ma_values)]
ma_limit = if (length(finite_ma_values) == 0L) 1 else max(1, max(abs(finite_ma_values)))
ma_x = log10(result_table$baseMean + 1)
ma_significant_color = "#2a78d6"
ma_nonsignificant_color = grDevices::adjustcolor(plot_muted, alpha.f=0.38)

png(
  filename=file.path(output_dir, "ma-oct4-targeting-vs-pooled-controls.png"),
  width=26,
  height=13,
  units="in",
  res=300,
  bg=plot_surface
)
par(
  mfrow=c(1, 2),
  mar=c(5.2, 5.2, 4.2, 1.5),
  oma=c(0, 0, 2.5, 0),
  bg=plot_surface,
  fg=plot_text,
  col.axis=plot_text,
  col.lab=plot_text,
  col.main=plot_text,
  family="sans"
)
ma_panels = list(
  "Raw log2 fold change"=result_table$log2FoldChange_raw,
  "Shrunken log2 fold change"=result_table$log2FoldChange_shrunken
)
for (panel_name in names(ma_panels)) {
  panel_values = ma_panels[[panel_name]]
  plot(
    ma_x,
    panel_values,
    type="n",
    xlab="log10(baseMean + 1)",
    ylab="log2 fold change: OCT4 targeting minus pooled control",
    main=panel_name,
    ylim=c(-ma_limit, ma_limit),
    axes=F
  )
  abline(h=axTicks(2), v=axTicks(1), col=plot_grid, lwd=0.7)
  abline(h=0, col=plot_axis, lwd=1.2)
  axis(1, col=plot_axis, col.ticks=plot_axis, col.axis=plot_text)
  axis(2, col=plot_axis, col.ticks=plot_axis, col.axis=plot_text, las=1)
  box(col=plot_axis)
  nonsignificant_rows = !result_table$significant_fdr & is.finite(panel_values)
  significant_rows = result_table$significant_fdr & is.finite(panel_values)
  points(ma_x[nonsignificant_rows], panel_values[nonsignificant_rows], pch=16, cex=0.28, col=ma_nonsignificant_color)
  points(ma_x[significant_rows], panel_values[significant_rows], pch=16, cex=0.42, col=ma_significant_color)
  legend(
    "topright",
    legend=c(sprintf("FDR < %.3g", alpha), "Not FDR-significant"),
    pch=16,
    col=c(ma_significant_color, plot_muted),
    bty="n",
    pt.cex=c(0.9, 0.8),
    text.col=plot_text
  )
}
mtext("K2LLPS differential expression", outer=T, side=3, line=0.5, cex=1.15, col=plot_text)
invisible(dev.off())

# Strict DEG set for visualization: padj threshold plus raw |log2FC| cutoff.
# Comparison remains OCT4-targeting (8 guides) vs pooled controls (guide18/19/20).
result_table$deg_strict = (
  !is.na(result_table$padj) &
  result_table$padj <= alpha &
  !is.na(result_table$log2FoldChange_raw) &
  abs(result_table$log2FoldChange_raw) >= deg_lfc_threshold
)
result_table$deg_strict_direction = ifelse(
  !result_table$deg_strict,
  "not_deg",
  ifelse(
    result_table$log2FoldChange_raw > 0,
    "up_in_OCT4_targeting",
    "down_in_OCT4_targeting"
  )
)
deg_table = result_table[result_table$deg_strict, , drop=F]
deg_up_table = deg_table[deg_table$deg_strict_direction == "up_in_OCT4_targeting", , drop=F]
deg_down_table = deg_table[deg_table$deg_strict_direction == "down_in_OCT4_targeting", , drop=F]
deg_up_count = NROW(deg_up_table)
deg_down_count = NROW(deg_down_table)
deg_total_count = NROW(deg_table)
if (deg_total_count == 0L) {
  stop_with(
    "No DEGs passed padj <= %s and |raw log2FC| >= %s",
    format(alpha, scientific=F),
    format(deg_lfc_threshold, scientific=F)
  )
}
deg_table_path = file.path(
  output_dir,
  "degs-oct4-targeting-vs-pooled-controls-lfc1-padj0.05.csv"
)
write.csv(deg_table, deg_table_path, row.names=F, quote=F, na="NA")

# Row-z heatmap on VST expression for the strict DEG set.
vst_deg_matrix = assay(vst_data)[deg_table$gene_id, metadata$sample_id, drop=F]
row_sd = matrixStats::rowSds(vst_deg_matrix)
keep_variable_deg = is.finite(row_sd) & row_sd > 0
if (sum(keep_variable_deg) < 2L) {
  stop_with("Fewer than 2 DEGs have positive VST variance for heatmap clustering")
}
if (any(!keep_variable_deg)) {
  message(sprintf(
    "Excluding %d DEG(s) with zero VST variance from the heatmap/violin z-score matrix",
    sum(!keep_variable_deg)
  ))
}
vst_deg_matrix = vst_deg_matrix[keep_variable_deg, , drop=F]
deg_table_plot = deg_table[match(rownames(vst_deg_matrix), deg_table$gene_id), , drop=F]
deg_z_matrix = t(scale(t(vst_deg_matrix)))
deg_z_matrix[!is.finite(deg_z_matrix)] = 0

row_cor = stats::cor(t(deg_z_matrix), use="pairwise.complete.obs")
if (any(!is.finite(row_cor))) {
  stop_with("Correlation distance for DEG heatmap contains non-finite values")
}
row_dist = stats::as.dist(1 - row_cor)
row_hclust = stats::hclust(row_dist, method="ward.D2")

# Column order follows the PC1-median gRNA-species order, then replicate.
target_guide_levels = c("guide01", "guide02", "guide03", "guide04", "guide05", "guide06", "guide14", "guide15")
control_guide_levels = c("guide18", "guide19", "guide20")
heatmap_species = ifelse(
  as.character(metadata$comparison_group) == "pooled_control",
  "pooled_control",
  as.character(metadata$guide)
)
if (!identical(sort(pc1_species_order), sort(c("pooled_control", target_guide_levels)))) {
  stop_with("PC1 species order does not match the expected 9 gRNA-species labels")
}
# Within each gRNA species: replicate, then underlying guide id (controls are
# guide18/19/20 under the single pooled_control species), then sample_id.
sample_order = metadata$sample_id[
  order(
    match(heatmap_species, pc1_species_order),
    suppressWarnings(as.integer(as.character(metadata$replicate))),
    match(as.character(metadata$guide), c(control_guide_levels, target_guide_levels)),
    metadata$sample_id
  )
]
deg_z_matrix = deg_z_matrix[, sample_order, drop=F]
# Keep heatmap annotations aligned with the 9-species violin/PC1 order.
# Do not add a separate per-library `guide` track: for OCT4 guides it duplicates
# gRNA_species, while for controls it splits pooled_control into guide18/19/20
# and looks inconsistent with the right-hand species blocks.
heatmap_annotation = data.frame(
  replicate = as.character(metadata[sample_order, "replicate"]),
  gRNA_species = heatmap_species[match(sample_order, metadata$sample_id)],
  comparison_group = as.character(metadata[sample_order, "comparison_group"]),
  row.names = sample_order,
  check.names=F
)
heatmap_annotation$gRNA_species = factor(heatmap_annotation$gRNA_species, levels=pc1_species_order)
heatmap_annotation$comparison_group = factor(
  heatmap_annotation$comparison_group,
  levels=c("pooled_control", "OCT4_targeting")
)
heatmap_annotation$replicate = factor(heatmap_annotation$replicate, levels=as.character(2:5))
annotation_colors = list(
  replicate = c("2"="#f0f0f0", "3"="#bdbdbd", "4"="#969696", "5"="#636363"),
  gRNA_species = c(
    pooled_control=unname(group_colors[["pooled_control"]]),
    guide01="#74c476", guide02="#41ab5d", guide03="#238b45", guide04="#006d2c",
    guide05="#00441b", guide06="#66c2a4", guide14="#2ca25f", guide15="#99d8c9"
  )[pc1_species_order],
  comparison_group = c(
    pooled_control=unname(group_colors[["pooled_control"]]),
    OCT4_targeting=unname(group_colors[["OCT4_targeting"]])
  )
)
z_limit = min(3, max(abs(deg_z_matrix), na.rm=T))
if (!is.finite(z_limit) || z_limit <= 0) z_limit = 3
heatmap_breaks = seq(-z_limit, z_limit, length.out=101)
heatmap_colors = grDevices::colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(100)

deg_heatmap_path = file.path(
  output_dir,
  "deg-heatmap-vst-rowz-oct4-targeting-vs-pooled-controls.png"
)
pheatmap::pheatmap(
  deg_z_matrix,
  color=heatmap_colors,
  breaks=heatmap_breaks,
  cluster_rows=row_hclust,
  cluster_cols=F,
  show_rownames=F,
  show_colnames=F,
  annotation_col=heatmap_annotation,
  annotation_colors=annotation_colors,
  clustering_method="ward.D2",
  border_color=NA,
  fontsize=9,
  main=sprintf(
    "DEG VST row Z-scores (%d genes; |raw LFC| >= %s, padj <= %s)\nrow hclust: ward.D2 / 1-cor; columns: PC1-median gRNA order, then replicate",
    NROW(deg_z_matrix),
    format(deg_lfc_threshold, scientific=F),
    format(alpha, scientific=F)
  ),
  filename=deg_heatmap_path,
  width=12,
  height=10
)

# 18 violins: up/down DEGs x (8 OCT4 guides + 1 aggregated pooled control).
# Use the same PC1-median species order as the PC1 score violin / heatmap columns.
violin_guide_levels = pc1_species_order
make_deg_violin_df = function(gene_ids, direction_label) {
  if (length(gene_ids) == 0L) {
    return(data.frame(
      gene_id=character(),
      sample_id=character(),
      row_z=numeric(),
      violin_group=character(),
      direction=character(),
      stringsAsFactors=F,
      check.names=F
    ))
  }
  sub_z = deg_z_matrix[gene_ids, , drop=F]
  long = data.frame(
    gene_id = rep(rownames(sub_z), times=NCOL(sub_z)),
    sample_id = rep(colnames(sub_z), each=NROW(sub_z)),
    row_z = as.numeric(sub_z),
    stringsAsFactors=F,
    check.names=F
  )
  long$guide = as.character(metadata[long$sample_id, "guide"])
  long$comparison_group = as.character(metadata[long$sample_id, "comparison_group"])
  long$violin_group = ifelse(
    long$comparison_group == "pooled_control",
    "pooled_control",
    long$guide
  )
  long$direction = direction_label
  long[, c("gene_id", "sample_id", "row_z", "violin_group", "direction"), drop=F]
}

up_gene_ids = deg_table_plot$gene_id[deg_table_plot$deg_strict_direction == "up_in_OCT4_targeting"]
down_gene_ids = deg_table_plot$gene_id[deg_table_plot$deg_strict_direction == "down_in_OCT4_targeting"]
violin_df = rbind(
  make_deg_violin_df(up_gene_ids, sprintf("Up in OCT4 targeting (n=%d)", length(up_gene_ids))),
  make_deg_violin_df(down_gene_ids, sprintf("Down in OCT4 targeting (n=%d)", length(down_gene_ids)))
)
if (NROW(violin_df) == 0L) {
  stop_with("No DEG values were available for violin plots")
}
violin_df$violin_group = factor(violin_df$violin_group, levels=violin_guide_levels)
violin_df$direction = factor(
  violin_df$direction,
  levels=c(
    sprintf("Up in OCT4 targeting (n=%d)", length(up_gene_ids)),
    sprintf("Down in OCT4 targeting (n=%d)", length(down_gene_ids))
  )
)
violin_df_path = file.path(
  output_dir,
  "deg-violin-data-oct4-targeting-vs-pooled-controls.csv"
)
write.csv(violin_df, violin_df_path, row.names=F, quote=F, na="NA")

violin_colors = c(
  pooled_control=unname(group_colors[["pooled_control"]]),
  guide01="#74c476", guide02="#41ab5d", guide03="#238b45", guide04="#006d2c",
  guide05="#00441b", guide06="#66c2a4", guide14="#2ca25f", guide15="#99d8c9"
)
deg_violin_path = file.path(
  output_dir,
  "deg-violin-oct4-targeting-vs-pooled-controls.png"
)
deg_violin_plot = ggplot2::ggplot(
  violin_df,
  ggplot2::aes(x=violin_group, y=row_z, fill=violin_group)
) +
  ggplot2::geom_hline(yintercept=0, color=plot_axis, linewidth=0.4) +
  ggplot2::geom_violin(scale="width", color=plot_text, linewidth=0.25, alpha=0.9) +
  ggplot2::stat_summary(fun=median, geom="point", size=1.4, color=plot_text) +
  ggplot2::facet_wrap(~direction, nrow=2, scales="free_x") +
  ggplot2::scale_fill_manual(values=violin_colors, guide="none") +
  ggplot2::labs(
    title="DEG VST row Z-score distributions",
    subtitle=sprintf(
      "18 violins = up/down DEGs x (8 OCT4 guides + pooled controls); x-axis: PC1-median gRNA order; DEGs: |raw LFC| >= %s and padj <= %s",
      format(deg_lfc_threshold, scientific=F),
      format(alpha, scientific=F)
    ),
    x=NULL,
    y="Row Z-score (VST)"
  ) +
  ggplot2::theme_bw(base_family="sans") +
  ggplot2::theme(
    plot.background=ggplot2::element_rect(fill=plot_surface, color=NA),
    panel.background=ggplot2::element_rect(fill=plot_surface, color=NA),
    panel.border=ggplot2::element_rect(color=plot_axis, fill=NA, linewidth=0.6),
    panel.grid.major.x=ggplot2::element_blank(),
    panel.grid.minor=ggplot2::element_blank(),
    panel.grid.major.y=ggplot2::element_line(color=plot_grid, linewidth=0.35),
    axis.ticks=ggplot2::element_line(color=plot_axis),
    axis.text.x=ggplot2::element_text(color=plot_text, size=9, angle=45, hjust=1),
    axis.text.y=ggplot2::element_text(color=plot_text, size=10),
    axis.title=ggplot2::element_text(color=plot_text, size=11),
    plot.title=ggplot2::element_text(color=plot_text, size=14, hjust=0.5),
    plot.subtitle=ggplot2::element_text(color=plot_muted, size=9, hjust=0.5),
    strip.background=ggplot2::element_rect(fill="#f0efea", color=plot_axis),
    strip.text=ggplot2::element_text(color=plot_text, size=10),
    plot.margin=ggplot2::margin(12, 14, 12, 12)
  )
ggplot2::ggsave(
  filename=deg_violin_path,
  plot=deg_violin_plot,
  width=12,
  height=9,
  units="in",
  dpi=300,
  bg=plot_surface
)

writeLines(
  capture.output(devtools::session_info(pkgs="attached")),
  file.path(output_dir, "sessionInfo.txt")
)

significant_up_count = sum(result_table$significant_up)
significant_down_count = sum(result_table$significant_down)
summary_lines = c(
  "K2LLPS DESeq2 pooled-comparison summary",
  sprintf("Expression input: %s", input_path),
  sprintf("Expression input MD5: %s", md5_value(input_path)),
  sprintf("Metadata input: %s", metadata_path),
  sprintf("Metadata input MD5: %s", md5_value(metadata_path)),
  sprintf("Design workbook: %s", design_workbook_path),
  sprintf("Design workbook MD5: %s", md5_value(design_workbook_path)),
  sprintf("Workbook sheet: %s", design_sheet_name),
  sprintf("Workbook/sample keys matched: %d", NROW(metadata)),
  sprintf("Design: %s", paste(deparse(design_formula), collapse=" ")),
  "Replicate interpretation: shared experimental block; modeled as a categorical factor",
  sprintf("Model matrix: %d rows x %d columns; rank %d; residual df %d", NROW(model_matrix), NCOL(model_matrix), model_rank, residual_degrees_freedom),
  sprintf("Coefficient: %s", coefficient_name),
  "Coefficient orientation: positive log2 fold change means higher expression in OCT4-targeting libraries",
  sprintf("Targeting guides/libraries: %d / %d", 8L, sum(metadata$comparison_group == "OCT4_targeting")),
  sprintf("Control guides/libraries: %d / %d", 3L, sum(metadata$comparison_group == "pooled_control")),
  sprintf("qPCR missing values: %d (%s)", length(missing_qpcr_samples), paste(missing_qpcr_samples, collapse=", ")),
  "qPCR role: post-treatment annotation only; not included in the DESeq2 design",
  sprintf("Count strategy: %s", count_strategy),
  sprintf("Dispersion fit: %s", dispersion_fit_type),
  sprintf("Samples retained: %d / 44", NCOL(DESeq2Obj)),
  sprintf("guide15_rep_4 retained: %s", "guide15_rep_4" %in% colnames(DESeq2Obj)),
  sprintf("Size-factor range: %.6f to %.6f", min(sizeFactors(DESeq2Obj)), max(sizeFactors(DESeq2Obj))),
  sprintf("Genes before prefilter: %d", genes_before_filter),
  sprintf("Genes after prefilter: %d", genes_after_filter),
  sprintf("Prefilter: total analysis count >= %s", format(min_total_count, scientific=F)),
  sprintf("Fractional cells in supplied counts: %d (%.4f%%)", fractional_cells, 100 * fractional_cells / length(count_matrix)),
  sprintf("Wald-test alpha: %.17g", alpha),
  sprintf("Strict DEG absolute raw LFC cutoff: %.17g", deg_lfc_threshold),
  "Independent filtering: enabled",
  "Cook's-distance filtering: enabled at results extraction",
  "LFC shrinkage: DESeq2 type='normal' on the named targeting-versus-control coefficient",
  sprintf("Genes with a Wald p-value: %d", sum(result_table$is_tested)),
  sprintf("Genes with an adjusted p-value: %d", sum(result_table$has_adjusted_pvalue)),
  sprintf("Genes removed by independent filtering: %d", sum(result_table$independent_filter_flag)),
  sprintf("Genes with missing p-value after nonzero expression: %d", sum(result_table$cooks_or_other_filter_flag)),
  sprintf("FDR-significant genes: %d", sum(result_table$significant_fdr)),
  sprintf("FDR-significant higher in OCT4 targeting: %d", significant_up_count),
  sprintf("FDR-significant lower in OCT4 targeting: %d", significant_down_count),
  sprintf(
    "Strict DEGs (|raw LFC| >= %s and padj <= %s): %d (up %d, down %d)",
    format(deg_lfc_threshold, scientific=F),
    format(alpha, scientific=F),
    deg_total_count,
    deg_up_count,
    deg_down_count
  ),
  "Strict DEG LFC source: raw (unshrunken) log2 fold change",
  sprintf("Strict DEG table: %s", deg_table_path),
  sprintf("DEG heatmap: %s", deg_heatmap_path),
  sprintf("DEG violin plot: %s", deg_violin_path),
  sprintf("Complete results: %s", complete_results_path),
  sprintf("Significant results: %s", significant_results_path),
  "Status: size factors, dispersions, Wald test, result extraction, LFC shrinkage, and diagnostics completed."
)
writeLines(summary_lines, file.path(output_dir, "run-summary.txt"))

message("Done. Differential-expression testing completed.")
message(sprintf(
  "FDR < %.3g: %d genes (%d higher, %d lower in OCT4 targeting)",
  alpha,
  sum(result_table$significant_fdr),
  significant_up_count,
  significant_down_count
))
message(sprintf(
  "Strict DEGs (|raw LFC| >= %s, padj <= %s): %d (up %d, down %d)",
  format(deg_lfc_threshold, scientific=F),
  format(alpha, scientific=F),
  deg_total_count,
  deg_up_count,
  deg_down_count
))
message("Complete results: ", complete_results_path)
message("Significant results: ", significant_results_path)
message("Strict DEG table: ", deg_table_path)
message("DEG heatmap: ", deg_heatmap_path)
message("DEG violin plot: ", deg_violin_path)
