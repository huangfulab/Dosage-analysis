#!/usr/bin/env Rscript

# Fail-fast verification for the K2LLPS guide-specific blocked DESeq2 analysis.

options(stringsAsFactors=F)
if (!requireNamespace("devtools", quietly=T)) stop("devtools is required for session information", call.=F)

script_dir = function() {
  file_arg = grep("^--file=", commandArgs(trailingOnly=F), value=T)
  if (length(file_arg) == 0L) return(normalizePath(getwd(), mustWork=T))

  script_file = sub("^--file=", "", file_arg[[1L]])
  return(dirname(normalizePath(script_file, mustWork=T)))
}
stop_with = function(...) stop(sprintf(...), call.=F)
assert_true = function(condition, ...) {
  if (length(condition) != 1L || is.na(condition) || !condition) stop_with(...)
}
assert_equal_numeric = function(observed, expected, label, tolerance=1e-8) {
  comparison = all.equal(
    as.numeric(observed),
    as.numeric(expected),
    tolerance=tolerance,
    check.attributes=F
  )
  if (!isTRUE(comparison)) stop_with("%s mismatch: %s", label, paste(comparison, collapse="; "))
}
finite_summary = function(values, summary_function) {
  values = values[is.finite(values)]
  if (length(values) == 0L) return(NA_real_)
  summary_function(values)
}

args = commandArgs(trailingOnly=T)
project_dir = script_dir()
pooled_output_dir = if (length(args) >= 1L) args[[1L]] else file.path(project_dir, "outputs")
guide_output_dir = if (length(args) >= 2L) args[[2L]] else file.path(pooled_output_dir, "guide-specific")

if (!requireNamespace("DESeq2", quietly=T)) stop_with("DESeq2 is required but is not installed")
suppressPackageStartupMessages(library(DESeq2))

pooled_output_dir = normalizePath(pooled_output_dir, mustWork=T)
guide_output_dir = normalizePath(guide_output_dir, mustWork=T)

target_guides = c("guide01", "guide02", "guide03", "guide04", "guide05", "guide06", "guide14", "guide15")
control_guides = c("guide18", "guide19", "guide20")
expected_guides = c(target_guides, control_guides)
expected_replicates = as.character(2:5)
expected_coefficients = paste0("contrast_group_", target_guides, "_vs_pooled_control")

pooled_files = c(
  "analysis-sample-metadata.csv",
  "sample_qc.csv",
  "DESeq2Obj-after-wald-test.rds",
  "de-results-oct4-targeting-vs-pooled-controls-complete.csv",
  "de-results-oct4-targeting-vs-pooled-controls-significant.csv",
  "dispersion-estimates.png",
  "expression-file-equivalence.txt",
  "expression-header-map.csv",
  "pca-coordinates.csv",
  "pca-oct4-targeting-vs-pooled-controls.png",
  "ma-oct4-targeting-vs-pooled-controls.png",
  "run-summary.txt",
  "sessionInfo.txt"
)
expected_guide_files = c(
  "DESeq2Obj-guide-specific-after-wald-test.rds",
  unlist(lapply(target_guides, function(target_guide) c(
    sprintf("de-results-%s-vs-pooled-controls-complete.csv", target_guide),
    sprintf("de-results-%s-vs-pooled-controls-significant.csv", target_guide)
  ))),
  "concordance-across-guides.csv",
  "pooled-source-md5.csv",
  "run-summary.txt",
  "sessionInfo.txt"
)
expected_guide_paths = file.path(guide_output_dir, expected_guide_files)
missing_guide_files = expected_guide_files[!file.exists(expected_guide_paths)]
assert_true(length(missing_guide_files) == 0L, "Missing guide-specific output file(s): %s", paste(missing_guide_files, collapse=", "))
empty_guide_files = expected_guide_files[file.info(expected_guide_paths)$size <= 0]
assert_true(length(empty_guide_files) == 0L, "Empty guide-specific output file(s): %s", paste(empty_guide_files, collapse=", "))

source_manifest = read.csv(file.path(guide_output_dir, "pooled-source-md5.csv"), check.names=F)
assert_true(
  identical(names(source_manifest), c("file", "size_bytes", "md5")),
  "pooled-source-md5.csv has an unexpected schema"
)
assert_true(identical(source_manifest$file, pooled_files), "Pooled source manifest file order or membership is incorrect")
pooled_paths = file.path(pooled_output_dir, pooled_files)
assert_true(all(file.exists(pooled_paths)), "A pooled source file recorded by the guide-specific analysis is missing")
current_source_manifest = data.frame(
  file=pooled_files,
  size_bytes=as.numeric(file.info(pooled_paths)$size),
  md5=unname(tools::md5sum(pooled_paths)),
  stringsAsFactors=F,
  check.names=F
)
assert_true(identical(source_manifest$file, current_source_manifest$file), "Pooled source manifest file names changed")
assert_equal_numeric(source_manifest$size_bytes, current_source_manifest$size_bytes, "Pooled source file sizes", tolerance=0)
assert_true(identical(source_manifest$md5, current_source_manifest$md5), "A pooled source file no longer matches its recorded MD5")

DESeq2Obj.original = readRDS(file.path(pooled_output_dir, "DESeq2Obj-after-wald-test.rds"))
DESeq2Obj = readRDS(file.path(guide_output_dir, "DESeq2Obj-guide-specific-after-wald-test.rds"))
assert_true(inherits(DESeq2Obj.original, "DESeqDataSet"), "Pooled fitted object is not a DESeqDataSet")
assert_true(inherits(DESeq2Obj, "DESeqDataSet"), "Guide-specific fitted object is not a DESeqDataSet")
assert_true(NCOL(DESeq2Obj.original) == 44L && NCOL(DESeq2Obj) == 44L, "Pooled and guide-specific fitted objects must retain all 44 samples")
assert_true(
  identical(paste(deparse(design(DESeq2Obj.original)), collapse=" "), "~replicate + comparison_group"),
  "Pooled fitted-object design must remain ~ replicate + comparison_group"
)
assert_true(identical(colnames(DESeq2Obj), colnames(DESeq2Obj.original)), "Guide-specific fitted sample order differs from pooled source")
assert_true(identical(rownames(DESeq2Obj), rownames(DESeq2Obj.original)), "Guide-specific fitted gene universe differs from pooled source")
assert_true(identical(counts(DESeq2Obj), counts(DESeq2Obj.original)), "Guide-specific fitted counts differ from pooled source counts")
assert_true("guide15_rep_4" %in% colnames(DESeq2Obj), "guide15_rep_4 was not retained")
assert_true(all(is.finite(sizeFactors(DESeq2Obj)) & sizeFactors(DESeq2Obj) > 0), "All guide-specific size factors must be finite and positive")
assert_true(all(is.finite(dispersions(DESeq2Obj)) & dispersions(DESeq2Obj) > 0), "All guide-specific dispersions must be finite and positive")
assert_equal_numeric(sizeFactors(DESeq2Obj.original), sizeFactors(DESeq2Obj), "Pooled versus guide-specific size factors")

metadata = as.data.frame(colData(DESeq2Obj), stringsAsFactors=F)
required_metadata = c(
  "count_column", "sample_id", "guide", "replicate", "hz_guide_number",
  "jy_guide_label", "qpcr_oct4_post_treatment", "qpcr_missing",
  "guide_class", "comparison_group", "contrast_group"
)
assert_true(
  all(required_metadata %in% names(metadata)),
  "Guide-specific metadata is missing column(s): %s",
  paste(setdiff(required_metadata, names(metadata)), collapse=", ")
)
assert_true(identical(rownames(metadata), colnames(DESeq2Obj)), "Guide-specific metadata order differs from sample order")
assert_true(identical(as.character(metadata$sample_id), colnames(DESeq2Obj)), "Guide-specific sample IDs differ from column names")
assert_true(setequal(as.character(metadata$guide), expected_guides), "Guide-specific guide set is incorrect")
assert_true(setequal(as.character(metadata$replicate), expected_replicates), "Guide-specific replicate set is incorrect")
assert_true(
  setequal(as.character(metadata$guide[as.character(metadata$comparison_group) == "OCT4_targeting"]), target_guides),
  "Pooled targeting-guide assignment is incorrect"
)
assert_true(
  setequal(as.character(metadata$guide[as.character(metadata$comparison_group) == "pooled_control"]), control_guides),
  "Pooled control-guide assignment is incorrect"
)
guide_replicate_table = table(as.character(metadata$guide), as.character(metadata$replicate))
assert_true(all(guide_replicate_table == 1L), "Every guide must occur exactly once in each replicate block")
assert_true(is.factor(colData(DESeq2Obj)$replicate), "replicate must be a factor")
assert_true(identical(levels(colData(DESeq2Obj)$replicate), expected_replicates), "Replicate factor levels must be 2, 3, 4, and 5")
assert_true(is.factor(colData(DESeq2Obj)$contrast_group), "contrast_group must be a factor")
assert_true(
  identical(levels(colData(DESeq2Obj)$contrast_group), c("pooled_control", target_guides)),
  "contrast_group levels or pooled-control reference are incorrect"
)
assert_true(sum(metadata$contrast_group == "pooled_control") == 12L, "Expected 12 pooled-control libraries")
assert_true(
  setequal(as.character(metadata$guide[metadata$contrast_group == "pooled_control"]), control_guides),
  "Pooled-control guide assignment is incorrect"
)
for (target_guide in target_guides) {
  assert_true(sum(metadata$contrast_group == target_guide) == 4L, "Expected four libraries for %s", target_guide)
  assert_true(
    all(as.character(metadata$guide[metadata$contrast_group == target_guide]) == target_guide),
    "contrast_group assignment is incorrect for %s",
    target_guide
  )
}
assert_true(
  identical(paste(deparse(design(DESeq2Obj)), collapse=" "), "~replicate + contrast_group"),
  "Guide-specific design must be ~ replicate + contrast_group"
)
assert_true(!any(grepl("qpcr", all.vars(design(DESeq2Obj)), ignore.case=T)), "qPCR must not occur in the guide-specific design")
model_matrix = model.matrix(design(DESeq2Obj), data=metadata)
assert_true(NROW(model_matrix) == 44L && NCOL(model_matrix) == 12L, "Guide-specific model matrix must be 44 rows x 12 columns")
assert_true(qr(model_matrix)$rank == NCOL(model_matrix), "Guide-specific model matrix is not full rank")
assert_true(NROW(model_matrix) - qr(model_matrix)$rank == 32L, "Guide-specific model must have 32 residual degrees of freedom")
assert_true(all(expected_coefficients %in% resultsNames(DESeq2Obj)), "One or more guide-specific coefficients are missing")

run_summary = readLines(file.path(guide_output_dir, "run-summary.txt"), warn=F)
alpha_line = grep("^Wald-test alpha: ", run_summary, value=T)
assert_true(length(alpha_line) == 1L, "Guide-specific run summary must contain exactly one Wald-test alpha")
alpha = suppressWarnings(as.numeric(sub("^Wald-test alpha: ", "", alpha_line)))
assert_true(is.finite(alpha) && alpha > 0 && alpha < 1, "Invalid Wald-test alpha in guide-specific run summary")
assert_true(any(run_summary == "Samples retained: 44 / 44"), "Guide-specific run summary does not confirm all 44 samples")
assert_true(any(run_summary == "guide15_rep_4 retained: TRUE"), "Guide-specific run summary does not confirm guide15_rep_4 retention")
assert_true(any(run_summary == "Model matrix: 44 rows x 12 columns; rank 12; residual df 32"), "Guide-specific run summary model dimensions are incorrect")

required_result_columns = c(
  "target_guide", "target_jy_guide_label", "pooled_control_guides", "coefficient",
  "gene_id", "gene_name", "gene_biotype", "baseMean",
  "mean_normalized_count_target_guide", "mean_normalized_count_pooled_control",
  "log2FoldChange_raw", "lfcSE_raw", "log2FoldChange_shrunken", "lfcSE_shrunken",
  "stat", "pvalue", "padj", "is_tested", "has_adjusted_pvalue",
  "independent_filter_flag", "cooks_or_other_filter_flag", "significant_fdr",
  "significant_up", "significant_down", "direction", "result_status"
)
normalized_counts = counts(DESeq2Obj, normalized=T)
control_samples = as.character(colData(DESeq2Obj)$contrast_group) == "pooled_control"
expected_control_means = rowMeans(normalized_counts[, control_samples, drop=F])
result_tables = vector("list", length(target_guides))
names(result_tables) = target_guides
contrast_counts = vector("list", length(target_guides))
names(contrast_counts) = target_guides

for (guide_index in seq_along(target_guides)) {
  target_guide = target_guides[[guide_index]]
  coefficient_name = expected_coefficients[[guide_index]]
  complete_path = file.path(guide_output_dir, sprintf("de-results-%s-vs-pooled-controls-complete.csv", target_guide))
  significant_path = file.path(guide_output_dir, sprintf("de-results-%s-vs-pooled-controls-significant.csv", target_guide))
  complete_results = read.csv(complete_path, check.names=F, na.strings="NA")
  significant_results = read.csv(significant_path, check.names=F, na.strings="NA", colClasses=c(gene_id="character"))

  assert_true(identical(names(complete_results), required_result_columns), "Complete result schema is incorrect for %s", target_guide)
  assert_true(identical(names(significant_results), required_result_columns), "Significant result schema is incorrect for %s", target_guide)
  assert_true(NROW(complete_results) == NROW(DESeq2Obj), "Complete result row count is incorrect for %s", target_guide)
  assert_true(!anyDuplicated(complete_results$gene_id), "Complete result gene IDs are duplicated for %s", target_guide)
  assert_true(setequal(complete_results$gene_id, rownames(DESeq2Obj)), "Complete result genes differ from the fitted object for %s", target_guide)
  assert_true(all(complete_results$target_guide == target_guide), "target_guide values are incorrect for %s", target_guide)
  assert_true(all(complete_results$pooled_control_guides == paste(control_guides, collapse=";")), "Pooled-control guide labels are incorrect for %s", target_guide)
  assert_true(all(complete_results$coefficient == coefficient_name), "Coefficient labels are incorrect for %s", target_guide)
  expected_guide_label = unique(as.character(metadata$jy_guide_label[as.character(metadata$guide) == target_guide]))
  assert_true(length(expected_guide_label) == 1L && all(complete_results$target_jy_guide_label == expected_guide_label), "Workbook guide label is incorrect for %s", target_guide)
  assert_true(
    identical(significant_results$gene_id, complete_results$gene_id[as.logical(complete_results$significant_fdr)]),
    "Significant results are not the exact complete-table subset for %s",
    target_guide
  )

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
  result_match = match(complete_results$gene_id, rownames(diffExprLst))
  assert_true(!anyNA(result_match), "Could not align recomputed DESeq2 results for %s", target_guide)
  assert_equal_numeric(complete_results$baseMean, diffExprLst$baseMean[result_match], sprintf("%s baseMean", target_guide))
  assert_equal_numeric(complete_results$log2FoldChange_raw, diffExprLst$log2FoldChange[result_match], sprintf("%s raw log2 fold change", target_guide))
  assert_equal_numeric(complete_results$lfcSE_raw, diffExprLst$lfcSE[result_match], sprintf("%s raw LFC standard error", target_guide))
  assert_equal_numeric(complete_results$stat, diffExprLst$stat[result_match], sprintf("%s Wald statistic", target_guide))
  assert_equal_numeric(complete_results$pvalue, diffExprLst$pvalue[result_match], sprintf("%s p-value", target_guide))
  assert_equal_numeric(complete_results$padj, diffExprLst$padj[result_match], sprintf("%s adjusted p-value", target_guide))
  assert_equal_numeric(complete_results$log2FoldChange_shrunken, diffExprLst.shrunk$log2FoldChange[result_match], sprintf("%s shrunken log2 fold change", target_guide))
  assert_equal_numeric(complete_results$lfcSE_shrunken, diffExprLst.shrunk$lfcSE[result_match], sprintf("%s shrunken LFC standard error", target_guide))

  expected_tested = !is.na(complete_results$pvalue)
  expected_adjusted = !is.na(complete_results$padj)
  expected_independent_filter = !is.na(complete_results$pvalue) & is.na(complete_results$padj)
  expected_cooks = is.na(complete_results$pvalue) & complete_results$baseMean > 0
  expected_significant = !is.na(complete_results$padj) & complete_results$padj < alpha
  expected_up = expected_significant & !is.na(complete_results$log2FoldChange_shrunken) & complete_results$log2FoldChange_shrunken > 0
  expected_down = expected_significant & !is.na(complete_results$log2FoldChange_shrunken) & complete_results$log2FoldChange_shrunken < 0
  expected_direction = ifelse(expected_up, "higher_in_target_guide", ifelse(expected_down, "lower_in_target_guide", "not_significant"))
  expected_status = ifelse(
    complete_results$baseMean == 0,
    "all_zero",
    ifelse(
      expected_cooks,
      "not_tested_cooks_or_other",
      ifelse(expected_independent_filter, "independent_filter", ifelse(expected_significant, "significant_fdr", "tested_not_significant"))
    )
  )
  assert_true(identical(as.logical(complete_results$is_tested), expected_tested), "is_tested is incorrect for %s", target_guide)
  assert_true(identical(as.logical(complete_results$has_adjusted_pvalue), expected_adjusted), "has_adjusted_pvalue is incorrect for %s", target_guide)
  assert_true(identical(as.logical(complete_results$independent_filter_flag), expected_independent_filter), "independent_filter_flag is incorrect for %s", target_guide)
  assert_true(identical(as.logical(complete_results$cooks_or_other_filter_flag), expected_cooks), "cooks_or_other_filter_flag is incorrect for %s", target_guide)
  assert_true(identical(as.logical(complete_results$significant_fdr), expected_significant), "significant_fdr is incorrect for %s", target_guide)
  assert_true(identical(as.logical(complete_results$significant_up), expected_up), "significant_up is incorrect for %s", target_guide)
  assert_true(identical(as.logical(complete_results$significant_down), expected_down), "significant_down is incorrect for %s", target_guide)
  assert_true(identical(complete_results$direction, expected_direction), "direction is incorrect for %s", target_guide)
  assert_true(identical(complete_results$result_status, expected_status), "result_status is incorrect for %s", target_guide)

  padj_sort = ifelse(is.na(complete_results$padj), Inf, complete_results$padj)
  effect_sort = ifelse(is.na(complete_results$log2FoldChange_shrunken), -Inf, abs(complete_results$log2FoldChange_shrunken))
  expected_order = order(!expected_significant, padj_sort, -effect_sort, complete_results$gene_id)
  assert_true(identical(expected_order, seq_len(NROW(complete_results))), "Complete results are not deterministically sorted for %s", target_guide)

  target_samples = as.character(colData(DESeq2Obj)$contrast_group) == target_guide
  expected_target_means = rowMeans(normalized_counts[, target_samples, drop=F])
  assert_equal_numeric(
    complete_results$mean_normalized_count_target_guide,
    expected_target_means[complete_results$gene_id],
    sprintf("%s target normalized-count means", target_guide)
  )
  assert_equal_numeric(
    complete_results$mean_normalized_count_pooled_control,
    expected_control_means[complete_results$gene_id],
    sprintf("%s pooled-control normalized-count means", target_guide)
  )

  summary_line = sprintf(
    "%s: coefficient=%s; Wald p-value=%d; adjusted p-value=%d; FDR-significant=%d; higher=%d; lower=%d",
    target_guide,
    coefficient_name,
    sum(expected_tested),
    sum(expected_adjusted),
    sum(expected_significant),
    sum(expected_up),
    sum(expected_down)
  )
  assert_true(any(run_summary == summary_line), "Guide-specific run summary totals are incorrect for %s", target_guide)

  result_tables[[target_guide]] = complete_results
  contrast_counts[[target_guide]] = c(
    significant=sum(expected_significant),
    higher=sum(expected_up),
    lower=sum(expected_down)
  )
}

gene_ids = rownames(DESeq2Obj)
result_by_gene = lapply(result_tables, function(result_table) {
  result_match = match(gene_ids, result_table$gene_id)
  assert_true(!anyNA(result_match), "Could not align a result table to the common gene universe")
  result_table[result_match, , drop=F]
})
is_tested_matrix = do.call(cbind, lapply(result_by_gene, function(x) as.logical(x$is_tested)))
has_adjusted_matrix = do.call(cbind, lapply(result_by_gene, function(x) as.logical(x$has_adjusted_pvalue)))
significant_matrix = do.call(cbind, lapply(result_by_gene, function(x) as.logical(x$significant_fdr)))
significant_up_matrix = do.call(cbind, lapply(result_by_gene, function(x) as.logical(x$significant_up)))
significant_down_matrix = do.call(cbind, lapply(result_by_gene, function(x) as.logical(x$significant_down)))
shrunken_lfc_matrix = do.call(cbind, lapply(result_by_gene, function(x) as.numeric(x$log2FoldChange_shrunken)))
colnames(is_tested_matrix) = target_guides
colnames(has_adjusted_matrix) = target_guides
colnames(significant_matrix) = target_guides
colnames(significant_up_matrix) = target_guides
colnames(significant_down_matrix) = target_guides
colnames(shrunken_lfc_matrix) = target_guides
collapse_guides = function(flag_values) paste(target_guides[as.logical(flag_values)], collapse=";")

guides_significant_fdr = rowSums(significant_matrix)
guides_significant_up = rowSums(significant_up_matrix)
guides_significant_down = rowSums(significant_down_matrix)
max_same_direction = pmax(guides_significant_up, guides_significant_down)
consistent_multiple = (guides_significant_up >= 2L & guides_significant_down == 0L) |
  (guides_significant_down >= 2L & guides_significant_up == 0L)
mixed_direction = guides_significant_up > 0L & guides_significant_down > 0L
concordance_status = ifelse(
  mixed_direction,
  "mixed_significant_directions",
  ifelse(
    guides_significant_up >= 2L,
    "consistent_higher_multiple_guides",
    ifelse(
      guides_significant_down >= 2L,
      "consistent_lower_multiple_guides",
      ifelse(guides_significant_fdr == 1L, "significant_one_guide_only", "not_significant_any_guide")
    )
  )
)
expected_concordance = data.frame(
  gene_id=gene_ids,
  gene_name=mcols(DESeq2Obj)[gene_ids, "gene_name"],
  gene_biotype=mcols(DESeq2Obj)[gene_ids, "gene_biotype"],
  guides_expected=length(target_guides),
  guides_with_finite_shrunken_lfc=rowSums(is.finite(shrunken_lfc_matrix)),
  guides_with_wald_pvalue=rowSums(is_tested_matrix),
  guides_with_adjusted_pvalue=rowSums(has_adjusted_matrix),
  guides_significant_fdr=guides_significant_fdr,
  guides_significant_up=guides_significant_up,
  guides_significant_down=guides_significant_down,
  max_same_direction_significant_guides=max_same_direction,
  fraction_guides_significant=guides_significant_fdr / length(target_guides),
  fraction_guides_same_direction=max_same_direction / length(target_guides),
  significant_guides=apply(significant_matrix, 1L, collapse_guides),
  significant_up_guides=apply(significant_up_matrix, 1L, collapse_guides),
  significant_down_guides=apply(significant_down_matrix, 1L, collapse_guides),
  median_log2FoldChange_shrunken=apply(shrunken_lfc_matrix, 1L, finite_summary, summary_function=median),
  min_log2FoldChange_shrunken=apply(shrunken_lfc_matrix, 1L, finite_summary, summary_function=min),
  max_log2FoldChange_shrunken=apply(shrunken_lfc_matrix, 1L, finite_summary, summary_function=max),
  consistent_across_multiple_guides=consistent_multiple,
  mixed_significant_directions=mixed_direction,
  concordance_status=concordance_status,
  stringsAsFactors=F,
  check.names=F
)
concordance_effect_sort = ifelse(
  is.na(expected_concordance$median_log2FoldChange_shrunken),
  -Inf,
  abs(expected_concordance$median_log2FoldChange_shrunken)
)
concordance_order = order(
  !expected_concordance$consistent_across_multiple_guides,
  -expected_concordance$max_same_direction_significant_guides,
  -expected_concordance$guides_significant_fdr,
  -concordance_effect_sort,
  expected_concordance$gene_id
)
expected_concordance = expected_concordance[concordance_order, , drop=F]
rownames(expected_concordance) = NULL

# Explicit text types retain empty guide lists and header-only significance IDs.
concordance = read.csv(
  file.path(guide_output_dir, "concordance-across-guides.csv"),
  check.names=F, na.strings="NA",
  colClasses=c(significant_guides="character", significant_up_guides="character", significant_down_guides="character")
)
assert_true(identical(names(concordance), names(expected_concordance)), "Concordance summary schema is incorrect")
assert_true(NROW(concordance) == NROW(DESeq2Obj), "Concordance summary must contain every retained gene")
assert_true(!anyDuplicated(concordance$gene_id), "Concordance summary gene IDs must be unique")
assert_true(identical(concordance$gene_id, expected_concordance$gene_id), "Concordance summary gene order or membership is incorrect")
assert_true(identical(concordance$gene_name, expected_concordance$gene_name), "Concordance gene names are incorrect")
assert_true(identical(concordance$gene_biotype, expected_concordance$gene_biotype), "Concordance gene biotypes are incorrect")
for (column_name in c(
  "guides_expected", "guides_with_finite_shrunken_lfc", "guides_with_wald_pvalue",
  "guides_with_adjusted_pvalue", "guides_significant_fdr", "guides_significant_up",
  "guides_significant_down", "max_same_direction_significant_guides",
  "fraction_guides_significant", "fraction_guides_same_direction",
  "median_log2FoldChange_shrunken", "min_log2FoldChange_shrunken", "max_log2FoldChange_shrunken"
)) {
  assert_equal_numeric(concordance[[column_name]], expected_concordance[[column_name]], sprintf("Concordance %s", column_name))
}
for (column_name in c("significant_guides", "significant_up_guides", "significant_down_guides", "concordance_status")) {
  assert_true(identical(concordance[[column_name]], expected_concordance[[column_name]]), "Concordance %s is incorrect", column_name)
}
assert_true(
  identical(as.logical(concordance$consistent_across_multiple_guides), expected_concordance$consistent_across_multiple_guides),
  "Concordance consistent_across_multiple_guides is incorrect"
)
assert_true(
  identical(as.logical(concordance$mixed_significant_directions), expected_concordance$mixed_significant_directions),
  "Concordance mixed_significant_directions is incorrect"
)
assert_true(
  any(run_summary == sprintf("Genes consistent across multiple guides: %d", sum(consistent_multiple))),
  "Guide-specific run summary has an incorrect consistent-gene count"
)
assert_true(
  any(run_summary == sprintf("Genes with mixed significant directions: %d", sum(mixed_direction))),
  "Guide-specific run summary has an incorrect mixed-direction count"
)

contrast_count_matrix = do.call(rbind, contrast_counts)
message(sprintf(
  "PASS: verified 44 samples, %d retained genes, eight guide-specific contrasts, and %d genes consistently significant in the same direction across at least two guides.",
  NROW(DESeq2Obj),
  sum(consistent_multiple)
))
message(sprintf(
  "Per-guide FDR-significant totals: %s",
  paste(sprintf("%s=%d", rownames(contrast_count_matrix), contrast_count_matrix[, "significant"]), collapse=", ")
))

devtools::session_info(pkgs="attached")
