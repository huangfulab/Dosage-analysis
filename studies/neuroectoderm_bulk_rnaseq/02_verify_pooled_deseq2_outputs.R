#!/usr/bin/env Rscript

# Fail-fast verification for the K2LLPS blocked pooled-control DESeq2 analysis.

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

args = commandArgs(trailingOnly=T)
project_dir = script_dir()
output_dir = if (length(args) >= 1L) args[[1L]] else file.path(project_dir, "outputs")
# The workbook is an external input; see the study README for placeholders.
design_workbook_path = if (length(args) >= 2L) args[[2L]] else file.path(project_dir, "inputs", "260320_OCT4_qPCR.xlsx")

if (!requireNamespace("DESeq2", quietly=T)) stop_with("DESeq2 is required but is not installed")
if (!requireNamespace("openxlsx", quietly=T)) stop_with("openxlsx is required but is not installed")
suppressPackageStartupMessages(library(DESeq2))

output_dir = normalizePath(output_dir, mustWork=T)
design_workbook_path = normalizePath(design_workbook_path, mustWork=T)

expected_files = c(
  "analysis-sample-metadata.csv",
  "sample_qc.csv",
  "DESeq2Obj-after-wald-test.rds",
  "de-results-oct4-targeting-vs-pooled-controls-complete.csv",
  "de-results-oct4-targeting-vs-pooled-controls-significant.csv",
  "dispersion-estimates.png",
  "expression-header-map.csv",
  "pca-coordinates.csv",
  "pca-oct4-targeting-vs-pooled-controls.png",
  "pca-loadings-oct4-targeting-vs-pooled-controls.csv",
  "pca-pc1-loadings-top-bottom-30-oct4-targeting-vs-pooled-controls.csv",
  "pca-pc1-loadings-top-bottom-30-oct4-targeting-vs-pooled-controls.png",
  "pca-pc1-eigenvector-arrows-oct4-targeting-vs-pooled-controls.csv",
  "pca-pc1-eigenvector-arrows-oct4-targeting-vs-pooled-controls.png",
  "pca-pc2-loadings-top-bottom-30-oct4-targeting-vs-pooled-controls.csv",
  "pca-pc2-loadings-top-bottom-30-oct4-targeting-vs-pooled-controls.png",
  "pca-pc2-eigenvector-arrows-oct4-targeting-vs-pooled-controls.csv",
  "pca-pc2-eigenvector-arrows-oct4-targeting-vs-pooled-controls.png",
  "pca-pc1-score-by-grna-species.csv",
  "pca-pc1-score-by-grna-species-medians.csv",
  "pca-pc1-score-violin-by-grna-species.png",
  "ma-oct4-targeting-vs-pooled-controls.png",
  "degs-oct4-targeting-vs-pooled-controls-lfc1-padj0.05.csv",
  "deg-heatmap-vst-rowz-oct4-targeting-vs-pooled-controls.png",
  "deg-violin-data-oct4-targeting-vs-pooled-controls.csv",
  "deg-violin-oct4-targeting-vs-pooled-controls.png",
  "run-summary.txt",
  "sessionInfo.txt"
)
expected_paths = file.path(output_dir, expected_files)
missing_files = expected_files[!file.exists(expected_paths)]
assert_true(length(missing_files) == 0L, "Missing output file(s): %s", paste(missing_files, collapse=", "))
empty_files = expected_files[file.info(expected_paths)$size <= 0]
assert_true(length(empty_files) == 0L, "Empty output file(s): %s", paste(empty_files, collapse=", "))

equivalence_report_path = file.path(output_dir, "expression-file-equivalence.txt")
assert_true(file.exists(equivalence_report_path), "Missing expression-file-equivalence.txt; run 00_verify_expression_files.R")
equivalence_report = readLines(equivalence_report_path, warn=F)
assert_true(
  any(trimws(equivalence_report) == "PRACTICALLY_EQUAL: TRUE"),
  "Expression-file equivalence report does not contain PRACTICALLY_EQUAL: TRUE"
)

metadata = read.csv(file.path(output_dir, "analysis-sample-metadata.csv"), check.names=F, na.strings="NA")
required_metadata = c(
  "count_column", "sample_id", "guide", "replicate", "hz_guide_number",
  "jy_guide_label", "qpcr_oct4_post_treatment", "qpcr_missing",
  "guide_class", "comparison_group"
)
assert_true(
  all(required_metadata %in% names(metadata)),
  "Analysis metadata is missing required columns: %s",
  paste(setdiff(required_metadata, names(metadata)), collapse=", ")
)
assert_true(NROW(metadata) == 44L, "Analysis metadata must contain 44 rows")
assert_true(!anyDuplicated(metadata$count_column), "count_column values must be unique")
assert_true(!anyDuplicated(metadata$sample_id), "sample_id values must be unique")
metadata$key = paste(metadata$guide, metadata$replicate, sep="::")
assert_true(!anyDuplicated(metadata$key), "Guide/replicate keys must be unique")
assert_true(
  identical(metadata$sample_id, paste0(metadata$guide, "_rep_", metadata$replicate)),
  "sample_id must equal '<guide>_rep_<replicate>'"
)

# Independently recover the biological sample IDs from the verified CSV-to-TSV
# header map. This catches a valid-looking but permuted count-column crosswalk.
header_map = read.csv(file.path(output_dir, "expression-header-map.csv"), check.names=F)
required_header_map = c("csv_header", "tsv_header")
assert_true(
  all(required_header_map %in% names(header_map)),
  "Expression header map is missing required columns: %s",
  paste(setdiff(required_header_map, names(header_map)), collapse=", ")
)
count_header_map = header_map[
  grepl("_count$", header_map$csv_header) & grepl("_count$", header_map$tsv_header),
  required_header_map,
  drop=F
]
assert_true(NROW(count_header_map) == 44L, "Expression header map must contain 44 count-column mappings")
assert_true(
  !anyNA(count_header_map$csv_header) && !anyNA(count_header_map$tsv_header),
  "Expression header map contains missing count headers"
)
assert_true(
  !anyDuplicated(count_header_map$csv_header) && !anyDuplicated(count_header_map$tsv_header),
  "Count headers in the expression header map must be unique"
)
assert_true(
  setequal(metadata$count_column, count_header_map$tsv_header),
  "Analysis metadata count columns differ from the verified expression header map"
)
count_header_map = count_header_map[match(metadata$count_column, count_header_map$tsv_header), , drop=F]
mapped_sample_ids = sub("_count$", "", count_header_map$csv_header)
assert_true(
  identical(metadata$sample_id, mapped_sample_ids),
  "Analysis metadata sample IDs do not match the verified CSV-to-TSV header map"
)

expected_target_guides = c("guide01", "guide02", "guide03", "guide04", "guide05", "guide06", "guide14", "guide15")
expected_control_guides = c("guide18", "guide19", "guide20")
expected_guides = c(expected_target_guides, expected_control_guides)
assert_true(setequal(unique(metadata$guide), expected_guides), "Guide membership does not match the expected 11 guides")
assert_true(setequal(unique(metadata$replicate), 2:5), "Replicate values must be exactly 2, 3, 4, and 5")
guide_replicate_table = table(metadata$guide, metadata$replicate)
assert_true(all(guide_replicate_table == 1L), "Every guide must occur exactly once in each replicate block")
assert_true(sum(metadata$comparison_group == "OCT4_targeting") == 32L, "Expected 32 OCT4-targeting libraries")
assert_true(sum(metadata$comparison_group == "pooled_control") == 12L, "Expected 12 pooled-control libraries")
assert_true(
  setequal(unique(metadata$guide[metadata$comparison_group == "OCT4_targeting"]), expected_target_guides),
  "OCT4-targeting guide assignment is incorrect"
)
assert_true(
  setequal(unique(metadata$guide[metadata$comparison_group == "pooled_control"]), expected_control_guides),
  "Pooled-control guide assignment is incorrect"
)
assert_true(
  setequal(unique(metadata$guide[metadata$guide_class == "nontargeting_control"]), c("guide18", "guide19")),
  "Nontargeting-control assignment is incorrect"
)
assert_true(
  identical(unique(metadata$guide[metadata$guide_class == "safe_harbor_control"]), "guide20"),
  "Safe-harbor-control assignment is incorrect"
)
assert_true(
  identical(metadata$sample_id[is.na(metadata$qpcr_oct4_post_treatment)], "guide04_rep_5"),
  "The only missing qPCR value must be guide04_rep_5"
)
assert_true(
  identical(metadata$sample_id[as.logical(metadata$qpcr_missing)], "guide04_rep_5"),
  "qpcr_missing must flag only guide04_rep_5"
)

# Rebuild the workbook keys independently and compare the linked labels and qPCR values.
design_sheet = openxlsx::read.xlsx(
  design_workbook_path,
  sheet="Sheet1",
  check.names=F,
  na.strings=c("NA", "")
)
required_design_columns = c(
  "HZ_gRNA_name", "JY_gRNA_name",
  "rep2_OCT4", "rep3_OCT4", "rep4_OCT4", "rep5_OCT4"
)
assert_true(all(required_design_columns %in% names(design_sheet)), "Workbook schema has changed")
workbook_guide = sprintf("guide%02d", as.integer(design_sheet$HZ_gRNA_name))
workbook_long = do.call(rbind, lapply(2:5, function(replicate_value) {
  data.frame(
    guide=workbook_guide,
    replicate=replicate_value,
    hz_guide_number=as.integer(design_sheet$HZ_gRNA_name),
    jy_guide_label=trimws(as.character(design_sheet$JY_gRNA_name)),
    qpcr_oct4_post_treatment=suppressWarnings(as.numeric(design_sheet[[paste0("rep", replicate_value, "_OCT4")]])),
    stringsAsFactors=F
  )
}))
workbook_long$key = paste(workbook_long$guide, workbook_long$replicate, sep="::")
assert_true(NROW(workbook_long) == 44L && !anyDuplicated(workbook_long$key), "Workbook must yield 44 unique guide/replicate keys")
workbook_match = match(metadata$key, workbook_long$key)
assert_true(!anyNA(workbook_match) && !anyDuplicated(workbook_match), "Workbook and analysis metadata do not match one-to-one")
assert_true(
  identical(metadata$jy_guide_label, workbook_long$jy_guide_label[workbook_match]),
  "JY guide labels do not match the workbook"
)
assert_true(
  identical(as.integer(metadata$hz_guide_number), workbook_long$hz_guide_number[workbook_match]),
  "HZ guide numbers do not match the workbook"
)
assert_equal_numeric(
  metadata$qpcr_oct4_post_treatment,
  workbook_long$qpcr_oct4_post_treatment[workbook_match],
  "qPCR values"
)

DESeq2Obj = readRDS(file.path(output_dir, "DESeq2Obj-after-wald-test.rds"))
assert_true(inherits(DESeq2Obj, "DESeqDataSet"), "DESeq2Obj-after-wald-test.rds is not a DESeqDataSet")
assert_true(NCOL(DESeq2Obj) == 44L, "The fitted DESeqDataSet must retain all 44 samples")
assert_true(identical(colnames(DESeq2Obj), metadata$sample_id), "Fitted DESeqDataSet sample order differs from analysis metadata")
assert_true("guide15_rep_4" %in% colnames(DESeq2Obj), "guide15_rep_4 was not retained")
assert_true(all(is.finite(sizeFactors(DESeq2Obj)) & sizeFactors(DESeq2Obj) > 0), "All size factors must be finite and positive")
assert_true(all(is.finite(dispersions(DESeq2Obj)) & dispersions(DESeq2Obj) > 0), "All dispersions must be finite and positive")
assert_true(is.factor(colData(DESeq2Obj)$replicate), "replicate must be a factor in colData")
assert_true(identical(levels(colData(DESeq2Obj)$replicate), as.character(2:5)), "Replicate factor levels must be 2, 3, 4, and 5")
assert_true(is.factor(colData(DESeq2Obj)$comparison_group), "comparison_group must be a factor in colData")
assert_true(
  identical(levels(colData(DESeq2Obj)$comparison_group), c("pooled_control", "OCT4_targeting")),
  "pooled_control must be the comparison_group reference level"
)
assert_true(
  identical(paste(deparse(design(DESeq2Obj)), collapse=" "), "~replicate + comparison_group"),
  "Fitted design must be ~ replicate + comparison_group"
)
assert_true(!any(grepl("qpcr", all.vars(design(DESeq2Obj)), ignore.case=T)), "qPCR must not occur in the fitted design")
model_matrix = model.matrix(design(DESeq2Obj), data=as.data.frame(colData(DESeq2Obj)))
assert_true(qr(model_matrix)$rank == NCOL(model_matrix), "Fitted model matrix is not full rank")
coefficient_name = "comparison_group_OCT4_targeting_vs_pooled_control"
assert_true(coefficient_name %in% resultsNames(DESeq2Obj), "Targeting-minus-control coefficient is missing")

run_summary = readLines(file.path(output_dir, "run-summary.txt"), warn=F)
alpha_line = grep("^Wald-test alpha: ", run_summary, value=T)
assert_true(length(alpha_line) == 1L, "run-summary.txt must contain exactly one Wald-test alpha")
alpha = suppressWarnings(as.numeric(sub("^Wald-test alpha: ", "", alpha_line)))
assert_true(is.finite(alpha) && alpha > 0 && alpha < 1, "Invalid Wald-test alpha in run-summary.txt")
assert_true(any(grepl(paste0("^Coefficient: ", coefficient_name, "$"), run_summary)), "Run summary coefficient does not match the fitted object")
assert_true(any(run_summary == "Samples retained: 44 / 44"), "Run summary does not confirm retention of all 44 samples")

# Sensitivity check for a replicate-specific targeting effect. The likelihood-
# ratio test compares the interaction model with the fitted additive design.
interaction_design = ~ replicate + comparison_group + replicate:comparison_group
interaction_model_matrix = model.matrix(
  interaction_design,
  data=as.data.frame(colData(DESeq2Obj))
)
assert_true(
  qr(interaction_model_matrix)$rank == NCOL(interaction_model_matrix),
  "Replicate-by-comparison-group interaction model matrix is not full rank"
)
assert_true(
  NROW(interaction_model_matrix) - NCOL(interaction_model_matrix) == 36L,
  "Interaction model must have 36 residual degrees of freedom"
)
# Preserve the pooled count strategy when rebuilding the sensitivity model.
count_strategy_line = grep("^Count strategy: ", run_summary, value=T)
assert_true(length(count_strategy_line) == 1L, "Missing pooled count strategy in run-summary.txt")
count_strategy = sub("^Count strategy: ", "", count_strategy_line)
assert_true(count_strategy %in% c("round", "glmGamPoi"), "Unknown pooled count strategy")
if (identical(count_strategy, "glmGamPoi")) {
  assert_true(requireNamespace("glmGamPoi", quietly=T), "glmGamPoi is required for fractional-count verification")
  interaction_se = SummarizedExperiment::SummarizedExperiment(
    assays=list(counts=counts(DESeq2Obj)),
    colData=colData(DESeq2Obj)
  )
  interaction_model = DESeqDataSet(interaction_se, design=interaction_design, skipIntegerMode=T)
} else {
  interaction_model = DESeqDataSetFromMatrix(
    countData=counts(DESeq2Obj),
    colData=as.data.frame(colData(DESeq2Obj)),
    design=interaction_design
  )
}
interaction_model = DESeq(
  interaction_model,
  fitType=if (identical(count_strategy, "glmGamPoi")) "glmGamPoi" else "parametric",
  test="LRT",
  reduced=~ replicate + comparison_group,
  minReplicatesForReplace=Inf,
  quiet=T
)
interaction_alpha = 0.05
interaction_results = results(
  interaction_model,
  alpha=interaction_alpha,
  independentFiltering=T,
  cooksCutoff=T
)
interaction_significant = !is.na(interaction_results$padj) & interaction_results$padj < interaction_alpha
interaction_significant_gene_ids = rownames(interaction_results)[interaction_significant]
expected_interaction_gene_ids = "ENSG00000143614"
assert_true(
  identical(interaction_significant_gene_ids, expected_interaction_gene_ids),
  "Expected only GATAD2B (%s) to have an FDR < 0.05 replicate-by-group interaction; found [%s]",
  expected_interaction_gene_ids,
  paste(interaction_significant_gene_ids, collapse=", ")
)

complete_results = read.csv(
  file.path(output_dir, "de-results-oct4-targeting-vs-pooled-controls-complete.csv"),
  check.names=F,
  na.strings="NA"
)
significant_results = read.csv(
  file.path(output_dir, "de-results-oct4-targeting-vs-pooled-controls-significant.csv"),
  check.names=F,
  na.strings="NA",
  colClasses=c(gene_id="character")
)
required_result_columns = c(
  "gene_id", "gene_name", "gene_biotype", "baseMean",
  "mean_normalized_count_oct4_targeting", "mean_normalized_count_pooled_control",
  "log2FoldChange_raw", "lfcSE_raw", "log2FoldChange_shrunken", "lfcSE_shrunken",
  "stat", "pvalue", "padj", "is_tested", "has_adjusted_pvalue",
  "independent_filter_flag", "cooks_or_other_filter_flag", "significant_fdr",
  "significant_up", "significant_down", "direction", "result_status"
)
assert_true(
  all(required_result_columns %in% names(complete_results)),
  "Complete results are missing required columns: %s",
  paste(setdiff(required_result_columns, names(complete_results)), collapse=", ")
)
assert_true(NROW(complete_results) == NROW(DESeq2Obj), "Complete result row count must equal the retained gene count")
assert_true(!anyDuplicated(complete_results$gene_id), "Complete result gene IDs must be unique")
assert_true(setequal(complete_results$gene_id, rownames(DESeq2Obj)), "Complete result genes differ from the fitted DESeqDataSet")
assert_true(
  identical(significant_results$gene_id, complete_results$gene_id[as.logical(complete_results$significant_fdr)]),
  "Significant results must exactly equal the significant subset of complete results"
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
assert_true(!anyNA(result_match), "Could not align complete results to recomputed DESeq2 results")
assert_equal_numeric(complete_results$baseMean, diffExprLst$baseMean[result_match], "baseMean")
assert_equal_numeric(complete_results$log2FoldChange_raw, diffExprLst$log2FoldChange[result_match], "raw log2 fold change")
assert_equal_numeric(complete_results$lfcSE_raw, diffExprLst$lfcSE[result_match], "raw LFC standard error")
assert_equal_numeric(complete_results$stat, diffExprLst$stat[result_match], "Wald statistic")
assert_equal_numeric(complete_results$pvalue, diffExprLst$pvalue[result_match], "p-value")
assert_equal_numeric(complete_results$padj, diffExprLst$padj[result_match], "adjusted p-value")
assert_equal_numeric(
  complete_results$log2FoldChange_shrunken,
  diffExprLst.shrunk$log2FoldChange[result_match],
  "shrunken log2 fold change"
)
assert_equal_numeric(
  complete_results$lfcSE_shrunken,
  diffExprLst.shrunk$lfcSE[result_match],
  "shrunken LFC standard error"
)

expected_significant = !is.na(complete_results$padj) & complete_results$padj < alpha
expected_up = expected_significant & !is.na(complete_results$log2FoldChange_shrunken) & complete_results$log2FoldChange_shrunken > 0
expected_down = expected_significant & !is.na(complete_results$log2FoldChange_shrunken) & complete_results$log2FoldChange_shrunken < 0
assert_true(identical(as.logical(complete_results$significant_fdr), expected_significant), "significant_fdr does not equal nonmissing padj < alpha")
assert_true(identical(as.logical(complete_results$significant_up), expected_up), "significant_up is inconsistent")
assert_true(identical(as.logical(complete_results$significant_down), expected_down), "significant_down is inconsistent")
assert_true(!any(as.logical(complete_results$significant_up) & as.logical(complete_results$significant_down)), "Up/down flags must be mutually exclusive")
assert_true(
  all(complete_results$direction[expected_up] == "higher_in_OCT4_targeting") &&
    all(complete_results$direction[expected_down] == "lower_in_OCT4_targeting"),
  "Result directions do not match the targeting-minus-control coefficient orientation"
)

padj_sort = ifelse(is.na(complete_results$padj), Inf, complete_results$padj)
effect_sort = ifelse(is.na(complete_results$log2FoldChange_shrunken), -Inf, abs(complete_results$log2FoldChange_shrunken))
expected_order = order(!expected_significant, padj_sort, -effect_sort, complete_results$gene_id)
assert_true(identical(expected_order, seq_len(NROW(complete_results))), "Complete results are not deterministically sorted")

normalized_counts = counts(DESeq2Obj, normalized=T)
targeting_samples = colData(DESeq2Obj)$comparison_group == "OCT4_targeting"
control_samples = colData(DESeq2Obj)$comparison_group == "pooled_control"
expected_targeting_means = rowMeans(normalized_counts[, targeting_samples, drop=F])
expected_control_means = rowMeans(normalized_counts[, control_samples, drop=F])
assert_equal_numeric(
  complete_results$mean_normalized_count_oct4_targeting,
  expected_targeting_means[complete_results$gene_id],
  "OCT4-targeting normalized-count means"
)
assert_equal_numeric(
  complete_results$mean_normalized_count_pooled_control,
  expected_control_means[complete_results$gene_id],
  "Pooled-control normalized-count means"
)

sample_qc = read.csv(file.path(output_dir, "sample_qc.csv"), check.names=F, na.strings="NA")
assert_true(NROW(sample_qc) == 44L, "sample_qc.csv must contain 44 rows")
assert_true(identical(sample_qc$sample_id, metadata$sample_id), "Sample QC order differs from analysis metadata")
assert_true(all(as.logical(sample_qc$included_in_analysis)), "Every sample must be marked included_in_analysis")
assert_true("guide15_rep_4" %in% sample_qc$sample_id, "guide15_rep_4 is missing from sample QC")

pca_coordinates = read.csv(file.path(output_dir, "pca-coordinates.csv"), check.names=F, na.strings="NA")
assert_true(NROW(pca_coordinates) == 44L, "PCA coordinates must contain 44 rows")
assert_true(!anyDuplicated(pca_coordinates$sample_id), "PCA sample IDs must be unique")
assert_true(setequal(pca_coordinates$sample_id, metadata$sample_id), "PCA samples differ from analysis metadata")
assert_true(all(is.finite(pca_coordinates$PC1) & is.finite(pca_coordinates$PC2)), "PCA coordinates must be finite")

pca_loadings = read.csv(
  file.path(output_dir, "pca-loadings-oct4-targeting-vs-pooled-controls.csv"),
  check.names=F,
  na.strings="NA"
)
required_pca_loading_columns = c(
  "gene_id", "gene_name", "gene_biotype", "row_variance",
  "PC1_loading", "PC2_loading", "abs_PC1_loading", "abs_PC2_loading",
  "PC_loading_radius", "PC1_direction", "PC2_direction"
)
assert_true(
  all(required_pca_loading_columns %in% names(pca_loadings)),
  "PCA loadings table is missing required columns: %s",
  paste(setdiff(required_pca_loading_columns, names(pca_loadings)), collapse=", ")
)
assert_true(NROW(pca_loadings) == 500L, "PCA loadings table must contain the top 500 variable genes")
assert_true(!anyDuplicated(pca_loadings$gene_id), "PCA loading gene IDs must be unique")
assert_true(all(is.finite(pca_loadings$PC1_loading) & is.finite(pca_loadings$PC2_loading)), "PCA loadings must be finite")
assert_true(
  isTRUE(all.equal(pca_loadings$PC1_loading, sort(pca_loadings$PC1_loading, decreasing=T), check.attributes=F)),
  "PCA loadings table must be sorted by decreasing PC1 loading"
)
assert_equal_numeric(
  pca_loadings$PC_loading_radius,
  sqrt(pca_loadings$PC1_loading ^ 2 + pca_loadings$PC2_loading ^ 2),
  "PCA loading radii"
)

# Reconstruct the same VST PCA used for the sample plot and confirm loadings/scores.
vst_data = varianceStabilizingTransformation(DESeq2Obj, blind=F)
vst_matrix = assay(vst_data)
row_variances = matrixStats::rowVars(vst_matrix)
names(row_variances) = rownames(vst_matrix)
finite_variance_genes = names(row_variances)[is.finite(row_variances)]
pca_selected_genes = names(sort(row_variances[finite_variance_genes], decreasing=T))[seq_len(500L)]
assert_true(setequal(pca_loadings$gene_id, pca_selected_genes), "PCA loadings gene set differs from top-500 variable genes")
pca_fit = prcomp(t(vst_matrix[pca_selected_genes, , drop=F]), center=T, scale.=F)
pca_score_check = data.frame(
  sample_id = rownames(pca_fit$x),
  PC1 = pca_fit$x[, 1L],
  PC2 = pca_fit$x[, 2L],
  stringsAsFactors=F,
  check.names=F
)
pca_score_check = pca_score_check[match(pca_coordinates$sample_id, pca_score_check$sample_id), , drop=F]
assert_equal_numeric(pca_score_check$PC1, pca_coordinates$PC1, "Reconstructed PCA PC1 scores")
assert_equal_numeric(pca_score_check$PC2, pca_coordinates$PC2, "Reconstructed PCA PC2 scores")
loading_match = match(pca_loadings$gene_id, rownames(pca_fit$rotation))
assert_true(!anyNA(loading_match), "PCA loadings table contains genes absent from the reconstructed rotation matrix")
assert_equal_numeric(
  pca_loadings$PC1_loading,
  pca_fit$rotation[loading_match, 1L],
  "Exported PC1 loadings"
)
assert_equal_numeric(
  pca_loadings$PC2_loading,
  pca_fit$rotation[loading_match, 2L],
  "Exported PC2 loadings"
)

required_pca_arrow_columns = c(
  "gene_id", "gene_name", "PC1_loading", "PC2_loading", "PC_loading_radius",
  "PC1_arrow", "PC2_arrow", "display_label"
)
pca_axis_checks = list(
  list(
    axis="PC1",
    loading_column="PC1_loading",
    extremes_csv="pca-pc1-loadings-top-bottom-30-oct4-targeting-vs-pooled-controls.csv",
    arrows_csv="pca-pc1-eigenvector-arrows-oct4-targeting-vs-pooled-controls.csv"
  ),
  list(
    axis="PC2",
    loading_column="PC2_loading",
    extremes_csv="pca-pc2-loadings-top-bottom-30-oct4-targeting-vs-pooled-controls.csv",
    arrows_csv="pca-pc2-eigenvector-arrows-oct4-targeting-vs-pooled-controls.csv"
  )
)
for (axis_check in pca_axis_checks) {
  pca_loadings_extremes = read.csv(
    file.path(output_dir, axis_check$extremes_csv),
    check.names=F,
    na.strings="NA"
  )
  assert_true(
    NROW(pca_loadings_extremes) == 60L,
    "%s top/bottom loadings table must contain 60 genes",
    axis_check$axis
  )
  assert_true(
    !anyDuplicated(pca_loadings_extremes$gene_id),
    "%s extreme loading gene IDs must be unique",
    axis_check$axis
  )
  assert_true(
    isTRUE(all.equal(
      pca_loadings_extremes[[axis_check$loading_column]],
      sort(pca_loadings_extremes[[axis_check$loading_column]], decreasing=F),
      check.attributes=F
    )),
    "%s top/bottom loadings must be sorted by increasing %s loading",
    axis_check$axis,
    axis_check$axis
  )

  axis_sorted_ids = pca_loadings$gene_id[order(-pca_loadings[[axis_check$loading_column]], pca_loadings$gene_id)]
  expected_extreme_ids = c(tail(axis_sorted_ids, 30L), head(axis_sorted_ids, 30L))
  assert_true(
    identical(sort(pca_loadings_extremes$gene_id), sort(expected_extreme_ids)),
    "%s top/bottom loadings must be the first and last 30 genes on that axis",
    axis_check$axis
  )

  extreme_match = match(pca_loadings_extremes$gene_id, pca_loadings$gene_id)
  assert_true(
    !anyNA(extreme_match),
    "%s extreme loadings contain genes absent from the full loadings table",
    axis_check$axis
  )
  assert_equal_numeric(
    pca_loadings_extremes$PC1_loading,
    pca_loadings$PC1_loading[extreme_match],
    sprintf("%s extreme PC1 loading values", axis_check$axis)
  )
  assert_equal_numeric(
    pca_loadings_extremes$PC2_loading,
    pca_loadings$PC2_loading[extreme_match],
    sprintf("%s extreme PC2 loading values", axis_check$axis)
  )

  pca_arrows = read.csv(
    file.path(output_dir, axis_check$arrows_csv),
    check.names=F,
    na.strings="NA"
  )
  assert_true(
    all(required_pca_arrow_columns %in% names(pca_arrows)),
    "%s eigenvector arrow table is missing required columns: %s",
    axis_check$axis,
    paste(setdiff(required_pca_arrow_columns, names(pca_arrows)), collapse=", ")
  )
  assert_true(NROW(pca_arrows) == 60L, "%s eigenvector arrow table must contain 60 genes", axis_check$axis)
  assert_true(!anyDuplicated(pca_arrows$gene_id), "%s eigenvector arrow gene IDs must be unique", axis_check$axis)
  assert_true(
    all(pca_arrows$gene_id %in% pca_loadings$gene_id),
    "%s arrow genes must come from the PCA loadings table",
    axis_check$axis
  )
  assert_true(
    identical(sort(pca_arrows$gene_id), sort(pca_loadings_extremes$gene_id)),
    "%s eigenvector arrows must use the same top/bottom gene set as the bar plot",
    axis_check$axis
  )
  expected_arrow_radius = sqrt(pca_arrows$PC1_loading ^ 2 + pca_arrows$PC2_loading ^ 2)
  assert_equal_numeric(
    pca_arrows$PC_loading_radius,
    expected_arrow_radius,
    sprintf("%s arrow loading radii", axis_check$axis)
  )
  extreme_for_arrows = pca_loadings_extremes[match(pca_arrows$gene_id, pca_loadings_extremes$gene_id), , drop=F]
  assert_equal_numeric(
    pca_arrows$PC1_loading,
    extreme_for_arrows$PC1_loading,
    sprintf("%s arrow PC1 loadings", axis_check$axis)
  )
  assert_equal_numeric(
    pca_arrows$PC2_loading,
    extreme_for_arrows$PC2_loading,
    sprintf("%s arrow PC2 loadings", axis_check$axis)
  )
  # Arrow coordinates are native loadings; no sample-plane rescaling is applied.
  assert_equal_numeric(
    pca_arrows$PC1_arrow,
    pca_arrows$PC1_loading,
    sprintf("%s native PC1 arrow coordinates", axis_check$axis)
  )
  assert_equal_numeric(
    pca_arrows$PC2_arrow,
    pca_arrows$PC2_loading,
    sprintf("%s native PC2 arrow coordinates", axis_check$axis)
  )
}

# Use the threshold recorded by the generating run, not the verifier environment.
deg_lfc_line = grep("^Strict DEG absolute raw LFC cutoff: ", run_summary, value=T)
assert_true(length(deg_lfc_line) == 1L, "Missing strict DEG LFC cutoff; rerun 01_prepare_deseq2.R")
deg_lfc_threshold = suppressWarnings(as.numeric(sub("^Strict DEG absolute raw LFC cutoff: ", "", deg_lfc_line)))
assert_true(is.finite(deg_lfc_threshold) && deg_lfc_threshold >= 0, "Invalid strict DEG LFC cutoff")
deg_table = read.csv(
  file.path(output_dir, "degs-oct4-targeting-vs-pooled-controls-lfc1-padj0.05.csv"),
  check.names=F,
  na.strings="NA"
)
expected_deg = (
  !is.na(complete_results$padj) &
  complete_results$padj <= alpha &
  !is.na(complete_results$log2FoldChange_raw) &
  abs(complete_results$log2FoldChange_raw) >= deg_lfc_threshold
)
expected_deg_ids = complete_results$gene_id[expected_deg]
assert_true(setequal(deg_table$gene_id, expected_deg_ids), "Strict DEG gene set mismatch")
assert_true(NROW(deg_table) == sum(expected_deg), "Strict DEG row count mismatch")
assert_true(all(abs(deg_table$log2FoldChange_raw) >= deg_lfc_threshold), "DEG table contains |raw LFC| below the recorded cutoff")
assert_true(all(!is.na(deg_table$padj) & deg_table$padj <= alpha), "DEG table contains padj > alpha")
expected_deg_up = sum(expected_deg & complete_results$log2FoldChange_raw > 0)
expected_deg_down = sum(expected_deg & complete_results$log2FoldChange_raw < 0)
assert_true(
  sum(deg_table$log2FoldChange_raw > 0) == expected_deg_up,
  "Strict DEG up-count mismatch"
)
assert_true(
  sum(deg_table$log2FoldChange_raw < 0) == expected_deg_down,
  "Strict DEG down-count mismatch"
)

violin_data = read.csv(
  file.path(output_dir, "deg-violin-data-oct4-targeting-vs-pooled-controls.csv"),
  check.names=F,
  na.strings="NA"
)
assert_true(NROW(violin_data) > 0L, "DEG violin data is empty")
assert_true(
  all(c("gene_id", "sample_id", "row_z", "violin_group", "direction") %in% names(violin_data)),
  "DEG violin data is missing required columns"
)
assert_true(all(is.finite(violin_data$row_z)), "DEG violin row_z values must be finite")
violin_groups = sort(unique(as.character(violin_data$violin_group)))
expected_violin_groups = sort(c(
  "pooled_control", "guide01", "guide02", "guide03", "guide04",
  "guide05", "guide06", "guide14", "guide15"
))
assert_true(identical(violin_groups, expected_violin_groups), "DEG violin groups mismatch")
assert_true(length(unique(violin_data$direction)) == 2L, "DEG violin plot must contain up and down panels")
# 9 groups x 2 directions = 18 violins.
assert_true(
  NROW(unique(violin_data[, c("direction", "violin_group")])) == 18L,
  "DEG violin data must define exactly 18 direction-by-group violins"
)

message(sprintf(
  paste0(
    "PASS: verified 44 samples, %d retained genes, and %d FDR-significant genes ",
    "(%d higher, %d lower in OCT4 targeting); strict DEGs (recorded raw LFC cutoff, padj<=alpha): %d ",
    "(up %d, down %d). Interaction sanity check: only GATAD2B (%s) is FDR-significant."
  ),
  NROW(complete_results),
  sum(expected_significant),
  sum(expected_up),
  sum(expected_down),
  sum(expected_deg),
  expected_deg_up,
  expected_deg_down,
  expected_interaction_gene_ids
))

devtools::session_info(pkgs="attached")
