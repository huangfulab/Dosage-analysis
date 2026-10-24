#!/usr/bin/env Rscript

# K2LLPS 3'-seq: guide-specific blocked DESeq2 contrasts against pooled controls.

options(stringsAsFactors=F)

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
finite_summary = function(values, summary_function) {
  values = values[is.finite(values)]
  if (length(values) == 0L) return(NA_real_)
  summary_function(values)
}

# Keep staging cleanup in a real function scope, including publication errors.
main = function() {
alpha = suppressWarnings(as.numeric(Sys.getenv("K2LLPS_ALPHA", unset="0.05")))
if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) {
  stop_with("K2LLPS_ALPHA must be a number strictly between 0 and 1")
}

args = commandArgs(trailingOnly=T)
project_dir = script_dir()
pooled_output_dir = if (length(args) >= 1L) args[[1L]] else file.path(project_dir, "outputs")
guide_output_dir = if (length(args) >= 2L) args[[2L]] else file.path(pooled_output_dir, "guide-specific")

if (!requireNamespace("DESeq2", quietly=T)) stop_with("DESeq2 is required but is not installed")
suppressPackageStartupMessages(library(DESeq2))
if (!requireNamespace("devtools", quietly=T)) {
  stop_with("devtools is required to capture the session information but is not installed")
}

pooled_output_dir = normalizePath(pooled_output_dir, winslash="/", mustWork=T)
source_directory = normalizePath(project_dir, winslash="/", mustWork=T)
assert_safe_guide_directory = function(path) {
  path = normalizePath(path, winslash="/", mustWork=F)
  # Publishing replaces this directory; never allow it to contain the source.
  protected = c(pooled_output_dir, source_directory)
  if (any(protected == path | startsWith(protected, paste0(sub("/$", "", path), "/"))) ||
      grepl("^([A-Za-z]:)?/$", path) || grepl("^//[^/]+/[^/]+/?$", path)) {
    stop_with("Guide-specific outputs must not replace or contain the pooled output directory or study source directory")
  }
  path
}
guide_output_dir = assert_safe_guide_directory(guide_output_dir)
# normalizePath cannot resolve nonexistent components (e.g. new/../outputs).
# Once the parent exists, canonicalize it and recheck before doing any fitting.
dir.create(dirname(guide_output_dir), recursive=T, showWarnings=F)
guide_parent = normalizePath(dirname(guide_output_dir), winslash="/", mustWork=T)
guide_output_dir = assert_safe_guide_directory(file.path(guide_parent, basename(guide_output_dir)))

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
pooled_paths = file.path(pooled_output_dir, pooled_files)
missing_pooled_files = pooled_files[!file.exists(pooled_paths)]
if (length(missing_pooled_files) > 0L) {
  stop_with("Missing pooled output file(s): %s", paste(missing_pooled_files, collapse=", "))
}
empty_pooled_files = pooled_files[file.info(pooled_paths)$size <= 0]
if (length(empty_pooled_files) > 0L) {
  stop_with("Empty pooled output file(s): %s", paste(empty_pooled_files, collapse=", "))
}

pooled_manifest = function() {
  file_information = file.info(pooled_paths)
  data.frame(
    file=pooled_files,
    size_bytes=as.numeric(file_information$size),
    md5=unname(tools::md5sum(pooled_paths)),
    stringsAsFactors=F,
    check.names=F
  )
}
assert_manifest_unchanged = function(observed, expected, stage_label) {
  if (!identical(observed$file, expected$file) ||
      !identical(observed$size_bytes, expected$size_bytes) ||
      !identical(observed$md5, expected$md5)) {
    changed = expected$file[
      observed$size_bytes != expected$size_bytes | observed$md5 != expected$md5
    ]
    stop_with(
      "Pooled outputs changed %s: %s",
      stage_label,
      paste(changed, collapse=", ")
    )
  }
}
source_manifest = pooled_manifest()

target_guides = c("guide01", "guide02", "guide03", "guide04", "guide05", "guide06", "guide14", "guide15")
control_guides = c("guide18", "guide19", "guide20")
expected_guides = c(target_guides, control_guides)
expected_replicates = as.character(2:5)
contrast_levels = c("pooled_control", target_guides)

DESeq2Obj.original = readRDS(file.path(pooled_output_dir, "DESeq2Obj-after-wald-test.rds"))
assert_true(inherits(DESeq2Obj.original, "DESeqDataSet"), "Pooled fitted object is not a DESeqDataSet")
assert_true(NCOL(DESeq2Obj.original) == 44L, "Pooled fitted object must contain 44 samples")
assert_true(NROW(DESeq2Obj.original) > 0L, "Pooled fitted object contains no retained genes")
assert_true(!anyDuplicated(rownames(DESeq2Obj.original)), "Pooled fitted object gene IDs must be unique")

metadata = as.data.frame(colData(DESeq2Obj.original), stringsAsFactors=F)
required_metadata = c(
  "count_column", "sample_id", "guide", "replicate", "hz_guide_number",
  "jy_guide_label", "qpcr_oct4_post_treatment", "qpcr_missing",
  "guide_class", "comparison_group"
)
missing_metadata = setdiff(required_metadata, names(metadata))
if (length(missing_metadata) > 0L) {
  stop_with("Pooled fitted-object metadata is missing column(s): %s", paste(missing_metadata, collapse=", "))
}
assert_true(identical(rownames(metadata), colnames(DESeq2Obj.original)), "Pooled fitted-object metadata order differs from sample order")
assert_true(identical(as.character(metadata$sample_id), colnames(DESeq2Obj.original)), "Pooled fitted-object sample IDs differ from column names")
assert_true(setequal(as.character(metadata$guide), expected_guides), "Pooled fitted-object guide set is incorrect")
assert_true(setequal(as.character(metadata$replicate), expected_replicates), "Pooled fitted-object replicate set is incorrect")
guide_replicate_table = table(as.character(metadata$guide), as.character(metadata$replicate))
assert_true(
  identical(sort(rownames(guide_replicate_table)), sort(expected_guides)) &&
    identical(sort(colnames(guide_replicate_table)), sort(expected_replicates)) &&
    all(guide_replicate_table == 1L),
  "Every guide must occur exactly once in each replicate block"
)
assert_true(all(c("gene_name", "gene_biotype") %in% names(mcols(DESeq2Obj.original))), "Pooled fitted object lacks gene annotations")

analysis_counts = counts(DESeq2Obj.original, normalized=F)
assert_true(!anyNA(analysis_counts) && all(is.finite(analysis_counts)) && all(analysis_counts >= 0), "Pooled fitted-object counts are invalid")
fractional_counts = any(abs(analysis_counts - round(analysis_counts)) > 1e-8)
if (fractional_counts && !requireNamespace("glmGamPoi", quietly=T)) {
  stop_with("The pooled fitted object contains fractional counts, but glmGamPoi is not installed")
}

metadata$guide = factor(as.character(metadata$guide), levels=expected_guides)
metadata$replicate = factor(as.character(metadata$replicate), levels=expected_replicates)
metadata$contrast_group = ifelse(as.character(metadata$guide) %in% target_guides, as.character(metadata$guide), "pooled_control")
metadata$contrast_group = factor(metadata$contrast_group, levels=contrast_levels)
rownames(metadata) = metadata$sample_id

assert_true(sum(metadata$contrast_group == "pooled_control") == 12L, "Expected 12 pooled-control libraries")
for (target_guide in target_guides) {
  assert_true(sum(metadata$contrast_group == target_guide) == 4L, "Expected four libraries for %s", target_guide)
}

design_formula = ~ replicate + contrast_group
model_matrix = model.matrix(design_formula, data=metadata)
model_rank = qr(model_matrix)$rank
residual_degrees_freedom = NROW(model_matrix) - model_rank
assert_true(NROW(model_matrix) == 44L && NCOL(model_matrix) == 12L, "Guide-specific model matrix must be 44 rows x 12 columns")
assert_true(model_rank == NCOL(model_matrix), "Guide-specific model matrix is not full rank")
assert_true(residual_degrees_freedom > 0L, "Guide-specific model has no residual degrees of freedom")

if (fractional_counts) {
  se = SummarizedExperiment::SummarizedExperiment(
    assays=list(counts=analysis_counts),
    colData=S4Vectors::DataFrame(metadata)
  )
  DESeq2Obj = DESeqDataSet(se, design=design_formula, skipIntegerMode=T)
  dispersion_fit_type = "glmGamPoi"
} else {
  DESeq2Obj = DESeqDataSetFromMatrix(
    countData=analysis_counts,
    colData=metadata,
    design=design_formula
  )
  dispersion_fit_type = "parametric"
}
mcols(DESeq2Obj)$gene_name = mcols(DESeq2Obj.original)$gene_name
mcols(DESeq2Obj)$gene_biotype = mcols(DESeq2Obj.original)$gene_biotype

DESeq2Obj = estimateSizeFactors(DESeq2Obj)
DESeq2Obj = estimateDispersions(DESeq2Obj, fitType=dispersion_fit_type)
DESeq2Obj = nbinomWaldTest(DESeq2Obj, betaPrior=F)

expected_coefficients = paste0("contrast_group_", target_guides, "_vs_pooled_control")
missing_coefficients = setdiff(expected_coefficients, resultsNames(DESeq2Obj))
if (length(missing_coefficients) > 0L) {
  stop_with("Missing guide-specific coefficient(s): %s", paste(missing_coefficients, collapse=", "))
}

normalized_counts = counts(DESeq2Obj, normalized=T)
control_samples = as.character(colData(DESeq2Obj)$contrast_group) == "pooled_control"
mean_control = rowMeans(normalized_counts[, control_samples, drop=F])
result_tables = vector("list", length(target_guides))
names(result_tables) = target_guides
contrast_summary_parts = vector("list", length(target_guides))
names(contrast_summary_parts) = target_guides

dir.create(dirname(guide_output_dir), recursive=T, showWarnings=F)
staging_dir = tempfile(
  pattern=paste0(".", basename(guide_output_dir), "-staging-"),
  tmpdir=dirname(guide_output_dir)
)
assert_true(dir.create(staging_dir, recursive=F, showWarnings=F), "Could not create guide-specific staging directory")
on.exit(if (dir.exists(staging_dir)) unlink(staging_dir, recursive=T, force=T), add=T)

saveRDS(DESeq2Obj, file.path(staging_dir, "DESeq2Obj-guide-specific-after-wald-test.rds"))

for (guide_index in seq_along(target_guides)) {
  target_guide = target_guides[[guide_index]]
  coefficient_name = expected_coefficients[[guide_index]]
  guide_label = unique(as.character(metadata$jy_guide_label[as.character(metadata$guide) == target_guide]))
  assert_true(length(guide_label) == 1L && !is.na(guide_label) && guide_label != "", "Could not determine workbook label for %s", target_guide)

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

  target_samples = as.character(colData(DESeq2Obj)$contrast_group) == target_guide
  mean_target = rowMeans(normalized_counts[, target_samples, drop=F])
  diffExprLst.df = as.data.frame(diffExprLst)
  diffExprLst.shrunk.df = as.data.frame(diffExprLst.shrunk)
  result_gene_ids = rownames(diffExprLst.df)

  result_table = data.frame(
    target_guide=target_guide,
    target_jy_guide_label=guide_label,
    pooled_control_guides=paste(control_guides, collapse=";"),
    coefficient=coefficient_name,
    gene_id=result_gene_ids,
    gene_name=mcols(DESeq2Obj)[result_gene_ids, "gene_name"],
    gene_biotype=mcols(DESeq2Obj)[result_gene_ids, "gene_biotype"],
    baseMean=diffExprLst.df$baseMean,
    mean_normalized_count_target_guide=mean_target[result_gene_ids],
    mean_normalized_count_pooled_control=mean_control[result_gene_ids],
    log2FoldChange_raw=diffExprLst.df$log2FoldChange,
    lfcSE_raw=diffExprLst.df$lfcSE,
    log2FoldChange_shrunken=diffExprLst.shrunk.df$log2FoldChange,
    lfcSE_shrunken=diffExprLst.shrunk.df$lfcSE,
    stat=diffExprLst.df$stat,
    pvalue=diffExprLst.df$pvalue,
    padj=diffExprLst.df$padj,
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
    "higher_in_target_guide",
    ifelse(result_table$significant_down, "lower_in_target_guide", "not_significant")
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

  complete_name = sprintf("de-results-%s-vs-pooled-controls-complete.csv", target_guide)
  significant_name = sprintf("de-results-%s-vs-pooled-controls-significant.csv", target_guide)
  write.csv(result_table, file.path(staging_dir, complete_name), row.names=F, quote=T, na="NA")
  write.csv(significant_table, file.path(staging_dir, significant_name), row.names=F, quote=T, na="NA")

  result_tables[[target_guide]] = result_table
  contrast_summary_parts[[target_guide]] = data.frame(
    target_guide=target_guide,
    coefficient=coefficient_name,
    libraries=sum(target_samples),
    genes_with_wald_pvalue=sum(result_table$is_tested),
    genes_with_adjusted_pvalue=sum(result_table$has_adjusted_pvalue),
    genes_significant_fdr=sum(result_table$significant_fdr),
    genes_significant_higher=sum(result_table$significant_up),
    genes_significant_lower=sum(result_table$significant_down),
    stringsAsFactors=F,
    check.names=F
  )
}
contrast_summary = do.call(rbind, contrast_summary_parts)
rownames(contrast_summary) = NULL

gene_ids = rownames(DESeq2Obj)
result_by_gene = lapply(result_tables, function(result_table) {
  result_match = match(gene_ids, result_table$gene_id)
  assert_true(!anyNA(result_match), "Could not align a guide result table to the common gene universe")
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

concordance = data.frame(
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
  is.na(concordance$median_log2FoldChange_shrunken),
  -Inf,
  abs(concordance$median_log2FoldChange_shrunken)
)
concordance_order = order(
  !concordance$consistent_across_multiple_guides,
  -concordance$max_same_direction_significant_guides,
  -concordance$guides_significant_fdr,
  -concordance_effect_sort,
  concordance$gene_id
)
concordance = concordance[concordance_order, , drop=F]
rownames(concordance) = NULL
write.csv(concordance, file.path(staging_dir, "concordance-across-guides.csv"), row.names=F, quote=T, na="NA")
write.csv(source_manifest, file.path(staging_dir, "pooled-source-md5.csv"), row.names=F, quote=T, na="NA")

writeLines(capture.output(devtools::session_info(pkgs="attached")), file.path(staging_dir, "sessionInfo.txt"))
summary_lines = c(
  "K2LLPS DESeq2 guide-specific contrast summary",
  sprintf("Pooled source directory: %s", pooled_output_dir),
  sprintf("Guide-specific output directory: %s", guide_output_dir),
  sprintf("Source pooled fitted object MD5: %s", source_manifest$md5[source_manifest$file == "DESeq2Obj-after-wald-test.rds"]),
  sprintf("Design: %s", paste(deparse(design_formula), collapse=" ")),
  "Replicate interpretation: shared experimental block; modeled as a categorical factor",
  sprintf("Model matrix: %d rows x %d columns; rank %d; residual df %d", NROW(model_matrix), NCOL(model_matrix), model_rank, residual_degrees_freedom),
  sprintf("Target guides: %s", paste(target_guides, collapse=", ")),
  sprintf("Pooled control guides: %s", paste(control_guides, collapse=", ")),
  sprintf("Target libraries per guide: %d", 4L),
  sprintf("Pooled control libraries: %d", sum(control_samples)),
  sprintf("Samples retained: %d / 44", NCOL(DESeq2Obj)),
  sprintf("guide15_rep_4 retained: %s", "guide15_rep_4" %in% colnames(DESeq2Obj)),
  sprintf("Genes retained from pooled prefilter: %d", NROW(DESeq2Obj)),
  sprintf("Fractional source counts: %s", fractional_counts),
  sprintf("Dispersion fit: %s", dispersion_fit_type),
  sprintf("Wald-test alpha: %.17g", alpha),
  "Independent filtering: enabled separately for each guide contrast",
  "Cook's-distance filtering: enabled at each results extraction",
  "LFC shrinkage: DESeq2 type='normal' on each named guide-versus-pooled-control coefficient",
  "Coefficient orientation: positive log2 fold change means higher expression in the named target guide",
  "Concordance rule: significant in at least two guides in one direction and no significant guide in the opposite direction",
  sprintf("Genes consistent across multiple guides: %d", sum(concordance$consistent_across_multiple_guides)),
  sprintf("Genes with mixed significant directions: %d", sum(concordance$mixed_significant_directions)),
  "Guide contrast totals:"
)
for (guide_index in seq_len(NROW(contrast_summary))) {
  guide_row = contrast_summary[guide_index, ]
  summary_lines = c(
    summary_lines,
    sprintf(
      "%s: coefficient=%s; Wald p-value=%d; adjusted p-value=%d; FDR-significant=%d; higher=%d; lower=%d",
      guide_row$target_guide,
      guide_row$coefficient,
      guide_row$genes_with_wald_pvalue,
      guide_row$genes_with_adjusted_pvalue,
      guide_row$genes_significant_fdr,
      guide_row$genes_significant_higher,
      guide_row$genes_significant_lower
    )
  )
}
summary_lines = c(
  summary_lines,
  sprintf("Concordance summary: %s", file.path(guide_output_dir, "concordance-across-guides.csv")),
  sprintf("Pooled source manifest: %s", file.path(guide_output_dir, "pooled-source-md5.csv")),
  "Status: guide-specific size factors, dispersions, Wald tests, result extraction, LFC shrinkage, and concordance completed."
)
writeLines(summary_lines, file.path(staging_dir, "run-summary.txt"))

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
staged_paths = file.path(staging_dir, expected_guide_files)
missing_staged_files = expected_guide_files[!file.exists(staged_paths)]
assert_true(length(missing_staged_files) == 0L, "Missing staged guide-specific file(s): %s", paste(missing_staged_files, collapse=", "))
empty_staged_files = expected_guide_files[file.info(staged_paths)$size <= 0]
assert_true(length(empty_staged_files) == 0L, "Empty staged guide-specific file(s): %s", paste(empty_staged_files, collapse=", "))
assert_manifest_unchanged(pooled_manifest(), source_manifest, "while generating guide-specific outputs")

guide_output_dir = assert_safe_guide_directory(guide_output_dir)
backup_dir = NULL
if (dir.exists(guide_output_dir)) {
  backup_dir = paste0(guide_output_dir, ".backup-", Sys.getpid())
  assert_true(!file.exists(backup_dir), "Guide-specific backup path already exists: %s", backup_dir)
  assert_true(file.rename(guide_output_dir, backup_dir), "Could not move the prior guide-specific output directory to a backup")
}
if (!file.rename(staging_dir, guide_output_dir)) {
  if (!is.null(backup_dir) && dir.exists(backup_dir)) file.rename(backup_dir, guide_output_dir)
  stop_with("Could not publish the staged guide-specific output directory")
}
if (!is.null(backup_dir) && dir.exists(backup_dir)) {
  unlink(backup_dir, recursive=T, force=T)
  assert_true(!dir.exists(backup_dir), "Published outputs, but could not remove the prior guide-specific backup")
}
assert_manifest_unchanged(pooled_manifest(), source_manifest, "after publishing guide-specific outputs")

message("Done. Guide-specific differential-expression testing completed.")
for (guide_index in seq_len(NROW(contrast_summary))) {
  guide_row = contrast_summary[guide_index, ]
  message(sprintf(
    "%s FDR < %.3g: %d genes (%d higher, %d lower)",
    guide_row$target_guide,
    alpha,
    guide_row$genes_significant_fdr,
    guide_row$genes_significant_higher,
    guide_row$genes_significant_lower
  ))
}
message(sprintf("Consistent across multiple guides: %d genes", sum(concordance$consistent_across_multiple_guides)))
message("Guide-specific outputs: ", guide_output_dir)
}

main()
devtools::session_info(pkgs="attached")
