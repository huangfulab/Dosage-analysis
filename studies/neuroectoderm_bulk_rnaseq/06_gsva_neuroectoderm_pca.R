#!/usr/bin/env Rscript

# K2LLPS 3'-seq: GSVA / ssGSEA neuroectoderm scores on the 44-sample VST matrix,
# plotted on the existing sample PCA coordinates.
#
# Usage:
#   Rscript 06_gsva_neuroectoderm_pca.R [output-directory]
#
# Gene-set provenance is the gene_set_provenance object below. The 24 symbols
# were downloaded from Enrichr CellMarker_2024 on 2026-09-08, not from a local
# Jielin xlsx. Named "neuroectoderm" alternatives were either 4 genes (too small
# for ssGSEA) or GO:0060563 (35 genes, missing SOX1/PAX6/NES).

options(stringsAsFactors=F)

script_dir = function() {
  file_arg = grep("^--file=", commandArgs(trailingOnly=F), value=T)
  if (length(file_arg) == 0L) return(normalizePath(getwd(), mustWork=T))
  script_file = sub("^--file=", "", file_arg[[1L]])
  return(dirname(normalizePath(script_file, mustWork=T)))
}
stop_with = function(...) stop(sprintf(...), call.=F)

args = commandArgs(trailingOnly=T)
project_dir = script_dir()
output_dir = if (length(args) >= 1L) args[[1L]] else file.path(project_dir, "outputs")

if (!requireNamespace("DESeq2", quietly=T)) stop_with("DESeq2 is required but is not installed")
if (!requireNamespace("GSVA", quietly=T)) stop_with("GSVA is required but is not installed")
if (!requireNamespace("ggplot2", quietly=T)) stop_with("ggplot2 is required but is not installed")
if (!requireNamespace("devtools", quietly=T)) {
  stop_with("devtools is required to capture the session information but is not installed")
}

suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(GSVA))

output_dir = normalizePath(output_dir, mustWork=T)

# Downloaded gene-set provenance. Symbols are the Enrichr CellMarker_2024 line
# "Neural Stem Cell Brain Human" in the original tab-separated order.
gene_set_provenance = list(
  retrieved_date="2026-09-08",
  retrieved_by="K2LLPS_DESeq2/06_gsva_neuroectoderm_pca.R",
  download_url="https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=CellMarker_2024",
  download_method="HTTP GET, Enrichr geneSetLibrary mode=text",
  library_name="CellMarker_2024",
  geneset_name="Neural Stem Cell Brain Human",
  n_symbols=24L,
  identifier_type="HGNC gene symbol",
  enrichr_record="Neural Stem Cell Brain Human<TAB><empty description><TAB>ADGRG1 ... VIM",
  primary_database="CellMarker",
  primary_database_url="http://bio-bigdata.hrbmu.edu.cn/CellMarker/",
  citations=c(
    "Zhang et al. 2019 Nucleic Acids Res 47:D721-D728; PMID 30289549; doi:10.1093/nar/gky900 (CellMarker)",
    "Hu et al. 2023 Nucleic Acids Res 51:D870-D876; PMID 36300619; doi:10.1093/nar/gkac947 (CellMarker 2.0)",
    "Kuleshov et al. 2016 Nucleic Acids Res 44:W90-W97; PMID 27141961; doi:10.1093/nar/gkw377 (Enrichr)"
  ),
  selection_note=paste(
    "Compact curated neural-lineage marker set including SOX1, SOX2, PAX6, NES, and HES5.",
    "Chosen over CellMarker Neuroectoderm (4 genes) and MSigDB GOBP_NEUROEPITHELIAL_CELL_DIFFERENTIATION / GO:0060563 (35 genes, missing SOX1/PAX6/NES)."
  )
)
ne_gene_symbols = c(
  "ADGRG1", "ASCL1", "BMI1", "CD24", "CD44", "CD81", "CXCR4", "DCX",
  "HES5", "ID1", "ITGA6", "ITGB1", "MAPT", "MCAM", "NES", "NGFR",
  "NOTCH1", "PAX6", "PROM1", "SLITRK6", "SOX1", "SOX2", "SOX9", "VIM"
)
if (
  length(ne_gene_symbols) != gene_set_provenance$n_symbols ||
  anyDuplicated(ne_gene_symbols) ||
  length(ne_gene_symbols) >= 100L
) {
  stop_with(
    "Embedded gene list must match provenance n_symbols=%d, be unique, and stay under 100",
    gene_set_provenance$n_symbols
  )
}

required_files = c(
  "DESeq2Obj-after-wald-test.rds",
  "analysis-sample-metadata.csv",
  "pca-coordinates.csv"
)
missing_files = required_files[!file.exists(file.path(output_dir, required_files))]
if (length(missing_files) > 0L) {
  stop_with("Missing pooled output file(s): %s", paste(missing_files, collapse=", "))
}

message("Loading fitted DESeq2 object and PCA coordinates")
DESeq2Obj = readRDS(file.path(output_dir, "DESeq2Obj-after-wald-test.rds"))
metadata = read.csv(
  file.path(output_dir, "analysis-sample-metadata.csv"),
  stringsAsFactors=F,
  check.names=F
)
pca_coordinates = read.csv(
  file.path(output_dir, "pca-coordinates.csv"),
  stringsAsFactors=F,
  check.names=F
)

if (NCOL(DESeq2Obj) != 44L) stop_with("Fitted object must contain 44 samples, found %d", NCOL(DESeq2Obj))
if (NROW(metadata) != 44L) stop_with("analysis-sample-metadata.csv must contain 44 rows, found %d", NROW(metadata))
if (NROW(pca_coordinates) != 44L) stop_with("pca-coordinates.csv must contain 44 rows, found %d", NROW(pca_coordinates))
if (!setequal(colnames(DESeq2Obj), metadata$sample_id)) {
  stop_with("Fitted object sample IDs do not match analysis-sample-metadata.csv")
}
if (!setequal(pca_coordinates$sample_id, metadata$sample_id)) {
  stop_with("PCA coordinate sample IDs do not match analysis-sample-metadata.csv")
}

annotation = data.frame(
  gene_id = rownames(DESeq2Obj),
  gene_name = as.character(mcols(DESeq2Obj)$gene_name),
  gene_biotype = as.character(mcols(DESeq2Obj)$gene_biotype),
  baseMean = unname(rowMeans(counts(DESeq2Obj, normalized=T))),
  stringsAsFactors=F,
  check.names=F
)
annotation$gene_name[is.na(annotation$gene_name)] = ""

map_symbol = function(symbol) {
  hits = annotation[annotation$gene_name == symbol, , drop=F]
  if (NROW(hits) == 0L) {
    return(data.frame(
      gene_symbol=symbol,
      mapped=F,
      gene_id=NA_character_,
      gene_biotype=NA_character_,
      n_ids=0L,
      stringsAsFactors=F,
      check.names=F
    ))
  }
  coding_hits = hits[hits$gene_biotype == "protein_coding", , drop=F]
  chosen = if (NROW(coding_hits) > 0L) coding_hits else hits
  chosen = chosen[order(-chosen$baseMean, chosen$gene_id), , drop=F]
  data.frame(
    gene_symbol=symbol,
    mapped=T,
    gene_id=chosen$gene_id[[1L]],
    gene_biotype=chosen$gene_biotype[[1L]],
    n_ids=NROW(hits),
    stringsAsFactors=F,
    check.names=F
  )
}

gene_set_map = do.call(rbind, lapply(ne_gene_symbols, map_symbol))
rownames(gene_set_map) = NULL
n_mapped = sum(gene_set_map$mapped)
unmapped_symbols = gene_set_map$gene_symbol[!gene_set_map$mapped]
if (n_mapped < 5L) {
  stop_with("Only %d / %d gene-set symbols mapped into the fitted object", n_mapped, length(ne_gene_symbols))
}
if (anyDuplicated(gene_set_map$gene_id[gene_set_map$mapped])) {
  stop_with("Mapped gene-set Ensembl IDs are not unique")
}

mapped_ids = gene_set_map$gene_id[gene_set_map$mapped]
gene_sets = list(NEUROECTODERM_NSC=mapped_ids)
message(sprintf("Mapped %d / %d gene-set symbols into the %d-gene VST matrix", n_mapped, length(ne_gene_symbols), NROW(DESeq2Obj)))
if (length(unmapped_symbols) > 0L) {
  message("Unmapped symbols: ", paste(unmapped_symbols, collapse=", "))
}

message("Computing variance-stabilizing transformation on the full fitted matrix")
vst_data = varianceStabilizingTransformation(DESeq2Obj, blind=F)
vst_matrix = assay(vst_data)
if (!identical(colnames(vst_matrix), colnames(DESeq2Obj))) {
  stop_with("VST sample order does not match the fitted object")
}
pca_check = plotPCA(vst_data, intgroup="comparison_group", returnData=T)
percent_variance = round(100 * attr(pca_check, "percentVar"), 1)
pca_check = pca_check[match(pca_coordinates$sample_id, rownames(pca_check)), , drop=F]
if (anyNA(rownames(pca_check))) stop_with("plotPCA reconstruction lost one or more sample IDs")
if (
  !isTRUE(all.equal(unname(pca_check$PC1), pca_coordinates$PC1, tolerance=1e-6, check.attributes=F)) ||
  !isTRUE(all.equal(unname(pca_check$PC2), pca_coordinates$PC2, tolerance=1e-6, check.attributes=F))
) {
  stop_with("Reconstructed PCA coordinates do not match pca-coordinates.csv")
}

message("Running ssGSEA")
ssgsea_param = ssgseaParam(
  exprData=vst_matrix,
  geneSets=gene_sets,
  minSize=1,
  maxSize=Inf,
  alpha=0.25,
  normalize=T,
  verbose=T
)
ssgsea_scores = gsva(ssgsea_param)
if (!is.matrix(ssgsea_scores)) ssgsea_scores = assay(ssgsea_scores)
if (!"NEUROECTODERM_NSC" %in% rownames(ssgsea_scores)) {
  stop_with("ssGSEA result is missing the NEUROECTODERM_NSC gene set")
}

message("Running GSVA")
gsva_param = gsvaParam(
  exprData=vst_matrix,
  geneSets=gene_sets,
  minSize=1,
  maxSize=Inf,
  kcdf="Gaussian",
  verbose=T
)
gsva_scores = gsva(gsva_param)
if (!is.matrix(gsva_scores)) gsva_scores = assay(gsva_scores)
if (!"NEUROECTODERM_NSC" %in% rownames(gsva_scores)) {
  stop_with("GSVA result is missing the NEUROECTODERM_NSC gene set")
}

sample_scores = pca_coordinates
sample_scores$ssgsea = unname(ssgsea_scores["NEUROECTODERM_NSC", sample_scores$sample_id])
sample_scores$gsva = unname(gsva_scores["NEUROECTODERM_NSC", sample_scores$sample_id])
if (anyNA(sample_scores$ssgsea) || anyNA(sample_scores$gsva)) {
  stop_with("Enrichment scores contain missing values")
}
if (NROW(sample_scores) != 44L) stop_with("Sample score table must contain 44 rows")
# Barbie normalization rescales ssGSEA scores but does not make them nonnegative.

plot_surface = "#fcfcfb"
plot_text = "#0b0b0b"
plot_muted = "#898781"
plot_grid = "#e1e0d9"
plot_axis = "#c3c2b7"

make_pca_score_plot = function(plot_data, score_column, method_label, color_subtitle, fill_mode) {
  plot_data$score = plot_data[[score_column]]
  if (!all(is.finite(plot_data$score))) {
    stop_with("%s enrichment scores contain non-finite values", method_label)
  }
  plot_data$comparison_group = factor(
    plot_data$comparison_group,
    levels=c("pooled_control", "OCT4_targeting")
  )
  if (identical(fill_mode, "signed")) {
    score_limit = max(abs(plot_data$score))
    if (score_limit == 0) score_limit = 1
    fill_scale = ggplot2::scale_fill_gradient2(
      low="#2166ac",
      mid="#ffffff",
      high="#b2182b",
      midpoint=0,
      limits=c(-score_limit, score_limit),
      name=sprintf("%s score", method_label)
    )
  } else if (identical(fill_mode, "range")) {
    score_min = min(plot_data$score)
    score_max = max(plot_data$score)
    # A constant score is valid; pad only the plotting limits, never the scores.
    if (score_min == score_max) {
      padding = max(0.01, abs(score_min) * 0.01)
      score_min = score_min - padding
      score_max = score_max + padding
    }
    fill_scale = ggplot2::scale_fill_gradientn(
      colours=c("grey50", "yellow", "red"),
      limits=c(score_min, score_max),
      name=sprintf("%s score", method_label)
    )
  } else {
    stop_with("fill_mode must be 'signed' or 'range'")
  }
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x=PC1, y=PC2, fill=score, shape=comparison_group)
  ) +
    ggplot2::geom_hline(yintercept=0, color=plot_axis, linewidth=0.45) +
    ggplot2::geom_vline(xintercept=0, color=plot_axis, linewidth=0.45) +
    ggplot2::geom_point(size=3.2, color=plot_text, stroke=0.35) +
    fill_scale +
    ggplot2::scale_shape_manual(
      values=c(pooled_control=21, OCT4_targeting=24),
      labels=c(pooled_control="Pooled control", OCT4_targeting="OCT4 targeting"),
      name=NULL
    ) +
    ggplot2::guides(
      fill=ggplot2::guide_colorbar(order=1, barwidth=8.5, barheight=0.6, title.position="top"),
      shape=ggplot2::guide_legend(order=2, override.aes=list(fill=plot_muted, color=plot_text))
    ) +
    ggplot2::labs(
      title=sprintf("K2LLPS PCA colored by %s neuroectoderm score", method_label),
      subtitle=color_subtitle,
      x=sprintf("PC1: %.1f%% variance", percent_variance[[1L]]),
      y=sprintf("PC2: %.1f%% variance", percent_variance[[2L]])
    ) +
    ggplot2::coord_cartesian(clip="off") +
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
      plot.title=ggplot2::element_text(color=plot_text, size=13, face="plain", hjust=0.5),
      plot.subtitle=ggplot2::element_text(color=plot_muted, size=8.5, hjust=0.5),
      legend.position="top",
      legend.box="vertical",
      legend.justification="left",
      legend.background=ggplot2::element_blank(),
      legend.key=ggplot2::element_blank(),
      legend.text=ggplot2::element_text(color=plot_text, size=10),
      legend.title=ggplot2::element_text(color=plot_text, size=10),
      plot.margin=ggplot2::margin(12, 18, 12, 12)
    )
}

ssgsea_plot = make_pca_score_plot(
  sample_scores,
  "ssgsea",
  "ssGSEA",
  sprintf(
    "Existing VST / top-500 PCA; ssGSEA on the full %d-gene matrix; fill = raw score, grey50-yellow-red over min-max; %d / %d symbols mapped",
    NROW(vst_matrix),
    n_mapped,
    length(ne_gene_symbols)
  ),
  "range"
)
gsva_plot = make_pca_score_plot(
  sample_scores,
  "gsva",
  "GSVA",
  sprintf(
    "Existing VST / top-500 PCA; GSVA on the full %d-gene matrix; fill = signed GSVA score; %d / %d symbols mapped",
    NROW(vst_matrix),
    n_mapped,
    length(ne_gene_symbols)
  ),
  "signed"
)

ggplot2::ggsave(
  filename=file.path(output_dir, "gsva-neuroectoderm-ssgsea-pca.png"),
  plot=ssgsea_plot,
  width=10,
  height=8,
  units="in",
  dpi=300,
  bg=plot_surface
)
ggplot2::ggsave(
  filename=file.path(output_dir, "gsva-neuroectoderm-gsva-pca.png"),
  plot=gsva_plot,
  width=10,
  height=8,
  units="in",
  dpi=300,
  bg=plot_surface
)

combined_path = file.path(output_dir, "gsva-neuroectoderm-ssgsea-gsva-pca.png")
grDevices::png(combined_path, width=18, height=8, units="in", res=300, bg=plot_surface)
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout=grid::grid.layout(1L, 2L)))
print(ssgsea_plot, vp=grid::viewport(layout.pos.row=1L, layout.pos.col=1L))
print(gsva_plot, vp=grid::viewport(layout.pos.row=1L, layout.pos.col=2L))
grDevices::dev.off()

gene_set_map$source_library = gene_set_provenance$library_name
gene_set_map$source_geneset = gene_set_provenance$geneset_name
gene_set_map$source_url = gene_set_provenance$download_url
gene_set_map$retrieved_date = gene_set_provenance$retrieved_date
write.csv(
  gene_set_map,
  file.path(output_dir, "gsva-neuroectoderm-gene-set.csv"),
  row.names=F,
  quote=F,
  na="NA"
)
provenance_lines = c(
  "K2LLPS neuroectoderm gene-set download provenance",
  sprintf("Retrieved date: %s", gene_set_provenance$retrieved_date),
  sprintf("Retrieved by: %s", gene_set_provenance$retrieved_by),
  sprintf("Download URL: %s", gene_set_provenance$download_url),
  sprintf("Download method: %s", gene_set_provenance$download_method),
  sprintf("Enrichr library: %s", gene_set_provenance$library_name),
  sprintf("Gene-set name: %s", gene_set_provenance$geneset_name),
  sprintf("Symbols: %d %s", gene_set_provenance$n_symbols, gene_set_provenance$identifier_type),
  sprintf("Enrichr record: %s", gene_set_provenance$enrichr_record),
  sprintf("Primary database: %s", gene_set_provenance$primary_database),
  sprintf("Primary database URL: %s", gene_set_provenance$primary_database_url),
  "Citations:",
  paste0("  - ", gene_set_provenance$citations),
  sprintf("Selection note: %s", gene_set_provenance$selection_note),
  sprintf("Embedded symbols: %s", paste(ne_gene_symbols, collapse=", "))
)
writeLines(provenance_lines, file.path(output_dir, "gsva-neuroectoderm-gene-set-provenance.txt"))
write.csv(
  sample_scores,
  file.path(output_dir, "gsva-neuroectoderm-sample-scores.csv"),
  row.names=F,
  quote=F,
  na="NA"
)

summary_lines = c(
  "K2LLPS GSVA / ssGSEA neuroectoderm summary",
  sprintf("Samples: %d", NROW(sample_scores)),
  sprintf("Expression matrix: VST of fitted DESeq2 object (%d genes x %d samples)", NROW(vst_matrix), NCOL(vst_matrix)),
  "PCA coordinates: existing outputs/pca-coordinates.csv (VST, top 500 variable genes)",
  sprintf("PC1/PC2 percent variance: %.1f / %.1f", percent_variance[[1L]], percent_variance[[2L]]),
  sprintf("Gene set: Enrichr %s / %s", gene_set_provenance$library_name, gene_set_provenance$geneset_name),
  sprintf("Gene-set download URL: %s", gene_set_provenance$download_url),
  sprintf("Gene-set retrieved: %s", gene_set_provenance$retrieved_date),
  sprintf("Gene-set citations: %s", paste(gene_set_provenance$citations, collapse=" | ")),
  sprintf("Gene-set symbols: %d", length(ne_gene_symbols)),
  sprintf("Mapped symbols: %d", n_mapped),
  sprintf("Unmapped symbols: %s", if (length(unmapped_symbols) == 0L) "none" else paste(unmapped_symbols, collapse=", ")),
  "Algorithms: ssgseaParam(normalize=TRUE, alpha=0.25); gsvaParam(kcdf=Gaussian)",
  sprintf("ssGSEA score range: %.4f to %.4f", min(sample_scores$ssgsea), max(sample_scores$ssgsea)),
  sprintf("GSVA score range: %.4f to %.4f", min(sample_scores$gsva), max(sample_scores$gsva)),
  "ssGSEA fill: raw score, grey50-yellow-red, limits min(score) to max(score)",
  "GSVA fill: signed score, blue-white-red, midpoint 0, symmetric limits at max(|score|)"
)
writeLines(summary_lines, file.path(output_dir, "gsva-neuroectoderm-run-summary.txt"))
writeLines(
  capture.output(devtools::session_info(pkgs="attached")),
  file.path(output_dir, "gsva-neuroectoderm-sessionInfo.txt")
)

message("Wrote:")
message("  ", file.path(output_dir, "gsva-neuroectoderm-gene-set.csv"))
message("  ", file.path(output_dir, "gsva-neuroectoderm-gene-set-provenance.txt"))
message("  ", file.path(output_dir, "gsva-neuroectoderm-sample-scores.csv"))
message("  ", file.path(output_dir, "gsva-neuroectoderm-ssgsea-pca.png"))
message("  ", file.path(output_dir, "gsva-neuroectoderm-gsva-pca.png"))
message("  ", file.path(output_dir, "gsva-neuroectoderm-ssgsea-gsva-pca.png"))
message(sprintf("ssGSEA range: %.4f to %.4f", min(sample_scores$ssgsea), max(sample_scores$ssgsea)))
message(sprintf("GSVA range: %.4f to %.4f", min(sample_scores$gsva), max(sample_scores$gsva)))
