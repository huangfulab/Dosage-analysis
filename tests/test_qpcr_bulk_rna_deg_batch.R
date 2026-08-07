#!/usr/bin/env Rscript

script_argument = grep("^--file=", commandArgs(trailingOnly = F), value = T)
script_path = normalizePath(sub("^--file=", "", script_argument), mustWork = T)
repository_directory = normalizePath(file.path(dirname(script_path), ".."), mustWork = T)
source(file.path(repository_directory, "R", "qpcr_bulk_rna_models.R"))

# Use common samples with linear, quadratic, and LL.4-like bulk-RNA fixtures
# to test batch alignment, model summaries, and percentage-point disagreement.
set.seed(20260728)
grna = rep(sprintf("gRNA_%02d", seq_len(11L)), each = 4L)
qpcr_center = rep(c(28, 38, 52, 66, 82, 100, 125, 158, 200, 252, 318), each = 4L)
qpcr_percent = qpcr_center + stats::rnorm(44L, mean = 0, sd = 2)
sample_data = data.frame(
  sample_id = sprintf("sample_%02d", seq_len(44L)),
  grna = grna,
  qpcr_percent = qpcr_percent,
  stringsAsFactors = F
)
long_data = do.call(rbind, list(
  data.frame(sample_data, gene = "GENE_LINEAR", bulk_rna_percent = 20 + 0.7 * qpcr_percent + stats::rnorm(44L, sd = 7)),
  data.frame(sample_data, gene = "GENE_QUADRATIC", bulk_rna_percent = 20 + 0.1 * qpcr_percent + 0.002 * qpcr_percent^2 + stats::rnorm(44L, sd = 7)),
  data.frame(sample_data, gene = "GENE_LL4", bulk_rna_percent = 30 + 175 / (1 + exp(-3.8 * (log(qpcr_percent) - log(118)))) + stats::rnorm(44L, sd = 4))
))

batch_fit = fit_qpcr_bulk_rna_deg_batch(
  long_data,
  sample_id_col = "sample_id",
  grna_col = "grna",
  qpcr_col = "qpcr_percent",
  gene_col = "gene",
  bulk_rna_col = "bulk_rna_percent"
)
stopifnot(
  inherits(batch_fit, "qpcr_bulk_rna_deg_batch"),
  length(batch_fit$genes) == 3L,
  all(vapply(batch_fit$aggregates, function(aggregate) NROW(aggregate$point_data) == 11L, logical(1))),
  all(vapply(batch_fit$aggregates, function(aggregate) all(aggregate$point_data$n_samples == 4L), logical(1)))
)

model_aic = collect_qpcr_bulk_rna_deg_results(batch_fit, "model_table")
stopifnot(
  NROW(model_aic) == 9L,
  all(table(model_aic$gene) == 3L),
  all(model_aic$model %in% c("linear", "quadratic", "ll4"))
)

concordance = qpcr_bulk_rna_deg_concordance(batch_fit)
stopifnot(
  identical(dim(concordance$matrix), c(3L, 11L)),
  NROW(concordance$long) == 33L,
  all.equal(concordance$long$difference_percentage_points,
    concordance$long$bulk_rna_percent - concordance$long$qpcr_percent,
    check.attributes = FALSE) == TRUE,
  concordance$limits[1L] == -concordance$limits[2L],
  max(abs(concordance$matrix), na.rm = TRUE) <= concordance$limits[2L]
)
model_summary = summarize_qpcr_bulk_rna_deg_models(batch_fit)
stopifnot(
  NROW(model_summary) == 9L,
  all(table(model_summary$gene) == 3L),
  all(model_summary$convergence_status == "converged"),
  all(vapply(split(model_summary, model_summary$gene), function(rows) sum(rows$aic_winner) >= 1L, logical(1)))
)
# Force a numerical AIC tie and a failed candidate to verify deterministic
# winner/tie labels rather than arbitrary selection.
tied_fit = batch_fit
tied_fit$fits[[1L]]$model_table$aic[1:2] = 10
tied_fit$fits[[1L]]$model_table$converged[3L] = FALSE
tied_fit$fits[[1L]]$model_table$convergence_message[3L] = "forced failure"
tied_summary = summarize_qpcr_bulk_rna_deg_models(tied_fit)
stopifnot(
  all(tied_summary$winner_label[tied_summary$gene == batch_fit$genes[1L] & tied_summary$model %in% c("linear", "quadratic")] == "tied_winner"),
  tied_summary$winner_label[tied_summary$gene == batch_fit$genes[1L] & tied_summary$model == "ll4"] == "failed"
)
heatmap_caller_file = tempfile(fileext = ".pdf")
grDevices::pdf(heatmap_caller_file)
heatmap_par_values = c("mfrow", "mar", "oma", "cex")
graphics::par(mfrow = c(1L, 1L), mar = c(4.7, 4.2, 3.3, 2.4), oma = c(0, 0, 0, 0), cex = 0.93)
heatmap_par_before = graphics::par(heatmap_par_values)
heatmap_devices_before = as.integer(grDevices::dev.list())
heatmap_file = tempfile(fileext = ".pdf")
heatmap_result = plot_qpcr_bulk_rna_deg_concordance_heatmap(batch_fit, heatmap_file)
stopifnot(
  file.exists(heatmap_file), file.info(heatmap_file)$size > 0,
  identical(as.integer(grDevices::dev.list()), heatmap_devices_before),
  isTRUE(all.equal(graphics::par(heatmap_par_values), heatmap_par_before))
)
unlink(heatmap_file)
grDevices::dev.off()
unlink(heatmap_caller_file)
overview_file = tempfile(fileext = ".pdf")
overview_result = plot_qpcr_bulk_rna_deg_model_overview(batch_fit, overview_file)
stopifnot(overview_result$pages == 3L, file.exists(overview_file), file.info(overview_file)$size > 0)
unlink(overview_file)

plot_file = tempfile(fileext = ".pdf")
plot_result = plot_qpcr_bulk_rna_deg_batch(batch_fit, file = plot_file)
stopifnot(
  file.exists(plot_file),
  file.info(plot_file)$size > 0,
  plot_result$pages == 3L
)
unlink(plot_file)
diagnostic_file = tempfile(fileext = ".pdf")
diagnostic_result = plot_qpcr_bulk_rna_deg_batch_diagnostics(batch_fit, file = diagnostic_file)
stopifnot(
  file.exists(diagnostic_file),
  file.info(diagnostic_file)$size > 0,
  diagnostic_result$pages == 3L
)
pdftotext = Sys.which("pdftotext")
if (nzchar(pdftotext)) {
  text_file = tempfile(fileext = ".txt")
  status = system2(pdftotext, c(shQuote(diagnostic_file), shQuote(text_file)))
  stopifnot(status == 0L)
  pdf_text = paste(readLines(text_file, warn = FALSE), collapse = "\n")
  stopifnot(all(vapply(batch_fit$genes, grepl, logical(1), x = pdf_text, fixed = TRUE)))
  unlink(text_file)
}
unlink(diagnostic_file)

# Missing a complete gRNA and changing shared qPCR both violate intentional
# batch-design invariants and therefore must stop rather than be auto-aligned.
missing_point_data = long_data[!(long_data$gene == "GENE_LL4" & long_data$grna == "gRNA_11"), , drop = F]
point_count_error = tryCatch(
  {
    fit_qpcr_bulk_rna_deg_batch(
      missing_point_data,
      sample_id_col = "sample_id",
      grna_col = "grna",
      qpcr_col = "qpcr_percent",
      gene_col = "gene",
      bulk_rna_col = "bulk_rna_percent"
    )
    F
  },
  error = function(error) {
    grepl("Expected 11 gRNA-level points", conditionMessage(error), fixed = T)
  }
)
stopifnot(point_count_error)

inconsistent_qpcr_data = long_data
inconsistent_qpcr_data$qpcr_percent[
  inconsistent_qpcr_data$gene == "GENE_LL4" &
    inconsistent_qpcr_data$sample_id == "sample_01"
] = 999
inconsistent_qpcr_error = tryCatch(
  {
    fit_qpcr_bulk_rna_deg_batch(
      inconsistent_qpcr_data,
      sample_id_col = "sample_id",
      grna_col = "grna",
      qpcr_col = "qpcr_percent",
      gene_col = "gene",
      bulk_rna_col = "bulk_rna_percent"
    )
    F
  },
  error = function(error) grepl("inconsistent gRNA or qPCR", conditionMessage(error), fixed = T)
)
stopifnot(inconsistent_qpcr_error)

cat("All batch qPCR--bulk RNA-seq model tests passed.\n")
