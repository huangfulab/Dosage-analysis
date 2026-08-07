#!/usr/bin/env Rscript

# Generate deterministic synthetic test data, not measured DEG expression or
# evidence for any gene-specific biological response.
arguments = commandArgs(trailingOnly = T)
output_file = if (length(arguments) == 0L) {
  file.path("data", "example_qpcr_bulk_rna_degs.csv")
} else if (length(arguments) == 1L) {
  arguments[[1L]]
} else {
  stop("Usage: Rscript scripts/generate_qpcr_bulk_rna_deg_example.R [output.csv]", call. = F)
}

# Repeat the same 44 samples for three deliberately different linear,
# quadratic, and LL.4-like gene fixtures to test batch alignment and ranking.
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

expression_by_gene = list(
  GENE_LINEAR = 18 + 0.72 * qpcr_percent + stats::rnorm(44L, sd = 7),
  GENE_QUADRATIC = 22 + 0.10 * qpcr_percent + 0.0022 * qpcr_percent^2 +
    stats::rnorm(44L, sd = 8),
  GENE_LL4 = 30 + 180 / (1 + exp(-3.7 * (log(qpcr_percent) - log(118)))) +
    stats::rnorm(44L, sd = 5)
)
long_data = do.call(rbind, lapply(names(expression_by_gene), function(gene) {
  data.frame(
    sample_data,
    gene = gene,
    bulk_rna_percent = expression_by_gene[[gene]],
    stringsAsFactors = F
  )
}))

output_directory = dirname(output_file)
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = T, showWarnings = F)
}
utils::write.csv(long_data, output_file, row.names = F)
cat(sprintf("Wrote %d long-form sample-by-gene rows for %d genes to %s\n", NROW(long_data), length(expression_by_gene), output_file))
