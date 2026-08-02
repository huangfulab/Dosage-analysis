#!/usr/bin/env Rscript

# Generate a deterministic synthetic fixture, not study data or biological evidence.
arguments = commandArgs(trailingOnly = T)
output_file = if (length(arguments) == 0L) {
  file.path("data", "example_qpcr_bulk_rna.csv")
} else if (length(arguments) == 1L) {
  arguments[[1L]]
} else {
  stop("Usage: Rscript scripts/generate_qpcr_bulk_rna_example.R [output.csv]", call. = F)
}

# Use the planned 11 gRNAs × 4 samples and a noisy LL.4-like relationship to
# exercise aggregation and all three candidate models under their defaults.
set.seed(20260727)
grna = sprintf("gRNA_%02d", seq_len(11L))
qpcr_means = c(28, 38, 52, 66, 82, 100, 125, 158, 200, 252, 318)
sample_data = do.call(rbind, lapply(seq_along(grna), function(index) {
  qpcr_values = qpcr_means[index] + stats::rnorm(4L, mean = 0, sd = 3)
  bulk_truth = 32 + 176 / (1 + exp(-3.7 * (log(qpcr_values) - log(118))))
  data.frame(
    sample_id = sprintf("%s_rep%d", grna[index], seq_len(4L)),
    grna = grna[index],
    qpcr_percent = qpcr_values,
    bulk_rna_percent = bulk_truth + stats::rnorm(4L, mean = 0, sd = 5),
    stringsAsFactors = F
  )
}))

output_directory = dirname(output_file)
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = T, showWarnings = F)
}
utils::write.csv(sample_data, output_file, row.names = F)
cat(sprintf("Wrote %d sample-level rows to %s\n", NROW(sample_data), output_file))
