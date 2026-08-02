#!/usr/bin/env Rscript

# Generate a deterministic toy dataset used to verify the complete workflow.
# This is not study data and must not be used for scientific interpretation.

arguments <- commandArgs(trailingOnly = TRUE)
output_file <- if (length(arguments) >= 1L) {
  arguments[[1L]]
} else {
  file.path("data", "example_dosage.csv")
}

# Fix simulation noise so this 32-targeted/8-control regression fixture is
# reproducible; its dose range exercises control-mean normalization.
set.seed(20260721)

targeted_dose <- seq(0.18, 0.95, length.out = 32L)
control_dose <- c(0.96, 0.98, 0.99, 1.00, 1.00, 1.01, 1.02, 1.04)
dose <- c(targeted_dose, control_dose)
group <- c(rep("Targeted", length(targeted_dose)), rep("Controls", length(control_dose)))
guide <- c(
  sprintf("target-guide-%02d", seq_along(targeted_dose)),
  sprintf("control-guide-%02d", seq_along(control_dose))
)

# Simulate contrasting LL.4-like and quadratic endpoints to exercise candidate
# discrimination. Constants and noise are test settings, not biological estimates.
hill_mean <- 0.45 + 1.35 /
  (1 + exp(-5.5 * (log(dose) - log(0.58))))
quadratic_mean <- 0.65 + 0.15 * dose + 1.05 * dose^2

example_data <- data.frame(
  guide = guide,
  group = group,
  dose = dose,
  response_hill = hill_mean + stats::rnorm(length(dose), sd = 0.035),
  response_quadratic = quadratic_mean + stats::rnorm(length(dose), sd = 0.035),
  stringsAsFactors = FALSE
)

output_directory <- dirname(output_file)
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
}
utils::write.csv(example_data, output_file, row.names = FALSE, quote = TRUE)
cat(sprintf("Wrote %d rows to %s\n", NROW(example_data), output_file))
