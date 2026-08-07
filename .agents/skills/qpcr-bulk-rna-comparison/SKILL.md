---
name: qpcr-bulk-rna-comparison
description: Build, revise, or validate a publication-ready qPCR versus bulk RNA-seq expression concordance analysis in this repository. Use when sample-level measurements must be aggregated by gRNA and compared with linear, quadratic, and drc::LL.4 models, including ANOVA, AIC tables, and a raw-percent scatterplot.
---

# qPCR–bulk RNA-seq comparison

Use this workflow for raw normalized expression percentages where the control mean is 100%.

## Run the analysis

1. For batch DEG analysis, keep one long-form row per sample and gene with a sample identifier, gRNA identifier, normalized qPCR value, gene, and normalized bulk RNA-seq value.
2. Require exactly four samples per gRNA and 11 gRNA-level points for every gene unless the study design explicitly differs. Use `fit_qpcr_bulk_rna_deg_batch()` to enforce the grouping before fitting.
3. Fit exactly three model families with `fit_qpcr_bulk_rna_models()`: linear `lm`, quadratic `lm`, and `drc::drm(..., fct = drc::LL.4(fixed = c(NA, min(y), max(y), NA)))` with bulk RNA-seq asymptotes fixed to the fitted gRNA-mean range.
4. Preserve the original percentage scale. Do not add log scaling, Spearman correlation, or spline models unless a separate analysis explicitly requires them.
5. Write the AIC/AICc/BIC table for all three models. Write individual ANOVA tables for the linear and quadratic fits, plus the nested linear-versus-quadratic ANOVA. Do not describe this nested test as a comparison involving LL.4.
6. Write one multi-page PDF with one gene per page; each page must show its 11 gRNA-level means, 100% reference lines, and all three fitted curves.

## Use the command-line runner

For one pre-aggregated target gene, use `analysis/run_qpcr_bulk_rna_comparison.R`. For all DEGs, run from the repository root:

```sh
Rscript analysis/run_qpcr_bulk_rna_deg_batch.R \
  --input path/to/input_long.csv \
  --sample-id-col sample_id \
  --grna-col grna \
  --qpcr-col qpcr_percent \
  --gene-col gene \
  --bulk-rna-col bulk_rna_percent \
  --output-dir results/qpcr_bulk_rna_degs
```

Inspect `gene_grna_level_means.csv` before interpreting any model: it is the 11-point dataset actually fitted for every gene. Retain the gene-labelled AIC, ANOVA, coefficient, prediction, multi-page PDF, model-object, manifest, and session-information outputs.

## Validate

Run `make test` after changing shared model code. For a complete batch smoke test, run `make qpcr-bulk-deg-example`; it creates deterministic long-form data, validates 11 points per gene, and writes a three-page PDF. Treat LL.4 preference by AIC cautiously: asymptotes are fixed, but 11 points still offer limited support for the remaining nonlinear parameters.
