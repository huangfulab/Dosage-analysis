# qPCR--bulk RNA-seq DEG comparison

This workflow compares normalized qPCR expression with normalized bulk RNA-seq
expression for every selected differentially expressed gene (DEG). It is a
batch analysis: each gene receives its own 44-sample to 11-point aggregation,
three model fits, and one page in one multi-page PDF.

## Long-form input

Supply one CSV row per sample and gene. For 44 samples and (G) genes, the
input has (44 \times G) rows. The qPCR measurement is repeated for each gene
because it is the common x-axis variable.

| Column | Meaning |
| --- | --- |
| `sample_id` | Stable identifier; unique within each gene |
| `grna` | gRNA species used to aggregate four samples |
| `qpcr_percent` | qPCR expression normalized so controls average 100% |
| `gene` | Bulk RNA-seq gene identifier |
| `bulk_rna_percent` | That gene's expression normalized so controls average 100% |

The runner checks **each gene independently**. It requires exactly four
complete samples per gRNA and exactly eleven resulting gRNA-level points. A
gene that fails either condition stops the run rather than being silently
omitted. It also requires every selected gene to contain the same sample IDs,
gRNA assignments, and qPCR values.

Optionally provide `--deg-file` with a CSV gene list. By default, every gene
represented in the long-form input is analysed.

## Prespecified analysis per gene

For every gene, the runner averages qPCR and bulk RNA-seq percentages within
each gRNA and fits exactly:

1. Linear model.
2. Quadratic model.
3. Four-parameter log-logistic model: `drc::drm(y ~ x, fct=drc::LL.4())`.

All fitting and plotting use the original percentage scale. No log-scale
transformation, spline, or correlation statistic is used. The outputs include
each model's AIC (plus AICc/BIC), individual linear and quadratic ANOVA tables,
and the nested linear-versus-quadratic ANOVA. The nested test does not compare
LL.4 with either linear model. Replicate summaries additionally report within-
gRNA mean, SD, SE, minimum, and maximum for each assay; SD/SE are NA when fewer
than two observations are available. Scatterplot SE bars are descriptive and
do not claim inferential independence among gRNAs.

## Run

From the repository root:

```bash
make qpcr-bulk-deg-example
```

For study data:

```bash
Rscript analysis/run_qpcr_bulk_rna_deg_batch.R \
  --input data/deg_expression_long.csv \
  --sample-id-col sample_id \
  --grna-col grna \
  --qpcr-col qpcr_percent \
  --gene-col gene \
  --bulk-rna-col bulk_rna_percent \
  --deg-file data/differentially_expressed_genes.csv \
  --deg-gene-col gene \
  --output-dir results/qpcr_bulk_rna_degs
```

Omit `--deg-file` to analyse every gene in `--input`.

## Outputs

| Output | Contents |
| --- | --- |
| `gene_grna_level_means.csv` | The 11 fitted gRNA-level points plus replicate mean/SD/SE/min/max summaries for every gene |
| `gene_grna_replicate_summaries.csv` | Explicit gene-labelled replicate summary table |
| `gene_model_aic.csv` | Gene-labelled AIC, AICc, BIC, and rankings for all three models |
| `gene_model_coefficients.csv` | Gene-labelled coefficients and standard errors |
| `gene_model_predictions.csv` | Raw-scale fitted values for converged model curves; failed LL.4 rows remain diagnostic/table output rather than a fitted curve |
| `gene_model_diagnostics.csv` | One row per gene/gRNA/model with observed, fitted, residual; lm diagnostics are included where defined |
| `gene_leave_one_grna_out.csv` | Gene-labelled omission/refit AIC, rank, RMSE, convergence, and failure messages |
| `gene_model_summary.csv` | Per-gene model support, AIC winner/tie labels, RMSE, and explicit failures |
| `batch_concordance_difference_matrix.csv` | Exact gene-by-gRNA matrix of bulk RNA-seq minus qPCR percentage points |
| `batch_concordance_difference_long.csv` | Long-form version of the same descriptive disagreement values |
| `gene_*anova.csv` | Gene-labelled linear, quadratic, and nested ANOVA tables |
| `qpcr_bulk_rna_all_degs.pdf` | One page per gene: points, descriptive SE bars, reference lines, and curves only for converged models; LL.4 failures are shown in tables/diagnostics, not as fabricated curves |
| `qpcr_bulk_rna_all_degs_diagnostics.pdf` | Separate multi-page observed-vs-fitted and residual-vs-fitted diagnostic panels |
| `qpcr_bulk_rna_batch_concordance_heatmap.pdf` | Symmetric diverging heatmap of descriptive percentage-point disagreement |
| `qpcr_bulk_rna_model_overview.pdf` | Per-gene AIC weights/winners and RMSE; failures and ties are labelled |
| `batch_model_objects.rds` | Complete fitted objects for all genes, diagnostics, and sensitivity results |
| `run_manifest.tsv` | Input columns, gene/page count, package version, and Git state |

The synthetic example is a software fixture only and must not be interpreted as
study data. An AIC-preferred LL.4 curve is descriptive, not evidence of a
biological saturation mechanism or causal qPCR-to-RNA-seq relationship. The
concordance heatmap is specifically bulk RNA-seq minus qPCR in percentage
points: it describes disagreement and is not mechanistic truth. Model overview
weights and winners likewise summarize predictive support only.
Diagnostic values describe the fitted observations and are not a substitute for
independent validation. Leave-one-gRNA-out results probe sensitivity to the
specified gRNA points; they are not external validation, and failures are
reported rather than omitted.
