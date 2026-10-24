# Neuroectoderm bulk RNA-seq: K2LLPS 3′-seq DESeq2 analysis

This study-specific analysis was added during manuscript review. It links the supplied K2LLPS gene-level expression matrix to the OCT4 guide design workbook and tests the pooled OCT4-targeting group against pooled nontargeting and safe-harbor controls.

## Location, provenance, and path conventions

This standalone workflow lives in `studies/neuroectoderm_bulk_rnaseq/` within **Dosage-analysis**, separate from the generic runners in the repository-root `analysis/` directory. Imported from [hyein-cbio/K2LLPS_DESeq2](https://github.com/hyein-cbio/K2LLPS_DESeq2), commit `404fea0c496b6eec66f7586046c4b9064f8d0560`, with input-path and documentation adaptations. The core model formulas are retained; validation and edge-case fixes made during relocation are listed below.

The study directory includes seven R scripts, this README, `.gitignore`, and `sample_metadata_draft.csv`. A synthetic integration test is supplied at repository-root `tests/test_neuroectoderm_workflow.R`. **Expression matrices, `260320_OCT4_qPCR.xlsx`, and generated outputs are not included.** The relocated scripts have been exercised with synthetic inputs, not the original study data. All study-specific numerical results described below refer to the original analysis, not a new study run in this repository.

All shell examples below assume this working directory, starting from the **Dosage-analysis repository root**:

```bash
cd studies/neuroectoderm_bulk_rnaseq
```

Default inputs are **placeholder locations** under this study directory; supply the real files before running:

```text
studies/neuroectoderm_bulk_rnaseq/
├── inputs/                              # not supplied; ignored by Git
│   ├── K2LLPS-expression-matrix.csv
│   ├── K2LLPS-expression-matrix.tsv
│   └── 260320_OCT4_qPCR.xlsx
├── sample_metadata_draft.csv             # supplied sample crosswalk
└── outputs/                             # generated; ignored by Git
    └── guide-specific/
```

Scripts locate default paths relative to their own file, not the caller's working directory. Explicit positional paths are interpreted relative to the caller's working directory unless absolute. Outputs therefore stay under this study's `outputs/` by default, without using the repository-root `data/`, `results/`, or `analysis/` directories. Most existing output names and rerun behaviour are retained; use a new output directory when preserving an earlier run. Scripts 00/01/05/06 can replace existing same-name files without prompting, and script 03 replaces the entire selected guide-specific directory after staging. Do not use these commands as backup operations. Pooled-DEG scatter outputs now have a `-pooled` suffix so they cannot overwrite the default union-DEG scatter outputs.

For external inputs, this alternative sequence uses `/path/to/...` as placeholders to replace (run from this study directory):

```bash
Rscript 00_verify_expression_files.R \
  /path/to/K2LLPS-expression-matrix.csv \
  /path/to/K2LLPS-expression-matrix.tsv outputs
Rscript 01_prepare_deseq2.R \
  /path/to/K2LLPS-expression-matrix.tsv sample_metadata_draft.csv outputs \
  /path/to/260320_OCT4_qPCR.xlsx
Rscript 02_verify_pooled_deseq2_outputs.R outputs /path/to/260320_OCT4_qPCR.xlsx
```

Use the same output directory for verification and downstream steps. For a custom output directory, pass it to every script; scripts 03/04 additionally accept a guide-specific output directory. Script 05 expects guide-specific results under `<output-directory>/guide-specific/` in its default union-DEG mode, so retain that layout when using it.

## Analysis design

The expression dataset contains 44 libraries: 11 guides measured in four shared experimental blocks (`rep_2` through `rep_5`). The workbook `260320_OCT4_qPCR.xlsx` maps each numeric HZ guide ID to its biological JY label and qPCR measurement.

The primary comparison is:

- **OCT4 targeting, 8 guides / 32 libraries:** `guide01`, `guide02`, `guide03`, `guide04`, `guide05`, `guide06`, `guide14`, `guide15`
- **Pooled controls, 3 guides / 12 libraries:** nontargeting `guide18` and `guide19`, plus safe-harbor `guide20`

The DESeq2 model is:

```r
~ replicate + comparison_group
```

`replicate` is categorical and represents the confirmed shared experimental block. `pooled_control` is the reference level, so a **positive log2 fold change means higher expression in the OCT4-targeting libraries**.

The workbook qPCR values are linked to every RNA-seq library for annotation and QC. They are post-treatment response measurements and are not included as model covariates. The expected single missing qPCR value is `guide04_rep_5`.

## Input handling

- Expression inputs (not included): `inputs/K2LLPS-expression-matrix.csv` and `inputs/K2LLPS-expression-matrix.tsv` for file verification; the TSV is used for DESeq2.
- Stable count-column/sample crosswalk: `sample_metadata_draft.csv`
- Guide design and qPCR workbook (not included): `inputs/260320_OCT4_qPCR.xlsx`, `Sheet1`.
- Only `*_count` columns are analyzed. The TSV must also contain matching `*_cpm` sample columns for schema validation; their values are not used in DESeq2.
- This is 3′-seq, so no gene-length normalization, TPM/FPKM conversion, or transcript-length offset is applied.
- Fractional counts are expected from featureCounts multi-mapper allocation.
- Workbook rows are never joined by position. HZ IDs are normalized to `guide%02d`, then linked by the explicit `guide + replicate` key.

## Dependencies

Beyond Bioconductor **DESeq2**, the pooled pipeline also needs:

- **openxlsx**, **devtools** (`devtools` also records session information for the standalone verifiers and scatter runner)
- **ggplot2**, **ggrepel** (PCA eigenvector arrow labels)
- **pheatmap** (DEG heatmap)
- **matrixStats** (via DESeq2 stack; used for row variances / SDs)
- **GSVA** (neuroectoderm ssGSEA / GSVA scores on the existing sample PCA)
- **drc** (fixed-asymptote LL.4 fitting in script 05)

Install required packages explicitly in your chosen environment; these scripts do not install them. Script 06 requires the modern GSVA parameter-object API (`ssgseaParam`, `gsvaParam`), not the older `gsva(expression, geneSets)` API.

The synthetic integration test was run with R 4.6.1, DESeq2 1.52.0, openxlsx 4.2.9, devtools 2.5.2, ggplot2 4.0.3, ggrepel 0.9.8, pheatmap 1.0.13, matrixStats 1.5.0, GSVA 2.6.6, and drc 3.0.1. This is an observed test environment, not a lockfile or a claim that other versions reproduce the study results. The optional `glmGamPoi` strategy was not runtime-tested because that package was unavailable.

## Verify the two expression files

```bash
Rscript 00_verify_expression_files.R
```

The verifier permits CSV-versus-TSV delimiters and renamed sample headers, but requires identical dimensions and annotations and compares every expression value. It writes:

- `outputs/expression-file-equivalence.txt`
- `outputs/expression-header-map.csv`

It exits with a failure status if the files are not practically equivalent: annotations must match exactly and every finite expression value must differ by no more than the absolute numeric tolerance (default `1e-14`). Optional positional arguments are `[csv] [tsv] [output-directory] [numeric-tolerance]`.

## Run differential expression

```bash
Rscript 01_prepare_deseq2.R
```

Optional positional arguments are:

```text
Rscript 01_prepare_deseq2.R [expression.tsv] [metadata.csv] [output-directory] [design-workbook.xlsx]
```

### Count strategy and thresholds

The default `round` strategy rounds fractional featureCounts weights to the nearest integer and uses DESeq2's standard parametric dispersion fit. It does not reconstruct counts from CPM.

To retain fractional values, install `glmGamPoi` separately and run:

```bash
K2LLPS_COUNT_STRATEGY=glmGamPoi Rscript 01_prepare_deseq2.R
```

Environment knobs:

| Variable | Default | Role |
|---|---:|---|
| `K2LLPS_MIN_TOTAL_COUNT` | `10` | Prefilter: total analysis count across samples |
| `K2LLPS_ALPHA` | `0.05` | Wald FDR threshold (`padj`) for significance tables and strict DEGs |
| `K2LLPS_DEG_LFC` | `1` | Absolute **raw** log2FC cutoff for the strict DEG heatmap/violin set |

```bash
K2LLPS_ALPHA=0.01 K2LLPS_DEG_LFC=1 Rscript 01_prepare_deseq2.R
```

The script estimates size factors and dispersions, runs a Wald test for `OCT4_targeting` minus `pooled_control`, and applies DESeq2 `type="normal"` LFC shrinkage.

- **Primary FDR table** (`de-results-...-significant.csv`): nonmissing `padj < K2LLPS_ALPHA`. Shrunken LFC is used for ranking and direction labels; p-values remain from the Wald test.
- **Strict visualization DEGs** (heatmap / DEG violins): nonmissing `padj <= K2LLPS_ALPHA` **and** `|raw log2FoldChange| >= K2LLPS_DEG_LFC`. Raw (unshrunken) LFC is used so the hard effect-size cutoff is not deflated by the normal prior. The original study run with defaults reported **2681 genes (2073 up, 608 down)**; this is not an expected result for arbitrary or synthetic inputs.

## PCA, loadings, and gRNA species order

Sample PCA and all loading figures use the same representation:

1. **VST** on the fitted DESeq2 object (`varianceStabilizingTransformation(..., blind=FALSE)`).
2. Top **500** genes by VST row variance (DESeq2 `plotPCA` default `ntop=500`).
3. `prcomp` on samples × those genes (center = TRUE, scale. = FALSE).

Gene loadings are the `rotation` columns (PC1/PC2). Positive PC1 loading pushes samples rightward on the sample PCA.

### PC1 / PC2 loading displays

For each of PC1 and PC2:

- Bar plot of the **top 30 and bottom 30** loadings on that axis.
- Native PC1–PC2 **arrow plot** for the **same gene set** as that bar plot (no sample points; `ggrepel` labels).

### Canonical gRNA species order (PC1 median)

Nine species are defined as:

- `pooled_control` = all libraries from `guide18`, `guide19`, and `guide20`
- eight individual OCT4 guides: `guide01`, `guide02`, `guide03`, `guide04`, `guide05`, `guide06`, `guide14`, `guide15`

**Final order source:** sample PC1 scores aggregated by species, sorted by **median PC1 ascending**.
That order is written to `outputs/pca-pc1-score-by-grna-species-medians.csv` and reused for:

- PC1 score violin (`pca-pc1-score-violin-by-grna-species.png`)
- DEG violin x-axis (both up and down panels)
- DEG heatmap **column** order (primary = this species order; secondary = replicate 2→5)

In the original study run the median order was:

`pooled_control → guide05 → guide03 → guide06 → guide15 → guide04 → guide02 → guide14 → guide01`

## Strict DEG heatmap and violins

- Expression matrix: VST, restricted to strict DEGs, **row Z-scored** across the 44 samples.
- **Rows:** `hclust(method = "ward.D2")` on distance `1 - correlation`.
- **Columns:** PC1-median gRNA species order, then replicate (full-rank blocks of four libraries per targeting guide; twelve pooled-control libraries as one species block).
- Heatmap column annotations: `replicate`, `gRNA_species`, `comparison_group` only. There is **no** separate per-library `guide` track, so controls stay a single `pooled_control` species block (guide18/19/20 are not split in the color bar).
- DEG violins: **18** violins = up/down × (8 OCT4 guides + `pooled_control`). Y-axis is the same row Z-score values; x-axis uses the PC1-median species order.

## Verify the generated analysis

```bash
Rscript 02_verify_pooled_deseq2_outputs.R
```

Optional positional arguments are `[output-directory] [design-workbook.xlsx]`.

The verifier independently checks workbook/sample linkage, the 8-versus-3 guide assignment, four shared blocks, model rank and reference level, coefficient direction, retention of all 44 libraries, result-table consistency, significance flags, normalized-count means, PCA loadings/arrows, strict DEG membership, DEG violin group structure, and required artifacts. It also confirms that `guide15_rep_4` remains included despite its lower library depth.

As a replicate-sensitivity sanity check, the verifier also fits `~ replicate + comparison_group + replicate:comparison_group` and performs an omnibus likelihood-ratio test against the additive model `~ replicate + comparison_group`. This jointly tests whether any of the three replicate-by-group interaction coefficients is nonzero. At its fixed FDR threshold of 0.05, independently of `K2LLPS_ALPHA`, the expected sole significant interaction gene is GATAD2B (`ENSG00000143614`). This interaction analysis is a verification check and does not replace the primary additive Wald-test results. **Script 02 is therefore an original-study regression verifier, not a general-purpose check for any valid 44-sample dataset.** Supplying schema-valid files alone does not guarantee that this biological expectation will hold. The synthetic integration test explicitly simulates this one interaction rather than bypassing the check.

The strict DEG LFC cutoff is read from the generating run's `run-summary.txt`, not hard-coded to 1 or taken from the verifier's environment. Old outputs without the explicit cutoff record must be regenerated into a new directory. The verifier also preserves the recorded count strategy when reconstructing its interaction model.

## Outputs

### Metadata and QC

- `sample_metadata_draft.csv`: stable source crosswalk from each TSV count column to its biological sample ID; never overwritten when it already exists.
- `outputs/analysis-sample-metadata.csv`: finalized 44-row sample table linked to workbook labels, guide class, comparison group, block, and qPCR values.
- `outputs/sample_qc.csv`: library sums, rounding change, depth ratio, size factor, Cook's-distance summaries, workbook annotations, and inclusion status.

### DESeq2 objects and results

- `outputs/DESeq2Obj-after-wald-test.rds`: fitted DESeqDataSet after the Wald test.
- `outputs/de-results-oct4-targeting-vs-pooled-controls-complete.csv`: every retained gene, including genes with missing p-values or adjusted p-values.
- `outputs/de-results-oct4-targeting-vs-pooled-controls-significant.csv`: FDR-significant subset (`padj < K2LLPS_ALPHA`), ranked by adjusted p-value and absolute shrunken effect.
- `outputs/degs-oct4-targeting-vs-pooled-controls-lfc1-padj0.05.csv`: strict DEGs for visualization (`|raw LFC| >= K2LLPS_DEG_LFC` and `padj <= K2LLPS_ALPHA`). This historical filename is fixed even when thresholds change; use `run-summary.txt` for the actual cutoffs.

Both primary result tables contain gene annotations, normalized-count means by group, raw and shrunken effects, Wald inference, filtering/test status, FDR flags, and direction. In the FDR table, `higher_in_OCT4_targeting` / `lower_in_OCT4_targeting` follow the **shrunken** effect sign.

### PCA and loadings

- `outputs/pca-coordinates.csv` / `pca-oct4-targeting-vs-pooled-controls.png`: sample PCA (VST, top 500 variable genes).
- `outputs/pca-loadings-oct4-targeting-vs-pooled-controls.csv`: PC1/PC2 loadings for those 500 genes.
- `outputs/pca-pc1-loadings-top-bottom-30-oct4-targeting-vs-pooled-controls.{csv,png}`: PC1 top/bottom 30 bar plot + table.
- `outputs/pca-pc2-loadings-top-bottom-30-oct4-targeting-vs-pooled-controls.{csv,png}`: PC2 top/bottom 30 bar plot + table.
- `outputs/pca-pc1-eigenvector-arrows-oct4-targeting-vs-pooled-controls.{csv,png}`: PC1 extreme genes as native loading arrows (ggrepel).
- `outputs/pca-pc2-eigenvector-arrows-oct4-targeting-vs-pooled-controls.{csv,png}`: PC2 extreme genes as native loading arrows (ggrepel).
- `outputs/pca-pc1-score-by-grna-species.csv`: per-sample PC1 with 9-species labels.
- `outputs/pca-pc1-score-by-grna-species-medians.csv`: **canonical species order** (sorted by median PC1).
- `outputs/pca-pc1-score-violin-by-grna-species.png`: PC1 score violins in that order.

### Strict DEG visualizations

- `outputs/deg-heatmap-vst-rowz-oct4-targeting-vs-pooled-controls.png`: VST row-Z heatmap; row cluster ward.D2 on `1 - cor`; columns by PC1-median species then replicate.
- `outputs/deg-violin-oct4-targeting-vs-pooled-controls.png`: 18 violins (up/down × 9 species), same species order.
- `outputs/deg-violin-data-oct4-targeting-vs-pooled-controls.csv`: long-form row-Z values for the violins.

### Other diagnostics and provenance

- `outputs/dispersion-estimates.png`
- `outputs/ma-oct4-targeting-vs-pooled-controls.png`
- `outputs/run-summary.txt`: inputs/checksums, design, coefficient, filtering, DE totals, and strict DEG counts.
- `outputs/sessionInfo.txt`: R and package versions.

All 44 samples are retained by default. QC diagnostics may motivate a separately documented sensitivity analysis, but the pipeline does not silently remove samples.

## Guide-specific contrasts

The pooled analysis above remains the primary aggregate comparison and is not modified by the guide-specific extension. The extension reads the verified pooled fitted object as a source of counts, annotations, metadata, sample order, and the retained prefilter universe (33,480 genes in the original study run), then reconstructs and refits a separate DESeq2 object with:

```r
~ replicate + contrast_group
```

`contrast_group` has `pooled_control` as the reference, followed by `guide01`, `guide02`, `guide03`, `guide04`, `guide05`, `guide06`, `guide14`, and `guide15`. Guides `guide18`, `guide19`, and `guide20` form the pooled-control level. The shared 44-sample model therefore estimates eight named guide-minus-control coefficients while adjusting for the categorical replicate blocks 2–5. Positive log2 fold change means higher expression in the named target guide.

Run the extension after the pooled verifier passes:

```bash
Rscript 02_verify_pooled_deseq2_outputs.R
Rscript 03_run_guide_specific_deseq2.R
Rscript 04_verify_guide_specific_outputs.R
```

Optional positional arguments are:

```text
Rscript 03_run_guide_specific_deseq2.R [pooled-output-directory] [guide-specific-output-directory]
Rscript 04_verify_guide_specific_outputs.R [pooled-output-directory] [guide-specific-output-directory]
```

`K2LLPS_ALPHA` controls the per-guide FDR threshold and defaults to 0.05. Each guide uses its own DESeq2 Benjamini-Hochberg adjusted p-values; shrunken effects determine direction and ranking, not significance.

### Guide-specific outputs

All new files are transactionally published under `outputs/guide-specific/`; the existing pooled filenames are never overwritten.

- `de-results-<guide>-vs-pooled-controls-complete.csv`: every retained gene for one guide-minus-pooled-control coefficient.
- `de-results-<guide>-vs-pooled-controls-significant.csv`: the exact ordered subset with nonmissing `padj < K2LLPS_ALPHA`.
- `concordance-across-guides.csv`: one row per retained gene with tested/significant guide counts, fixed-order guide lists, shrunken-effect summaries, mixed-direction flags, and recurrence status.
- `DESeq2Obj-guide-specific-after-wald-test.rds`: fitted DESeqDataSet for the separate guide-specific model.
- `pooled-source-md5.csv`: the pooled artifacts and hashes used by the run; the analysis and verifier fail if those files change.
- `run-summary.txt` and `sessionInfo.txt`: model provenance, per-guide totals, concordance totals, and package versions.

A gene is flagged `consistent_across_multiple_guides` only when it is FDR-significant in at least two target guides in the same direction and no target guide is FDR-significant in the opposite direction. The table retains all guide counts and lists so alternative recurrence thresholds can be applied without rerunning DESeq2. `mixed_significant_directions` identifies genes with significant positive and negative effects across different guides.

The guide-specific verifier independently reconstructs all eight raw and shrunken result sets, normalized-count means, flags, table ordering, and concordance fields. It also confirms 44/44 sample retention, inclusion of `guide15_rep_4`, the 44-by-12 full-rank design with 32 residual degrees of freedom, and byte-level agreement with the recorded pooled source outputs.

## Gene expression vs OCT4 qPCR dosage modeling

`05_gene_vs_qpcr_scatter.R` fits linear, quadratic, and fixed-asymptote LL.4 curves to 11 guide-mean points (8 OCT4 gRNAs + 3 controls). The default DEG batch is the **union of pairwise strict DEGs**: a gene is included if it is significant in at least one OCT4 gRNA versus the pooled-control average at `padj <= K2LLPS_ALPHA` and `|raw log2FoldChange| >= K2LLPS_DEG_LFC`. Direction is `up`/`down` when all contributing pairwise contrasts agree; mixed-direction genes are flagged `direction_status=mixed` and assigned a primary direction from the median shrunken LFC.

Unlike the generic raw-percentage qPCR workflow in the repository-root `analysis/`, this study script transforms bulk expression **before** aggregation: `n = DESeq2 normalized counts`, then `mean(log2(n + 1))` within each guide. Its default y-axis is `100 * guide_mean_log_expression / mean(control_guide_mean_log_expression)`; x is qPCR scaled to control-guide mean 100. These y-values are **relative log-expression, not percent RNA abundance**. If all control-guide counts for a gene are zero, y remains absolute mean log-expression and is labelled accordingly. qPCR averages omit the documented missing `guide04_rep_5` value (three qPCR measurements for that guide, four RNA-seq libraries); this does not drop the RNA library.

The union is reconstructed from the eight **complete** guide result tables, so a downstream `K2LLPS_ALPHA` change and the inclusive `padj <= alpha` boundary are not truncated by an earlier significant-only export. Pooled mode uses the already generated strict DEG list and its recorded thresholds instead.

Shared gene symbols that map to multiple Ensembl IDs are collapsed to a single row. The retained ID prefers `protein_coding`, then higher `baseMean`, then smaller minimum pairwise `padj`. The discarded IDs are recorded in `scatter-deg-duplicate-resolution.csv`.

The AIC tables do **not** auto-assign Figure 4A Linear/NA labels. They add the linear-model slope p-value (`pvalue_linear`) and the quadratic coefficients of `y = a x^2 + b x + c` so models can be called manually. `quadratic_vertex_x` and `quadratic_vertex_in_range` indicate whether the parabola turns inside the observed qPCR range.

```bash
Rscript 05_gene_vs_qpcr_scatter.R
Rscript 05_gene_vs_qpcr_scatter.R outputs --csv-only
Rscript 05_gene_vs_qpcr_scatter.R outputs --batch-deg-pooled
Rscript 05_gene_vs_qpcr_scatter.R outputs --batch-all
Rscript 05_gene_vs_qpcr_scatter.R outputs SOX1
```

### Scatter / modeling outputs

- `outputs/degs-union-oct4-guides-vs-pooled-controls-lfc1-padj0.05.csv`: union DEG membership, pairwise guide lists, and collapse metadata.
- `outputs/scatter-deg-duplicate-resolution.csv`: symbols that mapped to more than one Ensembl ID and the retained ID.
- `outputs/scatter-deg-up-vs-qpcr-model-aic.csv` / `scatter-deg-down-vs-qpcr-model-aic.csv`: per-direction fit tables.
- `outputs/scatter-deg-vs-qpcr-model-aic.csv`: combined up+down table.
- `outputs/scatter-deg-up-vs-qpcr-model-fits.pdf` / `scatter-deg-down-vs-qpcr-model-fits.pdf`: one page per gene.
- `outputs/scatter-all-genes-vs-qpcr-model-aic.csv` / `scatter-all-genes-vs-qpcr-model-fits.pdf`: optional all-gene batch.

The `--batch-deg-pooled` mode appends `-pooled` before `.csv`/`.pdf` on the three DEG AIC tables and two model-fit PDFs. A finite preview cap adds `-preview<N>` after that scope suffix. One-sided DEG lists produce a header-only AIC table for the empty direction, with the same schema as the non-empty direction.

## Neuroectoderm GSVA / ssGSEA on the sample PCA

`06_gsva_neuroectoderm_pca.R` scores the 44-sample variance-stabilized matrix with one compact gene set, then colors the **existing** sample PCA (`outputs/pca-coordinates.csv`) by those scores.

Gene set: Enrichr **CellMarker_2024 / Neural Stem Cell Brain Human** (24 symbols, including SOX1, SOX2, PAX6, NES, HES5). Downloaded 2026-09-08 from `https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=CellMarker_2024`. Named “neuroectoderm” alternatives were either 4 genes (too small for ssGSEA) or GO:0060563 (35 genes, missing SOX1/PAX6/NES). Primary database: CellMarker (Zhang et al. 2019 *NAR*; Hu et al. 2023 *NAR*). The script stores this as `gene_set_provenance` and writes `outputs/gsva-neuroectoderm-gene-set-provenance.txt`.

Two algorithms from Bioconductor **GSVA**:

- `ssgseaParam(..., normalize=TRUE, alpha=0.25)` — Barbie et al. 2009 ssGSEA
- `gsvaParam(..., kcdf="Gaussian")` — Hänzelmann et al. 2013 GSVA, Gaussian kernel for VST values

The expression input is the **full** fitted VST matrix (all retained genes), not the gene-set rows alone. PCA x/y are the stored PC1/PC2 coordinates; they are not recomputed for plotting. GSVA fill is the signed score, blue–white–red, midpoint 0, symmetric limits at `max(|score|)`. ssGSEA Barbie normalization does **not** guarantee nonnegative scores; that panel uses grey50–yellow–red over the observed `min(score)` to `max(score)`, retaining any negative values. Constant-score plotting limits are padded without changing the stored scores.

```bash
Rscript 06_gsva_neuroectoderm_pca.R
Rscript 06_gsva_neuroectoderm_pca.R outputs
```

### GSVA outputs

- `outputs/gsva-neuroectoderm-gene-set.csv`: 24 symbols, mapped Ensembl IDs, and download source columns.
- `outputs/gsva-neuroectoderm-gene-set-provenance.txt`: Enrichr URL, retrieval date, CellMarker/Enrichr citations.
- `outputs/gsva-neuroectoderm-sample-scores.csv`: 44 rows with PC1, PC2, raw `ssgsea`, and signed `gsva`.
- `outputs/gsva-neuroectoderm-ssgsea-pca.png` / `gsva-neuroectoderm-gsva-pca.png`: one PCA per algorithm.
- `outputs/gsva-neuroectoderm-ssgsea-gsva-pca.png`: the two panels side by side.
- `outputs/gsva-neuroectoderm-run-summary.txt` / `gsva-neuroectoderm-sessionInfo.txt`

## Relocation validation and limits

From the **repository root**, run:

```bash
Rscript tests/test_neuroectoderm_workflow.R
```

The test copies the scripts into a temporary directory, generates 1,201 synthetic genes across the supplied 44-sample crosswalk plus a synthetic workbook, and executes all seven scripts from a different working directory. It checks the pooled and guide-specific verifiers, nondefault LFC cutoff 2, capped union/pooled/all-gene scatter modes, single-gene figures, and GSVA/ssGSEA outputs. Regressions cover expression tolerance acceptance/rejection, padded workbook labels, a run with no significant guides, an empty DEG direction, inclusive downstream FDR boundaries, negative ssGSEA scores, and rejection of pooled-output ancestor paths. Existing repository `data/`, `results/`, and the supplied crosswalk are fingerprinted to ensure they remain unchanged.

Relocation fixes include the tolerance verdict, recorded-threshold verification, preserving fractional-count mode in the verifier, complete-table union selection, matching empty-direction schemas, separate pooled-scatter filenames, deterministic PCA jitter, acceptance of signed ssGSEA scores, a guard against publishing guide-specific results over the pooled directory or its ancestors, consistent workbook-label trimming, and explicit text-column types for empty verification outputs.

This demonstrates conditional software execution with synthetic fixtures, **not reproduction of the original study counts, biological findings, or figures**. Real inputs must satisfy the crosswalk, workbook schema, exact missing-qPCR expectation, and the original-study interaction regression. The pooled plot/verification path also expects at least 500 retained variable genes, sufficient signed PCA loading contributors, and strict DEGs in both directions. Review experimental-unit independence, pooled-control assumptions, batch structure, and multiple testing separately before interpreting the scientific results. No raw-data accession or fixed environment is provided here yet.
