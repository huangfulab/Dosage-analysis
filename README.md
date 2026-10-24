# OCT4 and NANOG dosage analysis

Reproducible R workflows for OCT4/NANOG dosage-response modelling and comparison
of qPCR with bulk RNA-seq expression. The repository contains standalone model
functions, command-line runners, synthetic data generators, and automated tests;
it is not an installable R package.

**Study data are not included.** Generated examples are synthetic software
fixtures, not experimental results or evidence for a biological mechanism.
The repository is not yet a complete manuscript reproduction archive: study-data
accessions, preprocessing, figure-to-command mappings, a locked environment,
and a versioned archive DOI remain to be finalized.

## Workflows

| Workflow | Runner | Analysis |
| --- | --- | --- |
| Dosage-response | `analysis/run_dosage_analysis.R` | One or more response columns; linear, quadratic, fixed-asymptote LL.4, and smoothing-spline candidates |
| Single-gene qPCR vs bulk RNA-seq | `analysis/run_qpcr_bulk_rna_comparison.R` | Sample-level input aggregated to gRNA means; linear, quadratic, and fixed-asymptote LL.4 candidates |
| Multi-gene qPCR vs bulk RNA-seq | `analysis/run_qpcr_bulk_rna_deg_batch.R` | Long-form sample-by-gene input; the same three models for each selected gene, with combined tables and multi-page PDFs |

The qPCR workflows use normalized expression percentages with a control mean
of 100%, without log transformation or additional control normalization. By
default, four samples per gRNA are averaged into eleven fitted gRNA-level
points for each gene. The batch workflow compares an externally selected gene
set; it does **not** perform differential-expression testing.

## Requirements and setup

- R with `Rscript` available on `PATH`.
- The CRAN package **`drc`** and its dependencies for LL.4 fitting. The remaining
  workflow code uses packages shipped with R, including `stats`, `graphics`,
  `grDevices`, `utils`, and `tools`.
- Git to clone the repository and capture commit/worktree provenance. Runners
  also work from extracted source archives, but Git metadata may be unavailable.
- `make` for the convenience targets; direct `Rscript` commands work without it.
- Optional: `pdftotext` (Poppler) enables an additional PDF-content assertion in
  the batch tests; that assertion is skipped when the executable is absent.

Clone the repository and run subsequent commands from its root:

```bash
git clone https://github.com/huangfulab/Dosage-analysis.git
cd Dosage-analysis
```

Install `drc` explicitly if needed; the runners do not install packages:

```bash
Rscript -e 'install.packages("drc", repos="https://cloud.r-project.org")'
```

Check your environment before running analyses:

```bash
Rscript --version
Rscript -e 'stopifnot(requireNamespace("drc", quietly=T)); print(utils::packageVersion("drc")); sessionInfo()'
```

The current development environment was inspected with R 4.6.1 and `drc` 3.0.1.
These are an environment reference, **not** a tested minimum-version range or a
locked dependency specification. No `renv.lock` or independently verified
clean-environment setup is supplied yet. Installing the latest CRAN packages
alone does not guarantee reproduction of an earlier run. The compatibility CI
uses explicit assumed minimum direct-package versions; see
[`ci/README.md`](ci/README.md) for the R 4.0–4.6 matrix, dependency floors,
validation status, and the distinction from a complete environment lockfile.

## Quick start without replacing existing data or results

The following commands use a POSIX shell (for example, Bash or Zsh). Run the
blocks in the **same shell**, from the repository root. First create a unique
workspace under `results/`; keep the printed path if you want to inspect it later:

```bash
mkdir -p results
demo_dir="$(mktemp -d results/readme-demo.XXXXXX)"
printf 'Demo workspace: %s\n' "$demo_dir"
```

All generated inputs and outputs below stay in that new workspace, leaving
`data/example_*.csv` and existing result directories untouched. If returning in
a new shell, create a new workspace first rather than reusing an unset variable.

### 1. Dosage-response example

```bash
Rscript scripts/generate_example_data.R "$demo_dir/example_dosage.csv"
Rscript analysis/run_dosage_analysis.R \
  --input "$demo_dir/example_dosage.csv" \
  --dose-col dose \
  --response-cols response_hill,response_quadratic \
  --group-col group \
  --control-label Controls \
  --normalize-dose true \
  --normalize-response true \
  --output-dir "$demo_dir/dosage"
```

This example explicitly divides dose and response by their respective control
means. The resulting control mean is **1**, not 100. Normalization is off by
default and is separate from the qPCR percentage-scale workflow.

### 2. Single-gene qPCR vs bulk RNA-seq example

```bash
Rscript scripts/generate_qpcr_bulk_rna_example.R "$demo_dir/example_qpcr_bulk_rna.csv"
Rscript analysis/run_qpcr_bulk_rna_comparison.R \
  --input "$demo_dir/example_qpcr_bulk_rna.csv" \
  --grna-col grna \
  --qpcr-col qpcr_percent \
  --bulk-rna-col bulk_rna_percent \
  --output-dir "$demo_dir/qpcr_bulk"
```

The fixture contains 44 sample-level rows. Inspect `grna_level_means.csv` in
the output directory to see the eleven points actually fitted.

### 3. Multi-gene qPCR vs bulk RNA-seq example

```bash
Rscript scripts/generate_qpcr_bulk_rna_deg_example.R "$demo_dir/example_qpcr_bulk_rna_degs.csv"
Rscript analysis/run_qpcr_bulk_rna_deg_batch.R \
  --input "$demo_dir/example_qpcr_bulk_rna_degs.csv" \
  --sample-id-col sample_id \
  --grna-col grna \
  --qpcr-col qpcr_percent \
  --gene-col gene \
  --bulk-rna-col bulk_rna_percent \
  --output-dir "$demo_dir/qpcr_bulk_degs"
```

The fixture contains three genes and 44 samples per gene (132 rows). Inspect
`gene_grna_level_means.csv` before interpreting the per-gene fits. The main PDF,
`qpcr_bulk_rna_all_degs.pdf`, contains one page per gene.

### Existing Make targets and overwrite safety

On a fresh checkout, `make example`, `make qpcr-bulk-example`, and
`make qpcr-bulk-deg-example` generate fixtures under `data/` and write to
`results/example`, `results/qpcr_bulk_rna_example`, and
`results/qpcr_bulk_rna_degs_example`, respectively.

**These targets are not safe rerun/backup commands:** the generators replace
their target input CSVs without prompting, even if the analysis then refuses
an existing non-empty output directory. Prefer the unique-workspace commands
above when preserving local inputs or previous runs.

All three runners default to `--overwrite false` and refuse non-empty output
directories. **`--overwrite true` recursively deletes the contents of the
selected output directory before analysis.** Use a new output path instead;
never point this option at a directory containing inputs or anything to retain.

## Using study data

Use CSV files with a header row and supply the actual column names through
CLI options. Full schemas and data-handling guidance are in
[`data/README.md`](data/README.md).

- **Dosage-response:** one row per declared analysis unit, one numeric dose
  column, and one or more numeric response columns. A group column and control
  label are needed when requesting control normalization. Missing/non-finite
  model values and constant responses are rejected by default; LL.4 requires
  strictly positive doses.
- **Single-gene qPCR:** sample-level `grna`, `qpcr_percent`, and
  `bulk_rna_percent` columns (or their explicitly named equivalents).
- **Multi-gene qPCR:** one row per sample and gene, with `sample_id`, `grna`,
  `qpcr_percent`, `gene`, and `bulk_rna_percent`. Every selected gene must have
  the same sample IDs, gRNA assignments, and qPCR values. qPCR values must be
  positive for LL.4. Incomplete or incorrectly grouped input stops the run.

The qPCR defaults require 44 complete samples per gene: eleven gRNAs with four
samples each. Override `--expected-replicates` and `--expected-points` only when
the experimental design differs, and document the change. The batch runner
accepts `--deg-file path/to/genes.csv` and `--deg-gene-col gene` to select genes;
without that file, it analyses all genes present in the input.

Document raw-data-to-CSV processing, units, biological/technical replicates,
normalization, exclusions, missingness, and batch structure before interpreting
study results. Do not commit confidential or controlled-access inputs.

See [`QPCR_BULK_RNA_WORKFLOW.md`](QPCR_BULK_RNA_WORKFLOW.md) for study-data
commands, model details, and the complete batch output inventory. Available
options can also be inspected without running an analysis:

```bash
Rscript analysis/run_dosage_analysis.R --help
Rscript analysis/run_qpcr_bulk_rna_comparison.R --help
Rscript analysis/run_qpcr_bulk_rna_deg_batch.R --help
```

## Outputs and provenance

| Output group | Contents |
| --- | --- |
| Fitted input | Dosage: `prepared_data.csv` and `normalization_factors.csv`; qPCR: sample input and gRNA-level means/replicate summaries |
| Model tables | AIC/AICc/BIC, coefficients, predictions, diagnostics, and linear-versus-quadratic nested tests; names vary by workflow |
| Sensitivity | qPCR leave-one-gRNA-out refits, convergence status, and failure messages |
| Figures | Model curves and diagnostic PDFs; batch qPCR additionally exports descriptive disagreement heatmaps and model-support overviews |
| Fitted objects | `model_objects.rds` or batch `batch_model_objects.rds` |
| Run record | `run_manifest.tsv`, `session_info.txt`, and `analysis_provenance.json` |

The run records capture settings, input/source MD5 fingerprints, Git commit
and dirty-worktree status, and session/package information. The JSON inventory
records sizes and MD5s of previously generated output files; it does not hash
itself. MD5 is used for integrity/provenance, not as a security guarantee.

A run record describes the code that **actually ran**, not the current checkout.
Older local results do not become current when code is updated. Before reusing
a result, compare its source fingerprints, settings, and environment with the
intended version; preserve earlier runs and regenerate into a separate folder.
PDF/RDS byte identity is distinct from numerical reproducibility, which needs
explicit reference values and tolerances.

`data/*` (except `data/README.md`) and `results/` are ignored by Git, including
the quick-start workspace. **Pushing code does not back up these files or put
them in a GitHub source archive.** Public data and reference results need an
explicit preservation plan, such as an appropriate dataset accession or an
approved archive deposit. Provenance files can include absolute local paths;
review them for private information before sharing.

## Validation and statistical interpretation

Run the seven model/runtime-state/CLI test suites from the repository root:

```bash
make test
```

The tests cover the dosage, single-gene qPCR, and batch qPCR workflows.
The GitHub Actions compatibility workflow additionally runs GSVA compatibility
and the synthetic neuroectoderm study integration test on every declared R
series, using its assumed minimum direct dependencies. A CI definition is not
proof of compatibility until those jobs pass; details are in
[`ci/README.md`](ci/README.md).

Passing them checks software behaviour; it does not establish manuscript
reproducibility or replace an independent clean-clone analysis with study data.

For every LL.4 fit, the lower/upper asymptotes are fixed to the observed
response minimum/maximum of the fitted rows; only slope and ED50 are estimated.
In the qPCR workflows, these bounds use the bulk RNA-seq gRNA means. Failures
are recorded rather than represented as successful curves.

Read [`METHODS.md`](METHODS.md) before interpreting AIC or fitted curves:

- Candidate comparisons require the same observations, response scale, and
  compatible likelihood definitions; spline information criteria are approximate.
- AIC preference is descriptive model support, not proof of saturation,
  causality, or a biological mechanism. Eleven gRNA means provide limited
  support for nonlinear fitting even with fixed asymptotes.
- The nested ANOVA compares linear with quadratic models, not LL.4.
- gRNA independence must be justified from the design. Repeated measurements,
  shared controls, or batch structure may require a different inferential model.
- SE bars, bulk-minus-qPCR percentage-point differences, and leave-one-gRNA-out
  sensitivity summaries are descriptive, not independent validation.

## Repository guide and publication status

- `R/`: reusable model, plotting, and CLI/provenance functions.
- `analysis/`: command-line analysis runners.
- `scripts/`: deterministic synthetic fixture generators.
- `tests/`: model and CLI test suites.
- [`METHODS.md`](METHODS.md): equations, normalization, assumptions, and limitations.
- [`QPCR_BULK_RNA_WORKFLOW.md`](QPCR_BULK_RNA_WORKFLOW.md): detailed qPCR workflow.
- [`REPRODUCIBILITY_CHECKLIST.md`](REPRODUCIBILITY_CHECKLIST.md): remaining manuscript
  reproduction and preservation requirements.
- [`CODE_AVAILABILITY.md`](CODE_AVAILABILITY.md): availability-statement template;
  placeholders must be completed before submission.

A final release should connect the exact source commit, approved inputs or
accessions, environment specification, execution commands, reference outputs,
and validation report. A Zenodo archive/DOI is not supplied by this README.
Once a release is archived, cite its **version-specific DOI** for the analysis
used in a manuscript; citation/authorship metadata remain to be finalized.

## License

Code and documentation are distributed under the [MIT License](LICENSE),
copyright 2026 Huangfu Lab. Study datasets require their own applicable
permissions and access terms.
