# R compatibility CI

`.github/workflows/r-compatibility.yaml` checks **every R 4.x minor release
series from 4.0 through 4.6**, using the latest patch release in each series on
Ubuntu 22.04. These are seven separate jobs, not a single R-major-version job.
The workflow runs on pushes, pull requests, manual dispatch, and a weekly
schedule. It does not publish releases or analysis results.

## Assumed minimum supported dependencies

These are **proposed minimum support baselines**, not claims that the code has
already passed on all seven environments or that no still-older package could
work. Each job must install the exact direct-package versions in
[`minimum-packages.csv`](minimum-packages.csv), load them, check their hard
dependency constraints, and pass the tests before the corresponding baseline
can be considered validated. A higher installed version is a failure, not an
acceptable substitute for a minimum-version test.

Common direct-package floors:

| Package | Assumed minimum |
| --- | --- |
| drc | 3.0-1 |
| ggplot2 | 3.4.0 (the study plots use `linewidth`) |
| ggrepel | 0.9.2 |
| openxlsx | 4.2.5.1 |
| pheatmap | 1.0.12 |
| devtools | 2.4.5 |

R/Bioconductor-specific floors:

| R series | Bioconductor | DESeq2 | GSVA | matrixStats | CRAN snapshot |
| --- | --- | --- | --- | --- | --- |
| 4.0 | 3.11 | 1.28.1 | 1.36.3 | 0.63.0 | 2022-11-08 |
| 4.1 | 3.13 | 1.32.0 | 1.40.1 | 0.63.0 | 2022-11-08 |
| 4.2 | 3.15 | 1.36.0 | 1.44.5 | 0.63.0 | 2022-11-08 |
| 4.3 | 3.17 | 1.40.2 | 1.48.3 | 1.0.0 | 2023-09-15 |
| 4.4 | 3.19 | 1.44.0 | 1.52.3 | 1.0.0 | 2024-06-15 |
| 4.5 | 3.21 | 1.48.2 | 2.2.1 | 1.4.1 | 2025-06-15 |
| 4.6 | 3.23 | 1.52.0 | 2.6.14 | 1.4.1 | 2026-07-01 |

The Bioconductor release/package series follow the earliest Bioconductor
release for each R series, using explicit maintenance patch versions. These
conservative floors are assumptions to test, **not the absolute first version
ever released**. The matrixStats floors account for the selected
MatrixGenerics versions; old Bioconductor code must not be combined silently
with a much newer matrixStats that has removed legacy `useNames` behaviour.

Only runtime/build dependencies (`Depends`, `Imports`, `LinkingTo`) are
installed, not the entire `Suggests` graph. The installer uses the selected
Bioconductor repositories and a dated Posit CRAN snapshot for transitive
dependencies, then asserts all nine direct pins and dependency constraints.
Transitive dependencies are **not** claimed to be their mathematical minimum
versions. All actually installed versions, repositories, R session details,
and test logs are uploaded as CI artifacts even when a job fails. This is not
a complete transitive lockfile or a manuscript reproduction environment.

Installation occurs only with `Rscript ci/dependencies.R --install` inside a
disposable GitHub Actions library. Local `--plan` is read-only; the ordinary
analysis runners and tests never install packages. There is no cross-job
library cache that could mask a missing dependency or replace a floor with a
previously installed newer version.

## Tests required in every R job

1. Exact-version and hard-dependency verification, plus CI utility tests.
2. `make test`: all six generic model and CLI test suites.
3. `Rscript tests/test_gsva_compat.R`: legacy and parameter-object argument
   mapping, early constructors without `verbose`, error propagation, and
   comparison with native calls to the installed GSVA implementation.
4. `Rscript tests/test_neuroectoderm_workflow.R`: the complete synthetic
   44-sample study workflow, including actual GSVA/ssGSEA scoring and plots.

Old-R jobs use the **real old GSVA package**, not mocks in place of enrichment
scoring. Mocks additionally exercise dispatch on a modern development machine.
The compatibility adapter changes API spelling only: ssGSEA uses alpha/tau
0.25 and normalization; GSVA uses a Gaussian kernel, tau 1, maxDiff TRUE, and
absRanking FALSE. Input VST, gene sets, PCA, model design, and cutoffs are
unchanged. The selected API and GSVA version are recorded in the study summary.
Upstream GSVA changes can still affect numerical results between releases;
the tests do not claim cross-version bit identity or identical study findings.

The optional `glmGamPoi` count strategy is not part of these baseline jobs;
the complete study test uses the default rounded-count strategy. No private
study data or existing local result files are required or uploaded.

## Coverage maintenance

[`r-versions.csv`](r-versions.csv) and `minimum-packages.csv` are the source of
truth. The matrix generator queries the R-hub stable-release API and **fails**
if a newly released R series has no declared baseline, instead of silently
omitting it. Extend both files, their tests, and this table for a new series.
The release lookup also fails on network errors; such a failure is not a green
compatibility result. Weekly/manual workflows become available according to
GitHub's default-branch workflow rules after the configuration is merged.

Read-only/offline local checks from the repository root:

```bash
python3 -m unittest discover -s ci -p 'test_*.py' -v
python3 ci/r_matrix.py --release-version 4.6.1
Rscript ci/dependencies.R --plan
Rscript tests/test_ci_dependencies.R
Rscript tests/test_gsva_compat.R
```

`--verify` is intended for the exact baseline environment; it should reject a
normal development library containing newer versions. Inspect actual GitHub
Actions results before saying that all R series and their floors pass.

## Upstream references

- [R setup and minor-series resolution](https://github.com/r-lib/actions/tree/v2/setup-r)
- [R/Bioconductor release mapping](https://bioconductor.org/about/release-announcements/)
- [Legacy GSVA matrix/list API source](https://github.com/bioc/GSVA/blob/RELEASE_3_12/R/gsva.R)
- [Early parameter-object API source](https://github.com/bioc/GSVA/blob/RELEASE_3_18/R/ssgseaParam.R)
- [Posit package snapshots](https://packagemanager.posit.co/)
