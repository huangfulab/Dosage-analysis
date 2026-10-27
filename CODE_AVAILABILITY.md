# Data and code availability

## Public code and archived releases

Dosage-response, qPCR--bulk RNA-seq comparison, and neuroectoderm bulk RNA-seq
analysis code is publicly available at
https://github.com/huangfulab/Dosage-analysis under the MIT License.
The existing source archives are:

| GitHub release | Zenodo version DOI | Tagged source commit |
| --- | --- | --- |
| [v1.0.0](https://github.com/huangfulab/Dosage-analysis/releases/tag/v1.0.0) | [10.5281/zenodo.23620777](https://doi.org/10.5281/zenodo.23620777) | `910af968793259a679718fc1c89177a0b38d4bfb` |
| [v1.0.1](https://github.com/huangfulab/Dosage-analysis/releases/tag/v1.0.1) | [10.5281/zenodo.23622499](https://doi.org/10.5281/zenodo.23622499) | `b53bd385c78c94c0152511e5ea3848c80a00fd57` |

Both releases were published on 2026-10-24. Their Zenodo records share the
concept DOI [10.5281/zenodo.23620776](https://doi.org/10.5281/zenodo.23620776).
The concept DOI identifies the version family; the two version DOIs identify
fixed source archives. No additional deposit is needed to cite these releases.

On 2026-10-27, both Zenodo ZIP checksums matched their record metadata, and
all 46 files in each ZIP matched the Git blob hashes at its tagged commit.
The annotated `v1.0.1` tag object `04fa911e53ad025a004fd70381f46f3afb0037b7`
resolves to the source commit shown above; the tag-object hash is not the
analysis commit. The [release diff](https://github.com/huangfulab/Dosage-analysis/compare/v1.0.0...v1.0.1)
contains only CI configuration/tests and documentation changes. Analysis
scripts are identical between these two releases.

`CITATION.cff` describes v1.0.1 and its version-specific DOI. It was added
after that release and is not included in either existing source archive.
The current `main` branch also contains subsequent CI changes and must not be
treated as the content of either archive.

## Manuscript citation and execution provenance

The selected public code archive for sharing and manuscript code availability
is **v1.0.1**, DOI **10.5281/zenodo.23622499**, at commit
`b53bd385c78c94c0152511e5ea3848c80a00fd57`. Use this version-specific DOI
in the manuscript's code availability statement and corresponding software
citation. DOI `10.5281/zenodo.23620777` remains a historical reference to
v1.0.0; it does not need to be withdrawn or redeposited.

**The exact source commit that generated the final manuscript figures and
tables has not been established by this archive audit.** Neither release
should be described as that execution version without matching original
run manifests, source/input fingerprints, parameter settings, and outputs.
The study README records import provenance from `hyein-cbio/K2LLPS_DESeq2`
at commit `404fea0c496b6eec66f7586046c4b9064f8d0560`; that import does not
by itself identify the final manuscript execution commit.

A statement limited to verified code availability is:

> Analysis source code is publicly available at
> https://github.com/huangfulab/Dosage-analysis. Source release v1.0.1
> (commit b53bd385c78c94c0152511e5ea3848c80a00fd57) is archived at
> https://doi.org/10.5281/zenodo.23622499.

This statement describes code availability. It does not assert that an
archived release reproduces all manuscript results.

## Public study data and input coverage

GEO [GSE349847](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE349847)
is public and contains 44 neuroectoderm bulk RNA-seq (3'-seq) samples and
the supplementary file `GSE349847_K2LLPS-expression-matrix.tsv.gz`.
It is a SubSeries of public SuperSeries
[GSE349851](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE349851),
which also links the Perturb-seq series GSE283614.

The GEO records were checked on 2026-10-27. This confirms public accessions
and the listed expression matrix, not byte-level identity with the local
files used for the manuscript. The source archives do not contain study
expression matrices, the `260320_OCT4_qPCR.xlsx` guide/qPCR workbook, or
generated study outputs. A public location for that workbook and a complete
inventory of all inputs used across the manuscript remain unverified.

Workflow instructions are in [`QPCR_BULK_RNA_WORKFLOW.md`](QPCR_BULK_RNA_WORKFLOW.md)
and [`studies/neuroectoderm_bulk_rnaseq/README.md`](studies/neuroectoderm_bulk_rnaseq/README.md).
These describe workflow commands and expected files, but do not provide a
complete final-manuscript figure/table-to-input-to-command-to-output map.

## Reproduction status

The [v1.0.1 release notes](https://github.com/huangfulab/Dosage-analysis/releases/tag/v1.0.1)
report local execution of the original 44-sample workflow and numerical
comparisons with preserved analysis results. Those local validation records
and outputs are not included in the public source release, and this audit
did not rerun or independently verify those numerical comparisons.

Public source preservation, reported local validation, CI/synthetic test
coverage, and independent reproduction of manuscript results are separate
claims. Remaining work includes confirming the execution commit and full
input coverage, documenting raw-data preprocessing and figure-specific
parameters, preserving reference outputs and a fixed environment, and
testing the complete manuscript workflow from a clean clone on an
independent machine. See [`REPRODUCIBILITY_CHECKLIST.md`](REPRODUCIBILITY_CHECKLIST.md).
