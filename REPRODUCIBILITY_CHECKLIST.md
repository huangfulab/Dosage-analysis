# Publication reproducibility checklist

## Present in the repository

- [x] Public-code license
- [x] Standalone functions with explicit inputs and structured outputs
- [x] No absolute paths, private workspace loads, or hidden global objects
- [x] Deterministic synthetic example
- [x] One-command example execution
- [x] Automated tests for all three model families
- [x] Input validation and explicit control normalization
- [x] Model equations, likelihood definition, and statistical cautions
- [x] Machine-readable metrics, coefficients, predictions, and fitted objects
- [x] Git commit, dirty-worktree state, run parameters, and R session capture

## Completed public code preservation

- [x] Public GitHub repository with versioned releases v1.0.0 and v1.0.1
- [x] Permanent source archives: v1.0.0 DOI `10.5281/zenodo.23620777` and
      v1.0.1 DOI `10.5281/zenodo.23622499`
- [x] Release/tag-to-commit mapping and both archive file inventories checked
      against Git on 2026-10-27
- [x] Software citation metadata in `CITATION.cff` for archived v1.0.1
- [x] Replace availability-template placeholders with verified information
      and explicit unresolved items in `CODE_AVAILABILITY.md`
- [x] Confirm public GEO GSE349847 (44-sample neuroectoderm bulk RNA-seq and
      K2LLPS expression matrix) and its parent SuperSeries GSE349851

These completed items establish public code/data availability within the
stated scope. They do not confirm the exact manuscript execution version or
complete manuscript reproducibility. The v1.0.1 release notes report local
study-data comparisons, but their records are not in the public archive and
were not independently revalidated in this documentation audit.

## Manuscript reproduction details still to verify or document

- [ ] Deposit or provide an accession for every study input used by the code
      (including the guide/qPCR workbook; the GEO series alone does not prove
      complete input coverage)
- [ ] Match the final manuscript run records and source fingerprints to an
      exact execution commit; distinguish it from release and import commits
- [ ] Add the raw-data-to-analysis-CSV preparation workflow
- [ ] Document experimental units, replicates, exclusions, and missing values
- [ ] Verify whether guide-level observations are independent; use an
      appropriate hierarchical model if the design requires one
- [ ] Add an exact figure/table-to-command-to-output map
- [ ] State all final parameter choices rather than relying on defaults
- [ ] Add expected checksums or reference summaries for manuscript outputs
- [ ] Run the workflow from a clean clone on an independent machine
- [ ] Confirm the manuscript's final code citation against its execution
      provenance and reconcile it with the Key Resources Table

## Scope of the remaining checks

Uploading functions to GitHub is useful but does not by itself reproduce a
paper. Editors or reviewers can reasonably ask for the data, preprocessing,
software versions, exact execution commands, final parameter settings, and the
link from code outputs to reported figures and conclusions. A mutable GitHub
branch also does not identify the exact execution version. A release, commit,
and permanent archive identify a source snapshot; connecting that snapshot
to the paper requires run-level evidence and figure/table mappings.

The statistical review may additionally question whether the sample size
supports four-parameter curves, whether AIC candidates share the same response
and likelihood definition, and whether gRNA-level observations are independent.
These limitations must be addressed in the analysis design and manuscript, not
only in code comments.
