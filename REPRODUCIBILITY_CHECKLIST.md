# Publication reproducibility checklist

## Present in this draft

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

## Required before manuscript submission

- [ ] Deposit or provide an accession for every study input used by the code
- [ ] Add the raw-data-to-analysis-CSV preparation workflow
- [ ] Document experimental units, replicates, exclusions, and missing values
- [ ] Verify whether guide-level observations are independent; use an
      appropriate hierarchical model if the design requires one
- [ ] Add an exact figure/table-to-command-to-output map
- [ ] State all final parameter choices rather than relying on defaults
- [ ] Add expected checksums or reference summaries for manuscript outputs
- [ ] Run the workflow from a clean clone on an independent machine
- [ ] Create a versioned release and permanent archive DOI
- [ ] Fill every placeholder in `CODE_AVAILABILITY.md`
- [ ] Add citation metadata and final authorship information

## Likely editorial or reviewer objections if left incomplete

Uploading functions to GitHub is useful but does not by itself reproduce a
paper. Editors or reviewers can reasonably ask for the data, preprocessing,
software versions, exact execution commands, final parameter settings, and the
link from code outputs to reported figures and conclusions. A mutable GitHub
branch also does not identify the exact reviewed version; a release, commit,
and permanent archive make that version unambiguous.

The statistical review may additionally question whether the sample size
supports four-parameter curves, whether AIC candidates share the same response
and likelihood definition, and whether gRNA-level observations are independent.
These limitations must be addressed in the analysis design and manuscript, not
only in code comments.
