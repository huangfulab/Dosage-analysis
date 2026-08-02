# Input data

The repository does not yet contain study data.

`scripts/generate_example_data.R` creates `data/example_dosage.csv`, a
deterministic synthetic dataset used only to test the software. The generated
file is ignored by Git.

`scripts/generate_qpcr_bulk_rna_example.R` creates
`data/example_qpcr_bulk_rna.csv`, a deterministic 44-row synthetic dataset
for the qPCR--bulk RNA-seq workflow. It contains four sample-level rows for
each of eleven gRNA species and is also ignored by Git.

`scripts/generate_qpcr_bulk_rna_deg_example.R` creates
`data/example_qpcr_bulk_rna_degs.csv`, a deterministic long-form fixture with
three genes and the same 44 samples for each gene.

## Study-data schema

Prepare one CSV with:

- one row per declared analysis unit;
- one finite numeric dose column;
- one or more finite numeric response columns;
- an optional group column, required for control normalization; and
- stable identifiers for guide, sample, replicate, target, and batch whenever
  those fields exist.

The final manuscript archive must document how raw deposited data become this
table. Include all inclusion and exclusion rules, aggregation operations,
normalizations, and mappings from identifiers to experimental groups.

Do not add confidential, controlled-access, or personally identifying data to
Git. Deposit such data in an appropriate repository and document the accession
and access procedure.

## qPCR--bulk RNA-seq DEG input schema

The batch qPCR--bulk RNA-seq runner expects long-form data: one row per sample
and gene. The final fitted data contain one mean per gRNA species for every
gene.

| Field | Requirement |
| --- | --- |
| Sample ID column | Non-missing identifier; unique within each gene and identical across every selected gene |
| gRNA column | Non-missing species identifier; exactly four complete rows per species by default |
| qPCR column | Finite numeric normalized expression in percent of controls; every value must be greater than zero for LL.4 |
| Gene column | Non-missing DEG identifier; one row per sample within each gene |
| Bulk RNA-seq column | Finite numeric normalized expression in percent of controls |

The default batch runner deliberately requires 44 complete rows and eleven
gRNA species **for every gene**. Change `--expected-replicates` and
`--expected-points` only if the experimental design differs, and document that
change in the manuscript.
