.PHONY: test example qpcr-bulk-example qpcr-bulk-deg-example

test:
	Rscript tests/test_models.R
	Rscript tests/test_cli.R
	Rscript tests/test_qpcr_bulk_rna_models.R
	Rscript tests/test_qpcr_bulk_rna_cli.R
	Rscript tests/test_qpcr_bulk_rna_deg_batch.R
	Rscript tests/test_qpcr_bulk_rna_deg_batch_cli.R

example:
	Rscript scripts/generate_example_data.R data/example_dosage.csv
	Rscript analysis/run_dosage_analysis.R --input data/example_dosage.csv --dose-col dose --response-cols response_hill,response_quadratic --group-col group --control-label Controls --normalize-dose true --normalize-response true --output-dir results/example

qpcr-bulk-example:
	Rscript scripts/generate_qpcr_bulk_rna_example.R data/example_qpcr_bulk_rna.csv
	Rscript analysis/run_qpcr_bulk_rna_comparison.R --input data/example_qpcr_bulk_rna.csv --grna-col grna --qpcr-col qpcr_percent --bulk-rna-col bulk_rna_percent --output-dir results/qpcr_bulk_rna_example

qpcr-bulk-deg-example:
	Rscript scripts/generate_qpcr_bulk_rna_deg_example.R data/example_qpcr_bulk_rna_degs.csv
	Rscript analysis/run_qpcr_bulk_rna_deg_batch.R --input data/example_qpcr_bulk_rna_degs.csv --sample-id-col sample_id --grna-col grna --qpcr-col qpcr_percent --gene-col gene --bulk-rna-col bulk_rna_percent --output-dir results/qpcr_bulk_rna_degs_example
