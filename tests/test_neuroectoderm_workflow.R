#!/usr/bin/env Rscript

# Run from the Dosage-analysis repository root. These are synthetic software
# checks, not study results. All generated inputs/results live outside the repo.
main = function() {
  repository_dir = normalizePath(getwd(), mustWork=T)
  study_dir = file.path(repository_dir, "studies", "neuroectoderm_bulk_rnaseq")
  required_packages = c(
    "DESeq2", "openxlsx", "devtools", "ggplot2", "ggrepel", "pheatmap",
    "matrixStats", "GSVA", "drc"
  )
  missing = required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly=T)]
  if (length(missing) > 0L) stop("Missing test dependencies: ", paste(missing, collapse=", "))

  # Load only the label helper, not the runner's top-level analysis statements.
  label_expression = Filter(function(expr) {
    is.call(expr) && identical(expr[[1L]], as.name("=")) &&
      identical(expr[[2L]], as.name("direction_text"))
  }, as.list(parse(file.path(study_dir, "05_gene_vs_qpcr_scatter.R"))))
  stopifnot(length(label_expression) == 1L)
  label_env = new.env(parent=baseenv())
  eval(label_expression[[1L]], envir=label_env)
  stopifnot(
    identical(label_env$direction_text(NULL), "single-gene"),
    identical(label_env$direction_text("all"), "all genes"),
    identical(label_env$direction_text("not_deg"), "not DEG"),
    identical(label_env$direction_text("up"), "UP DEG")
  )

  # Fingerprint existing study fixtures/results: this test must never change them.
  protected_files = c(
    list.files(file.path(repository_dir, "data"), full.names=T, recursive=T, all.files=T),
    list.files(file.path(repository_dir, "results"), full.names=T, recursive=T, all.files=T),
    file.path(study_dir, "sample_metadata_draft.csv")
  )
  protected_hashes = tools::md5sum(protected_files)
  test_dir = tempfile("neuroectoderm-workflow-", tmpdir=Sys.getenv("TMPDIR", unset="/tmp"))
  dir.create(test_dir)
  message("Synthetic validation workspace: ", test_dir)
  successful = F
  on.exit({
    setwd(repository_dir)
    if (successful) unlink(test_dir, recursive=T) else message("Failure artifacts retained: ", test_dir)
  }, add=T)
  bundle_dir = file.path(test_dir, "relocated-study")
  dir.create(bundle_dir)
  scripts = list.files(study_dir, pattern="^[0-6][0-6]_.*[.]R$", full.names=T)
  stopifnot(length(scripts) == 7L)
  stopifnot(all(file.copy(c(scripts, file.path(study_dir, "sample_metadata_draft.csv")), bundle_dir)))
  input_dir = file.path(bundle_dir, "inputs")
  output_dir = file.path(bundle_dir, "outputs")
  dir.create(input_dir)
  metadata = read.csv(file.path(bundle_dir, "sample_metadata_draft.csv"))
  expected_guides = sprintf("guide%02d", c(1:6, 14:15, 18:20))
  stopifnot(
    NROW(metadata) == 44L,
    !anyDuplicated(metadata$count_column), !anyDuplicated(metadata$sample_id),
    identical(metadata$sample_id, paste0(metadata$guide, "_rep_", metadata$replicate)),
    setequal(metadata$guide, expected_guides), setequal(metadata$replicate, 2:5),
    all(table(metadata$guide, metadata$replicate) == 1L)
  )
  pooled_design = model.matrix(~factor(replicate) + factor(guide %in% expected_guides[1:8]), metadata)
  guide_design = model.matrix(~factor(replicate) + factor(ifelse(guide %in% expected_guides[1:8], guide, "control")), metadata)
  stopifnot(qr(pooled_design)$rank == 5L, qr(guide_design)$rank == 12L)

  # A strongly specified interaction fixture preserves the study verifier's
  # original sole-GATAD2B expectation without disabling or bypassing that check.
  set.seed(20261025)
  gene_count = 1201L
  gene_ids = c(sprintf("SYNTH%05d", seq_len(gene_count - 1L)), "ENSG00000143614")
  neural_symbols = c(
    "ADGRG1", "ASCL1", "BMI1", "CD24", "CD44", "CD81", "CXCR4", "DCX",
    "HES5", "ID1", "ITGA6", "ITGB1", "MAPT", "MCAM", "NES", "NGFR",
    "NOTCH1", "PAX6", "PROM1", "SLITRK6", "SOX1", "SOX2", "SOX9", "VIM"
  )
  gene_names = c(neural_symbols, "DUPLICATE", "DUPLICATE", sprintf("TOY%04d", 27:1200), "GATAD2B")
  baseline = exp(runif(gene_count, log(80), log(1500)))
  # Low-ranked neural markers deliberately exercise valid negative ssGSEA scores.
  baseline[seq_along(neural_symbols)] = runif(length(neural_symbols), 2, 10)
  lfc = rep(0, gene_count)
  lfc[1:70] = 2.5
  lfc[71:140] = -2.5
  is_target = metadata$guide %in% expected_guides[1:8]
  mean_counts = outer(baseline, rep(1, 44L)) * 2 ^ outer(lfc, as.integer(is_target))
  block_effect = c("2"=0.9, "3"=1, "4"=1.1, "5"=1.05)
  mean_counts = sweep(mean_counts, 2L, block_effect[as.character(metadata$replicate)], "*")
  interaction_effect = c("2"=8, "3"=0.125, "4"=8, "5"=0.125)
  mean_counts[gene_count, is_target] = mean_counts[gene_count, is_target] * interaction_effect[as.character(metadata$replicate[is_target])]
  dispersion = 0.015 + 2 / baseline
  counts = matrix(rnbinom(length(mean_counts), mu=as.numeric(mean_counts), size=rep(1 / dispersion, times=44L)), nrow=gene_count)
  # Balance non-interaction genes within each group/block to avoid chance
  # interaction discoveries in a regression fixture with one specified hit.
  # Guide-level/sample noise remains for the guide-specific dispersion model.
  regular = seq_len(gene_count - 1L)
  for (block in 2:5) {
    for (target_group in c(F, T)) {
      columns = which(metadata$replicate == block & is_target == target_group)
      observed_sum = rowSums(counts[regular, columns, drop=F])
      desired_sum = rowSums(mean_counts[regular, columns, drop=F])
      values = counts[regular, columns, drop=F]
      empty = observed_sum == 0
      values[empty, ] = mean_counts[regular[empty], columns, drop=F]
      observed_sum[empty] = desired_sum[empty]
      counts[regular, columns] = sweep(values, 1L, desired_sum / observed_sum, "*")
    }
  }
  # Exercise fractional count rounding while retaining an equivalent CSV/TSV.
  counts = counts + 0.25
  annotation = data.frame(gene_id=gene_ids, gene_name=gene_names, gene_biotype="protein_coding")
  annotation$gene_biotype[[26L]] = "lncRNA"
  sox1_index = which(gene_names == "SOX1")
  annotation$gene_biotype[sox1_index] = NA_character_
  annotation$gene_name[[300L]] = 'TOY, "quoted" marker'
  tsv = cbind(annotation, as.data.frame(counts), as.data.frame(counts / 10))
  names(tsv) = c(names(annotation), metadata$count_column, sub("_count$", "_cpm", metadata$count_column))
  csv = tsv
  names(csv) = c(names(annotation), paste0(metadata$sample_id, "_count"), paste0(metadata$sample_id, "_cpm"))
  write.table(tsv, file.path(input_dir, "K2LLPS-expression-matrix.tsv"), sep="\t", quote=T, qmethod="double", row.names=F)
  write.csv(csv, file.path(input_dir, "K2LLPS-expression-matrix.csv"), row.names=F)
  workbook = data.frame(
    HZ_gRNA_name=as.integer(sub("guide", "", expected_guides)),
    JY_gRNA_name=paste0("  synthetic_", expected_guides, "  ")
  )
  workbook$JY_gRNA_name[[1L]] = '  synthetic, "guide01"  '
  qpcr_mean = c(20, 30, 40, 50, 60, 70, 80, 90, 98, 100, 102)
  for (block in 2:5) workbook[[paste0("rep", block, "_OCT4")]] = qpcr_mean + (block - 3.5)
  workbook$rep5_OCT4[workbook$HZ_gRNA_name == 4L] = NA_real_
  openxlsx::write.xlsx(list(Sheet1=workbook), file.path(input_dir, "260320_OCT4_qPCR.xlsx"))

  # Deliberately execute from a different cwd to test script-relative defaults.
  setwd(test_dir)
  run_index = 0L
  run = function(script, arguments=character(), env=character(), expected_status=0L, expected_text=NULL) {
    run_index <<- run_index + 1L
    log_path = file.path(test_dir, sprintf("%02d-%s.log", run_index, sub("[.]R$", "", script)))
    defaults = c(K2LLPS_COUNT_STRATEGY="round", K2LLPS_MIN_TOTAL_COUNT="10", K2LLPS_ALPHA="0.05", K2LLPS_DEG_LFC="1")
    defaults[names(env)] = env
    status = suppressWarnings(system2(
      file.path(R.home("bin"), "Rscript"),
      shQuote(c(file.path(bundle_dir, script), arguments)),
      stdout=log_path, stderr=log_path,
      env=paste(names(defaults), defaults, sep="=")
    ))
    log = readLines(log_path, warn=F)
    if (!identical(as.integer(status), as.integer(expected_status)) ||
        (!is.null(expected_text) && !any(grepl(expected_text, log, fixed=T)))) {
      message(paste(tail(log, 60L), collapse="\n"))
      stop("Unexpected validation outcome: ", script, "; log: ", log_path)
    }
    message("PASS: ", script, if (length(arguments)) paste0(" ", paste(arguments, collapse=" ")) else " [default paths]")
    invisible(log)
  }
  run("00_verify_expression_files.R")
  run("01_prepare_deseq2.R", env=c(K2LLPS_DEG_LFC="2"))
  run("02_verify_pooled_deseq2_outputs.R", expected_text="PASS: verified 44 samples")
  pooled_metadata = read.csv(file.path(output_dir, "analysis-sample-metadata.csv"))
  stopifnot(all(pooled_metadata$jy_guide_label[pooled_metadata$guide == "guide01"] == 'synthetic, "guide01"'))
  run("03_run_guide_specific_deseq2.R")
  run("04_verify_guide_specific_outputs.R", expected_text="PASS: verified 44 samples")
  no_hits_dir = file.path(test_dir, "no-significant-guides")
  run("03_run_guide_specific_deseq2.R", c(output_dir, no_hits_dir), env=c(K2LLPS_ALPHA="1e-300"))
  run("04_verify_guide_specific_outputs.R", c(output_dir, no_hits_dir), expected_text="PASS: verified 44 samples")
  no_hits = read.csv(file.path(no_hits_dir, "concordance-across-guides.csv"))
  stopifnot(all(no_hits$guides_significant_fdr == 0L))
  run("05_gene_vs_qpcr_scatter.R", c(output_dir, "--batch-deg", "2", "--csv-only"))
  run("05_gene_vs_qpcr_scatter.R", c(output_dir, "--batch-deg-pooled", "2", "--csv-only"))
  run("05_gene_vs_qpcr_scatter.R", c(output_dir, "--batch-all", "2"))
  run("05_gene_vs_qpcr_scatter.R", c(output_dir, "SOX1"))
  run("06_gsva_neuroectoderm_pca.R")
  mapping = read.csv(file.path(output_dir, "gsva-neuroectoderm-gene-set.csv"))
  stopifnot(identical(mapping$gene_id[mapping$gene_symbol == "SOX1"], gene_ids[sox1_index]))
  scores = read.csv(file.path(output_dir, "gsva-neuroectoderm-sample-scores.csv"))
  stopifnot(NROW(scores) == 44L, all(is.finite(scores$ssgsea)), any(scores$ssgsea < 0))
  run("03_run_guide_specific_deseq2.R", c(output_dir, bundle_dir), expected_status=1L,
    expected_text="must not replace or contain the pooled output directory")
  run("03_run_guide_specific_deseq2.R", c(output_dir, output_dir), expected_status=1L,
    expected_text="must not replace or contain the pooled output directory")
  run("03_run_guide_specific_deseq2.R", c(output_dir, file.path(output_dir, "..")), expected_status=1L,
    expected_text="must not replace or contain the pooled output directory")
  pooled_files = list.files(output_dir, recursive=F, full.names=T)
  pooled_hashes = tools::md5sum(pooled_files[!file.info(pooled_files)$isdir])
  missing_alias = file.path(bundle_dir, "new-component", "..", "outputs")
  run("03_run_guide_specific_deseq2.R", c(output_dir, missing_alias), expected_status=1L,
    expected_text="must not replace or contain the pooled output directory")
  stopifnot(identical(pooled_hashes, tools::md5sum(names(pooled_hashes))))
  # A regular-file destination forces publication to fail after staging.
  blocked_destination = file.path(test_dir, "blocked-guide-output")
  writeLines("sentinel", blocked_destination)
  run("03_run_guide_specific_deseq2.R", c(output_dir, blocked_destination), expected_status=1L,
    expected_text="Could not publish the staged guide-specific output directory")
  stopifnot(identical(readLines(blocked_destination), "sentinel"),
    !any(grepl("^.blocked-guide-output-staging-", list.files(test_dir, all.files=T))))

  external_pooled = file.path(test_dir, "external-pooled")
  dir.create(external_pooled)
  source_files = c(file.path(bundle_dir, basename(scripts)),
    file.path(bundle_dir, "sample_metadata_draft.csv"))
  source_hashes = tools::md5sum(source_files)
  run("03_run_guide_specific_deseq2.R", c(external_pooled, bundle_dir), expected_status=1L,
    expected_text="study source directory")
  stopifnot(identical(source_hashes, tools::md5sum(source_files)))

  # Regressions for tolerance and inconsistent expression files, isolated from
  # the pooled source outputs whose hashes are required to remain unchanged.
  annotation_csv = file.path(test_dir, "annotation-only.csv")
  annotation_tsv = file.path(test_dir, "annotation-only.tsv")
  write.csv(annotation[1:3, ], annotation_csv, row.names=F)
  write.table(annotation[1:3, ], annotation_tsv, sep="\t", quote=T, row.names=F)
  run("00_verify_expression_files.R", c(annotation_csv, annotation_tsv,
    file.path(test_dir, "annotation-only-check")), expected_status=1L, expected_text="PRACTICALLY_EQUAL: FALSE")
  small_csv = csv[1:3, 1:5]
  small_tsv = small_csv
  names(small_tsv)[4:5] = c("renamed_1_count", "renamed_2_count")
  small_tsv[[4L]][[1L]] = small_tsv[[4L]][[1L]] + 1e-7
  csv_path = file.path(test_dir, "tolerance.csv")
  tsv_path = file.path(test_dir, "tolerance.tsv")
  write.csv(small_csv, csv_path, row.names=F)
  write.table(small_tsv, tsv_path, sep="\t", quote=F, row.names=F)
  run("00_verify_expression_files.R", c(csv_path, tsv_path, file.path(test_dir, "tolerance-pass"), "0.000001"))
  run("00_verify_expression_files.R", c(csv_path, tsv_path, file.path(test_dir, "tolerance-fail"), "0"), expected_status=1L)

  # A union with no DOWN genes must still combine compatible zero-row schemas.
  # Work on a copy so the independently verified guide outputs stay untouched.
  one_sided = file.path(test_dir, "one-sided")
  dir.create(one_sided)
  stopifnot(all(file.copy(file.path(output_dir, c("DESeq2Obj-after-wald-test.rds", "analysis-sample-metadata.csv")), one_sided)))
  dir.create(file.path(one_sided, "guide-specific"))
  stopifnot(file.copy(file.path(output_dir, "guide-specific", "pooled-source-md5.csv"),
    file.path(one_sided, "guide-specific", "pooled-source-md5.csv")))
  for (guide in expected_guides[1:8]) {
    filename = sprintf("de-results-%s-vs-pooled-controls-complete.csv", guide)
    d = read.csv(file.path(output_dir, "guide-specific", filename))
    d$padj[d$log2FoldChange_raw <= 0 | is.na(d$log2FoldChange_raw)] = NA_real_
    write.csv(d, file.path(one_sided, "guide-specific", filename), row.names=F)
  }
  stale_down_pdf = file.path(one_sided, "scatter-deg-down-vs-qpcr-model-fits-preview2.pdf")
  writeLines("stale direction plot", stale_down_pdf)
  run("05_gene_vs_qpcr_scatter.R", c(one_sided, "--batch-deg", "2"))
  stopifnot(!file.exists(stale_down_pdf))
  up = read.csv(file.path(one_sided, "scatter-deg-up-vs-qpcr-model-aic-preview2.csv"))
  down = read.csv(file.path(one_sided, "scatter-deg-down-vs-qpcr-model-aic-preview2.csv"))
  stopifnot(NROW(up) > 0L, NROW(down) == 0L, identical(names(up), names(down)))

  # A padj exactly at the requested cutoff must be recoverable from complete
  # tables, even when the upstream significant-only export would omit that row.
  boundary_id = gene_ids[[30L]]
  for (guide in expected_guides[1:8]) {
    path = file.path(one_sided, "guide-specific", sprintf("de-results-%s-vs-pooled-controls-complete.csv", guide))
    d = read.csv(path)
    d$padj = NA_real_
    hit = d$gene_id == boundary_id
    d$padj[hit] = 0.1
    d$log2FoldChange_raw[hit] = 2
    d$log2FoldChange_shrunken[hit] = 1.5
    write.csv(d, path, row.names=F)
  }
  run("05_gene_vs_qpcr_scatter.R", c(one_sided, "--batch-deg", "--csv-only"), env=c(K2LLPS_ALPHA="0.1"))
  boundary = read.csv(file.path(one_sided, "scatter-deg-vs-qpcr-model-aic.csv"))
  stopifnot(identical(boundary$gene_id, boundary_id), boundary$n_guides_deg == 8L)
  # A new pooled object must not be paired with old stage-03 guide results.
  stale_path = file.path(one_sided, "DESeq2Obj-after-wald-test.rds")
  stale_dds = readRDS(stale_path)
  attr(stale_dds, "synthetic_revision") = "rerun-stage-01"
  saveRDS(stale_dds, stale_path)
  run("05_gene_vs_qpcr_scatter.R", c(one_sided, "--batch-deg", "--csv-only"),
    expected_status=1L, expected_text="Stale guide-specific outputs")
  stopifnot(identical(protected_hashes, tools::md5sum(protected_files)))
  message("PASS: existing repository data/results unchanged; synthetic execution and regressions completed.")
  successful = T
}

main()
sessionInfo()
