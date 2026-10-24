#!/usr/bin/env Rscript

# Run from the repository root. Load only the adapter, never the study runner.
expressions = as.list(parse("studies/neuroectoderm_bulk_rnaseq/06_gsva_neuroectoderm_pca.R"))
adapter = Filter(function(expr) {
  is.call(expr) && identical(expr[[1L]], as.name("=")) &&
    identical(expr[[2L]], as.name("call_gsva_compat"))
}, expressions)
stopifnot(length(adapter) == 1L)
eval(adapter[[1L]])

expression = matrix(seq_len(24L), nrow=6L,
  dimnames=list(paste0("g", 1:6), paste0("s", 1:4)))
gene_sets = list(example=c("g1", "g3", "g5"))
legacy = function(expr, gset.idx.list, method, min.sz, max.sz, parallel.sz,
    verbose, BPPARAM, tau, ssgsea.norm=NULL, kcdf=NULL, mx.diff=NULL, abs.ranking=NULL) {
  stopifnot(identical(expr, expression), identical(gset.idx.list, gene_sets),
    min.sz == 1, is.infinite(max.sz), parallel.sz == 1L,
    identical(BPPARAM, "serial-test"), identical(verbose, F))
  if (identical(method, "ssgsea")) {
    stopifnot(tau == 0.25, identical(ssgsea.norm, T), is.null(kcdf))
  } else {
    stopifnot(tau == 1, identical(kcdf, "Gaussian"), identical(mx.diff, T), identical(abs.ranking, F))
  }
  method
}
for (method in c("ssgsea", "gsva")) {
  stopifnot(identical(call_gsva_compat(expression, gene_sets, method, legacy,
    verbose=F, BPPARAM="serial-test"), method))
}

# Early parameter-object releases have no constructor-level verbose argument.
ssgsea_param = function(exprData, geneSets, minSize, maxSize, alpha, normalize) {
  stopifnot(identical(exprData, expression), identical(geneSets, gene_sets),
    minSize == 1, is.infinite(maxSize), alpha == 0.25, identical(normalize, T))
  "ssgsea-parameter"
}
gsva_param = function(exprData, geneSets, minSize, maxSize, kcdf, tau, maxDiff, absRanking) {
  stopifnot(identical(exprData, expression), identical(geneSets, gene_sets),
    minSize == 1, is.infinite(maxSize), identical(kcdf, "Gaussian"), tau == 1,
    identical(maxDiff, T), identical(absRanking, F))
  "gsva-parameter"
}
modern = function(param, verbose, BPPARAM) {
  stopifnot(identical(verbose, F), identical(BPPARAM, "serial-test"))
  param
}
stopifnot(
  identical(call_gsva_compat(expression, gene_sets, "ssgsea", modern, ssgsea_param,
    verbose=F, BPPARAM="serial-test"), "ssgsea-parameter"),
  identical(call_gsva_compat(expression, gene_sets, "gsva", modern, gsva_param,
    verbose=F, BPPARAM="serial-test"), "gsva-parameter")
)
newer_param = function(exprData, geneSets, minSize, maxSize, alpha, normalize, verbose) {
  stopifnot(identical(verbose, F))
  ssgsea_param(exprData, geneSets, minSize, maxSize, alpha, normalize)
}
stopifnot(identical(call_gsva_compat(expression, gene_sets, "ssgsea", modern, newer_param,
  verbose=F, BPPARAM="serial-test"), "ssgsea-parameter"))
# An analysis failure must propagate, not silently switch APIs or algorithms.
failed = tryCatch(call_gsva_compat(expression, gene_sets, "ssgsea",
  function(...) stop("expected scoring failure"), ssgsea_param,
  verbose=F, BPPARAM="serial-test"), error=identity)
stopifnot(inherits(failed, "error"), identical(conditionMessage(failed), "expected scoring failure"))

# Check the installed package on every matrix entry. Old R entries run the real
# legacy implementation; mocks above only test dispatch/options on newer GSVA.
stopifnot(requireNamespace("GSVA", quietly=T), requireNamespace("BiocParallel", quietly=T))
set.seed(20261025)
expression = matrix(rnorm(100L * 12L), nrow=100L,
  dimnames=list(paste0("g", 1:100), paste0("s", 1:12)))
gene_sets = list(example=paste0("g", 1:20))
exports = getNamespaceExports("GSVA")
for (method in c("ssgsea", "gsva")) {
  constructor_name = if (method == "ssgsea") "ssgseaParam" else "gsvaParam"
  constructor = if (constructor_name %in% exports) getExportedValue("GSVA", constructor_name) else NULL
  observed = call_gsva_compat(expression, gene_sets, method, GSVA::gsva, constructor, verbose=F)
  stopifnot(is.matrix(observed), identical(dim(observed), c(1L, 12L)), all(is.finite(observed)),
    identical(colnames(observed), colnames(expression)), identical(rownames(observed), "example"))
  if (is.null(constructor)) {
    if (method == "ssgsea") {
      expected = GSVA::gsva(expression, gene_sets, method="ssgsea", min.sz=1, max.sz=Inf,
        tau=0.25, ssgsea.norm=T, parallel.sz=1L, verbose=F, BPPARAM=BiocParallel::SerialParam())
    } else {
      expected = GSVA::gsva(expression, gene_sets, method="gsva", min.sz=1, max.sz=Inf,
        kcdf="Gaussian", tau=1, mx.diff=T, abs.ranking=F, parallel.sz=1L,
        verbose=F, BPPARAM=BiocParallel::SerialParam())
    }
  } else {
    arguments = list(exprData=expression, geneSets=gene_sets, minSize=1, maxSize=Inf)
    if (method == "ssgsea") {
      arguments$alpha = 0.25
      arguments$normalize = T
    } else {
      arguments$kcdf = "Gaussian"
      arguments$tau = 1
      arguments$maxDiff = T
      arguments$absRanking = F
    }
    if ("verbose" %in% names(formals(constructor))) arguments$verbose = F
    expected = GSVA::gsva(do.call(constructor, arguments), verbose=F, BPPARAM=BiocParallel::SerialParam())
  }
  stopifnot(isTRUE(all.equal(observed, expected, tolerance=1e-12)))
  message("PASS: ", method, " ", if (is.null(constructor)) "legacy" else "parameter-object", " API")
}
message("GSVA compatibility dispatch, settings, failure propagation, and installed-package tests passed.")
devtools::session_info(pkgs="attached")
