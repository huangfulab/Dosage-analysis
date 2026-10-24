#!/usr/bin/env Rscript

# Run from the repository root. Check shared LL.4 capture behavior without
# changing the caller's redirected message stream, on success and fit failure.
source("R/dosage_models.R")
source("R/qpcr_bulk_rna_models.R")
study_expressions = as.list(parse("studies/neuroectoderm_bulk_rnaseq/05_gene_vs_qpcr_scatter.R"))
study_environment = new.env(parent=baseenv())
fit_expression = Filter(function(expr) {
  is.call(expr) && identical(expr[[1L]], as.name("=")) && identical(expr[[2L]], as.name("fit_ll4"))
}, study_expressions)
stopifnot(length(fit_expression) == 1L)
eval(fit_expression[[1L]], envir=study_environment)

check_message_capture = function(run, reuse_lower_connection=F) {
  previous = getConnection(sink.number(type="message"))
  path = tempfile()
  filler = if (reuse_lower_connection) file(tempfile(), open="wt") else NULL
  connection = file(path, open="wt")
  if (!is.null(filler)) {
    filler_path = summary(filler)$description
    close(filler)
    unlink(filler_path)
  }
  on.exit({
    sink(previous, type="message")
    close(connection)
    unlink(path)
  }, add=T)
  sink(connection, type="message")
  connections_before = showConnections(all=T)
  run()
  stopifnot(identical(sink.number(type="message"), as.integer(connection)),
    identical(showConnections(all=T), connections_before))
  message("caller-message-after-fit")
  flush(connection)
  stopifnot(identical(readLines(path), "caller-message-after-fit"))
}

normal = data.frame(dose=seq(0.1, 1, length.out=11L), response=seq(1, 3, length.out=11L)^2)
failing = data.frame(dose=exp(seq(log(1e-9), log(1e9), length.out=6L)),
  response=c(-1e300, 1e300, -1e300, 1e300, -1e300, 1e300))
# Keep any optimizer-created default graphics files outside the repository.
working_directory = getwd()
workspace = tempfile("runtime-state-")
dir.create(workspace)
setwd(workspace)
for (data in list(normal, failing)) {
  for (reuse_lower_connection in c(F, T)) {
    check_message_capture(function() {
      suppressWarnings(runmodelsaic_drc(data, "dose", "response"))
    }, reuse_lower_connection)
    check_message_capture(function() {
      aggregate = aggregate_grna_expression(
        data.frame(grna=seq_len(NROW(data)), qpcr=data$dose, bulk=data$response),
        "grna", "qpcr", "bulk", expected_replicates=NULL, expected_points=NULL)
      suppressWarnings(fit_qpcr_bulk_rna_models(aggregate))
    }, reuse_lower_connection)
    check_message_capture(function() {
      suppressWarnings(study_environment$fit_ll4(data.frame(x=data$dose, y=data$response),
        min(data$response), max(data$response)))
    }, reuse_lower_connection)
  }
}
setwd(working_directory)
unlink(workspace, recursive=T)
message("LL.4 message capture preserves the caller's connection and routing.")
