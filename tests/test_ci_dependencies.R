#!/usr/bin/env Rscript

source("ci/dependencies.R")
expect_error = function(expression, text) {
  error = tryCatch({ force(expression); NULL }, error=identity)
  stopifnot(inherits(error, "error"), grepl(text, conditionMessage(error), fixed=T))
}
for (minor in 0:6) {
  plan = minimum_plan(paste0("4.", minor))
  stopifnot(NROW(plan$packages) == 9L, !anyDuplicated(plan$packages$package))
  wanted = rbind(plan$compatibility_packages, plan$packages)
  if (minor == 0L) {
    stopifnot(identical(plan$compatibility_packages$package, "locfit"),
      identical(plan$compatibility_packages$version, "1.5-9.4"))
  } else if (minor == 6L) {
    stopifnot(identical(plan$compatibility_packages$package, "cowplot"),
      identical(plan$compatibility_packages$version, "1.1.3"))
  } else stopifnot(NROW(plan$compatibility_packages) == 0L)
  installed = cbind(Package=wanted$package, Version=wanted$version)
  verify_dependency_versions(plan, installed=installed, load=F)
  installed[1L, "Version"] = "99.0.0"
  expect_error(verify_dependency_versions(plan, installed=installed, load=F), "not silently upgrade")
  expect_error(verify_dependency_versions(plan, installed=installed[-1L, ], load=F), "Missing minimum dependencies")
}
for (series in c("4.4", "4.5", "4.6")) {
  plan = minimum_plan(series)
  stopifnot(plan$packages$version[plan$packages$package == "ggplot2"] == "3.4.3")
  if (series != "4.4") stopifnot(plan$packages$version[plan$packages$package == "matrixStats"] == "1.5.0")
}
# A current source release may have no Archive entry; fallback must preserve
# the exact version and must discard failed/partial downloads.
urls = exact_source_urls("matrixStats", "1.5.0", "CRAN", "unused")
stopifnot(length(urls) == 2L, all(grepl("matrixStats_1.5.0.tar.gz", urls, fixed=T)))
archive = tempfile()
selected = download_exact_source(urls, archive, download=function(url, destination, mode) {
  writeLines("partial", destination)
  if (identical(url, urls[[1L]])) stop("404")
  writeLines("exact source", destination)
  0L
})
stopifnot(identical(selected, urls[[2L]]), identical(readLines(archive), "exact source"))
no_pause = function(attempt) NULL
expect_error(download_exact_source(urls, archive,
  download=function(...) 1L, pause=no_pause), "Could not download exact package source")
stopifnot(!file.exists(archive))
# A connect error on every URL is retried; the next attempt can succeed.
calls = 0L
pauses = 0L
selected = download_exact_source(urls, archive, download=function(url, destination, mode) {
  calls <<- calls + 1L
  if (calls <= length(urls)) stop("SSL connect error")
  writeLines("exact source", destination)
  0L
}, pause=function(attempt) { pauses <<- pauses + 1L })
stopifnot(identical(selected, urls[[1L]]), calls == length(urls) + 1L, pauses == 1L,
  identical(readLines(archive), "exact source"))
previous_method = getOption("download.file.method")
previous_extra = getOption("download.file.extra")
options(download.file.method="curl", download.file.extra="--retry 5")
stopifnot(identical(use_ci_libcurl_downloads(), "libcurl"),
  identical(getOption("download.file.method"), "libcurl"),
  is.null(getOption("download.file.extra")))
options(download.file.method=previous_method, download.file.extra=previous_extra)
expect_error(minimum_plan("4.7"), "NROW(row) == 1L")
installed = cbind(
  Package=c("one", "two"), Version=c("1.0", "2.0"),
  Depends=c("R (>= 4.0.0)", NA_character_),
  Imports=c("two (>= 2.0)", NA_character_), LinkingTo=NA_character_
)
verify_hard_dependencies(installed)
installed[1L, "Imports"] = "two (>= 3.0)"
expect_error(verify_hard_dependencies(installed), "requires two (>= 3.0)")
installed[1L, "Imports"] = "missingPackage"
expect_error(verify_hard_dependencies(installed), "missing hard dependency")
installed[1L, "Imports"] = "two (nonsense)"
expect_error(verify_hard_dependencies(installed), "Unrecognized hard dependency")

# The R 4.6 snapshot's cowplot 1.2.0 conflicts with the ggplot2 floor.
# A transitive compatibility pin must resolve that conflict without raising it.
installed = cbind(
  Package=c("cowplot", "ggplot2"), Version=c("1.2.0", "3.4.3"),
  Depends=c("R (>= 3.5.0)", NA_character_),
  Imports=c("ggplot2 (>= 3.5.2)", NA_character_), LinkingTo=NA_character_
)
expect_error(verify_hard_dependencies(installed), "cowplot requires ggplot2 (>= 3.5.2)")
installed[1L, "Version"] = minimum_plan("4.6")$compatibility_packages$version
installed[1L, "Imports"] = "ggplot2 (>= 3.4.0)"
verify_hard_dependencies(installed)

# Verify the install guard without accessing the network or modifying a library.
check_install_guard = function() {
  previous = Sys.getenv("GITHUB_ACTIONS", unset=NA_character_)
  on.exit(if (is.na(previous)) Sys.unsetenv("GITHUB_ACTIONS") else Sys.setenv(GITHUB_ACTIONS=previous))
  Sys.unsetenv("GITHUB_ACTIONS")
  expect_error(install_minimum_dependencies(minimum_plan("4.0")), "restricted to disposable GitHub Actions")
}
check_install_guard()
message("CI baseline, exact-version, dependency-constraint, and no-local-install tests passed.")
sessionInfo()
