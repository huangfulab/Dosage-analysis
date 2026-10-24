#!/usr/bin/env Rscript

source("ci/dependencies.R")
expect_error = function(expression, text) {
  error = tryCatch({ force(expression); NULL }, error=identity)
  stopifnot(inherits(error, "error"), grepl(text, conditionMessage(error), fixed=T))
}
for (minor in 0:6) {
  plan = minimum_plan(paste0("4.", minor))
  stopifnot(NROW(plan$packages) == 9L, !anyDuplicated(plan$packages$package))
  installed = cbind(Package=plan$packages$package, Version=plan$packages$version)
  verify_dependency_versions(plan, installed=installed, load=F)
  installed[1L, "Version"] = "99.0.0"
  expect_error(verify_dependency_versions(plan, installed=installed, load=F), "not silently upgrade")
  expect_error(verify_dependency_versions(plan, installed=installed[-1L, ], load=F), "Missing minimum dependencies")
}
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

# Verify the install guard without accessing the network or modifying a library.
check_install_guard = function() {
  previous = Sys.getenv("GITHUB_ACTIONS", unset=NA_character_)
  on.exit(if (is.na(previous)) Sys.unsetenv("GITHUB_ACTIONS") else Sys.setenv(GITHUB_ACTIONS=previous))
  Sys.unsetenv("GITHUB_ACTIONS")
  expect_error(install_minimum_dependencies(minimum_plan("4.0")), "restricted to disposable GitHub Actions")
}
check_install_guard()
message("CI baseline, exact-version, dependency-constraint, and no-local-install tests passed.")
devtools::session_info(pkgs="attached")
