#!/usr/bin/env Rscript

# CI utility, not an analysis runner. Installation requires BOTH --install and
# GitHub Actions so ordinary local tests never alter the user's R libraries.
minimum_plan = function(r_series, root="ci") {
  versions = utils::read.csv(file.path(root, "r-versions.csv"), stringsAsFactors=F, colClasses="character")
  packages = utils::read.csv(file.path(root, "minimum-packages.csv"), stringsAsFactors=F,
    check.names=F, colClasses="character")
  row = versions[versions$r == r_series, , drop=F]
  stopifnot(NROW(row) == 1L, r_series %in% names(packages))
  plan = data.frame(package=packages$package, version=packages[[r_series]],
    source=packages$source, stringsAsFactors=F)
  # The shared 2022 snapshot's locfit requires R >= 4.1. DESeq2 on R 4.0
  # therefore needs the last R-4.0-compatible source, not a newer direct floor.
  compatibility = data.frame(package="locfit", version="1.5-9.4", source="CRAN", stringsAsFactors=F)
  if (r_series != "4.0") compatibility = compatibility[FALSE, , drop=F]
  # The R 4.6 snapshot's cowplot 1.2.0 (via drc/car/doBy) requires
  # ggplot2 >= 3.5.2. Keep our direct floor 3.4.3 and pin a compatible cowplot.
  if (r_series == "4.6") {
    compatibility = data.frame(package="cowplot", version="1.1.3", source="CRAN", stringsAsFactors=F)
  }
  list(r=r_series, bioconductor=row$bioconductor, cran_snapshot=row$cran_snapshot,
    packages=plan, compatibility_packages=compatibility)
}

verify_dependency_versions = function(plan, installed=utils::installed.packages(), load=T) {
  wanted = rbind(plan$compatibility_packages, plan$packages)
  selected = match(wanted$package, installed[, "Package"])
  if (anyNA(selected)) stop("Missing minimum dependencies: ", paste(wanted$package[is.na(selected)], collapse=", "))
  actual = installed[selected, "Version"]
  equal = vapply(seq_along(actual), function(i) {
    package_version(actual[[i]]) == package_version(wanted$version[[i]])
  }, logical(1))
  if (!all(equal)) {
    stop("CI must use exact assumed minimum versions, not silently upgrade: ",
      paste(sprintf("%s=%s (expected %s)", wanted$package[!equal], actual[!equal], wanted$version[!equal]), collapse=", "))
  }
  if (load) {
    loaded = vapply(wanted$package, requireNamespace, logical(1), quietly=T)
    if (!all(loaded)) stop("Minimum packages failed to load: ", paste(wanted$package[!loaded], collapse=", "))
  }
  invisible(wanted)
}

# Validate hard dependency constraints after all exact pins are installed.
# Pinning a direct package must never leave its dependency graph inconsistent.
verify_hard_dependencies = function(installed=utils::installed.packages()) {
  fields = intersect(c("Depends", "Imports", "LinkingTo"), colnames(installed))
  for (i in seq_len(NROW(installed))) {
    dependencies = installed[i, fields]
    dependencies = dependencies[!is.na(dependencies) & nzchar(dependencies)]
    for (dependency in trimws(unlist(strsplit(dependencies, ",", fixed=T)))) {
      parsed = regmatches(dependency, regexec(
        "^([A-Za-z][A-Za-z0-9.]*)\\s*(?:\\(\\s*(>=|<=|==|>|<|!=)\\s*([0-9][0-9.-]*)\\s*\\))?$",
        dependency, perl=T
      ))[[1L]]
      if (length(parsed) == 0L) stop("Unrecognized hard dependency: ", dependency)
      name = parsed[[2L]]
      if (name == "R") {
        actual = getRversion()
      } else {
        index = match(name, installed[, "Package"])
        if (is.na(index)) stop(installed[i, "Package"], " has a missing hard dependency: ", name)
        actual = package_version(installed[index, "Version"])
      }
      if (length(parsed) >= 4L && nzchar(parsed[[3L]])) {
        expected = package_version(parsed[[4L]])
        valid = switch(parsed[[3L]], ">="=actual >= expected, "<="=actual <= expected,
          "=="=actual == expected, ">"=actual > expected, "<"=actual < expected, "!="=actual != expected)
        if (!valid) stop(installed[i, "Package"], " requires ", dependency, "; installed ", actual)
      }
    }
  }
  invisible(T)
}

write_ci_environment = function(log_dir) {
  dir.create(log_dir, recursive=T, showWarnings=F)
  inventory = as.data.frame(utils::installed.packages(), stringsAsFactors=F)
  utils::write.table(inventory, file.path(log_dir, "installed-packages.tsv"), sep="\t", row.names=F, quote=T)
  writeLines(capture.output(sessionInfo()), file.path(log_dir, "session-info.txt"))
}

exact_source_urls = function(package, version, source, bioc) {
  if (source == "Bioconductor") {
    return(sprintf("%s/bioc/src/contrib/%s_%s.tar.gz", bioc, package, version))
  }
  c(sprintf("https://cran.r-project.org/src/contrib/Archive/%s/%s_%s.tar.gz", package, package, version),
    sprintf("https://cran.r-project.org/src/contrib/%s_%s.tar.gz", package, version))
}

# Pause only on GitHub Actions. Unit tests inject pause=function(attempt) NULL.
ci_retry_pause = function(attempt) {
  if (identical(Sys.getenv("GITHUB_ACTIONS"), "true")) Sys.sleep(min(2^(attempt - 1L), 30))
  invisible(attempt)
}

# libcurl inside download.file does not honour curl's --retry flags.
# The external curl binary does, including the SSL connect errors seen in CI.
configure_ci_downloads = function() {
  method = if (nzchar(Sys.which("curl"))) "curl" else "libcurl"
  extra = if (identical(method, "curl")) "--retry 5 --retry-connrefused --retry-delay 2 --retry-all-errors" else ""
  options(download.file.method=method, download.file.extra=extra)
  invisible(list(method=method, extra=extra))
}

# A still-current CRAN release has no Archive entry. Try both locations for
# the SAME exact version; never interpret a missing archive as permission to
# upgrade the floor. The injected downloader permits offline fallback tests.
download_exact_source = function(urls, destination, download=utils::download.file, attempts=5L, pause=ci_retry_pause) {
  for (attempt in seq_len(attempts)) {
    for (url in urls) {
      unlink(destination)
      downloaded = tryCatch({
        status = suppressWarnings(download(url, destination, mode="wb"))
        ok = isTRUE(status == 0L) && file.exists(destination) && isTRUE(file.info(destination)$size > 0)
        if (!isTRUE(ok)) message("Source unavailable: ", url, "; download status ", paste(status))
        ok
      }, error=function(error) {
        message("Source unavailable: ", url, "; ", conditionMessage(error))
        FALSE
      })
      if (isTRUE(downloaded)) return(invisible(url))
    }
    if (attempt < attempts) pause(attempt)
  }
  unlink(destination)
  stop("Could not download exact package source: ", paste(urls, collapse=", "))
}

exact_version_installed = function(package, version, library) {
  installed = utils::installed.packages(lib.loc=library, noCache=TRUE)
  index = match(package, installed[, "Package"])
  if (is.na(index)) return(FALSE)
  package_version(installed[index, "Version"]) == package_version(version)
}

# remotes::install_local can return after a dependency download failure
# without leaving the requested package installed. Retry that package only.
install_exact_source = function(package, version, urls, source_archive, library, repositories,
    attempts=3L, pause=ci_retry_pause) {
  for (attempt in seq_len(attempts)) {
    source_url = download_exact_source(urls, source_archive)
    message("Installing assumed minimum: ", package, " ", version, " from ", source_url,
      " (attempt ", attempt, " of ", attempts, ")")
    install_error = tryCatch({
      remotes::install_local(source_archive, lib=library, dependencies=NA, upgrade="never",
        build_vignettes=F, repos=repositories, quiet=F)
      NULL
    }, error=function(error) error)
    unlink(source_archive)
    if (isTRUE(exact_version_installed(package, version, library))) return(invisible(source_url))
    detail = if (is.null(install_error)) "installer returned without the required version" else conditionMessage(install_error)
    message("Attempt ", attempt, " did not install ", package, " ", version, ": ", detail)
    if (attempt < attempts) pause(attempt)
  }
  stop("Failed to install ", package, " ", version, " after ", attempts, " attempts")
}

install_minimum_dependencies = function(plan) {
  if (!identical(Sys.getenv("GITHUB_ACTIONS"), "true")) {
    stop("--install is restricted to disposable GitHub Actions libraries")
  }
  library = Sys.getenv("CI_R_LIBRARY")
  if (!nzchar(library) || !identical(normalizePath(.libPaths()[[1L]], mustWork=T), normalizePath(library, mustWork=T))) {
    stop("CI_R_LIBRARY must be the first R library; refusing a user/global library install")
  }
  cran = paste0("https://packagemanager.posit.co/cran/", plan$cran_snapshot)
  bioc = paste0("https://bioconductor.org/packages/", plan$bioconductor)
  repositories = c(BioCsoft=paste0(bioc, "/bioc"), BioCann=paste0(bioc, "/data/annotation"),
    BioCexp=paste0(bioc, "/data/experiment"), CRAN=cran)
  options(repos=repositories, timeout=600, Ncpus=2L, pkgType="source")
  configure_ci_downloads()
  log_dir = Sys.getenv("CI_LOG_DIR", unset=tempdir())
  dir.create(log_dir, recursive=T, showWarnings=F)
  on.exit(write_ci_environment(log_dir), add=T)
  packages = rbind(plan$compatibility_packages, plan$packages)
  utils::write.csv(packages, file.path(log_dir, "assumed-minimum-packages.csv"), row.names=F)
  writeLines(c(paste("R series:", plan$r), paste("Bioconductor:", plan$bioconductor),
    paste("CRAN snapshot:", cran)), file.path(log_dir, "repositories.txt"))

  # remotes is installer tooling only; never install Suggests or upgrade all
  # packages. Archived direct dependencies are installed from exact sources.
  bootstrap = tempfile(fileext=".tar.gz")
  on.exit(unlink(bootstrap), add=T)
  download_exact_source("https://cran.r-project.org/src/contrib/Archive/remotes/remotes_2.4.2.tar.gz", bootstrap)
  utils::install.packages(bootstrap, repos=NULL, type="source", lib=library)
  if (!requireNamespace("remotes", quietly=T)) stop("Could not bootstrap remotes")
  # Compatibility pins must precede direct packages so dependency resolution
  # sees a supported installed version before attempting the snapshot version.
  source_archive = tempfile(fileext=".tar.gz")
  on.exit(unlink(source_archive), add=T)
  for (i in seq_len(NROW(packages))) {
    package = packages[i, ]
    urls = exact_source_urls(package$package, package$version, package$source, bioc)
    install_exact_source(package$package, package$version, urls, source_archive, library, repositories)
  }
  # Dependencies installed for a later direct package may have replaced an
  # earlier pin. Any such replacement is a failure, not a minimum-test pass.
  verify_dependency_versions(plan)
  verify_hard_dependencies()
  message("PASS: exact direct dependency floors and their hard dependency constraints verified")
}

if (sys.nframe() == 0L) {
  args = commandArgs(trailingOnly=T)
  mode = if (length(args)) args[[1L]] else "--plan"
  if (!mode %in% c("--plan", "--install", "--verify", "--report")) stop("Use --plan, --install, --verify, or --report")
  r_series = paste(R.version$major, strsplit(R.version$minor, ".", fixed=T)[[1L]][[1L]], sep=".")
  if (mode == "--report") {
    write_ci_environment(Sys.getenv("CI_LOG_DIR", unset=tempdir()))
  } else {
    plan = minimum_plan(r_series)
    if (mode == "--plan") print(plan)
    if (mode == "--install") install_minimum_dependencies(plan)
    if (mode == "--verify") {
      verify_dependency_versions(plan)
      verify_hard_dependencies()
    }
  }
}
