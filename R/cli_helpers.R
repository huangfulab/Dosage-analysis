# Shared dependency-free command-line and reproducibility helpers.

# Parse only alternating --option value pairs so CLI failures are deterministic
# and unsupported flag/repeated-option forms cannot be silently misinterpreted.
parse_cli_options <- function(arguments, usage, allowed_options) {
  if (length(arguments) == 0L || "--help" %in% arguments) {
    usage()
    quit(save = "no", status = 0L)
  }
  if (length(arguments) %% 2L != 0L) {
    stop("Every command-line option must have a value. Use `--help` for usage.", call. = FALSE)
  }
  keys <- arguments[seq.int(1L, length(arguments), by = 2L)]
  values <- arguments[seq.int(2L, length(arguments), by = 2L)]
  if (any(!startsWith(keys, "--"))) {
    stop("Options must start with `--`.", call. = FALSE)
  }
  keys <- substring(keys, 3L)
  if (anyDuplicated(keys)) {
    stop("Each command-line option may be supplied only once.", call. = FALSE)
  }
  unknown <- setdiff(keys, allowed_options)
  if (length(unknown) > 0L) {
    stop(
      sprintf("Unknown command-line option(s): --%s", paste(unknown, collapse = ", --")),
      call. = FALSE
    )
  }
  result <- as.list(values)
  names(result) <- keys
  result
}

option_or <- function(options, name, default = NULL) {
  if (name %in% names(options)) options[[name]] else default
}

parse_cli_boolean <- function(value, name) {
  normalized <- tolower(value)
  if (!normalized %in% c("true", "false")) {
    stop(sprintf("`--%s` must be `true` or `false`.", name), call. = FALSE)
  }
  identical(normalized, "true")
}

parse_cli_positive_integer <- function(value, name) {
  parsed <- suppressWarnings(as.numeric(value))
  if (length(parsed) != 1L || !is.finite(parsed) || parsed < 1L ||
      parsed != floor(parsed) || parsed > .Machine$integer.max) {
    stop(sprintf("`--%s` must be one positive integer.", name), call. = FALSE)
  }
  as.integer(parsed)
}

parse_cli_numeric_list <- function(value, name, integer = FALSE) {
  parsed <- suppressWarnings(as.numeric(strsplit(value, ",", fixed = TRUE)[[1L]]))
  if (length(parsed) == 0L || any(!is.finite(parsed))) {
    stop(sprintf("`--%s` must be a comma-separated numeric list.", name), call. = FALSE)
  }
  if (integer && any(parsed != as.integer(parsed))) {
    stop(sprintf("`--%s` must contain integers.", name), call. = FALSE)
  }
  if (integer) as.integer(parsed) else parsed
}

# Safety boundary: overwrite = TRUE recursively removes every immediate child
# of `path`; callers must select and validate the output directory beforehand.
ensure_output_directory <- function(path, overwrite = FALSE, protected_paths=character()) {
  if (file.exists(path) && !dir.exists(path)) {
    stop(sprintf("Output path exists but is not a directory: %s", path), call.=F)
  }
  # Resolve aliases with initially missing components before checking contents
  # or deleting anything. Never remove an input or executing source file.
  dir.create(path, recursive=T, showWarnings=F)
  if (!dir.exists(path)) stop("Could not create output directory: ", path, call.=F)
  canonical = normalizePath(path, winslash="/", mustWork=T)
  protected = c(
    vapply(protected_paths, normalizePath, character(1), winslash="/", mustWork=T),
    file.path(vapply(dirname(protected_paths), normalizePath, character(1), winslash="/", mustWork=T),
      basename(protected_paths))
  )
  if (any(protected == canonical | startsWith(protected, paste0(sub("/$", "", canonical), "/"))) ||
      grepl("^([A-Za-z]:)?/$", canonical) || grepl("^//[^/]+/[^/]+/?$", canonical)) {
    stop("Output directory must not contain input or source files.", call.=F)
  }
  entries <- list.files(path, all.files=TRUE, no..=TRUE, include.dirs=TRUE)
  if (length(entries) > 0L && !isTRUE(overwrite)) {
    stop("Output directory is not empty; use `--overwrite true` to replace its contents.", call.=F)
  }
  if (isTRUE(overwrite) && length(entries) > 0L) {
    unlink(file.path(path, entries), recursive=TRUE, force=TRUE)
  }
  invisible(path)
}

# MD5 values are provenance/integrity identifiers, not security guarantees.
file_md5 <- function(path) {
  if (!file.exists(path)) stop(sprintf("Cannot hash missing file: %s", path), call. = FALSE)
  unname(as.character(tools::md5sum(path)))
}

# Capture the executing R and package environment needed to reproduce a run.
write_base_session_info <- function(path) {
  report <- capture.output(sessionInfo())
  writeLines(report, path)
  invisible(report)
}

.json_escape <- function(value) {
  value <- enc2utf8(as.character(value))
  value = gsub("\\", "\\\\", value, fixed=T)
  value = gsub("\"", "\\\"", value, fixed=T)
  # JSON forbids literal control characters; R strings cannot contain NUL.
  for (code in 1:31) {
    value = gsub(intToUtf8(code), sprintf("\\u%04x", code), value, fixed=T)
  }
  value
}

# Encode the small, known provenance data model without adding a JSON package
# dependency to a command-line workflow intended to run with base R alone.
.json_encode <- function(value) {
  if (is.null(value) || length(value) == 0L) return("null")
  if (is.data.frame(value)) {
    rows <- lapply(seq_len(NROW(value)), function(index) as.list(value[index, , drop = FALSE]))
    return(.json_encode(rows))
  }
  if (is.list(value)) {
    names_value <- names(value)
    if (!is.null(names_value) && all(nzchar(names_value))) {
      fields <- vapply(seq_along(value), function(index) {
        paste0("\"", .json_escape(names_value[[index]]), "\":", .json_encode(value[[index]]))
      }, character(1L))
      return(paste0("{", paste(fields, collapse = ","), "}"))
    }
    return(paste0("[", paste(vapply(value, .json_encode, character(1L)), collapse = ","), "]"))
  }
  if (is.logical(value)) {
    encoded <- ifelse(is.na(value), "null", ifelse(value, "true", "false"))
  } else if (is.numeric(value) || is.integer(value)) {
    encoded <- ifelse(is.finite(value), format(value, scientific=FALSE, trim=TRUE, digits=17L, decimal.mark="."), "null")
  } else {
    encoded <- ifelse(is.na(value), "null", paste0("\"", .json_escape(value), "\""))
  }
  if (length(encoded) == 1L) encoded[[1L]] else paste0("[", paste(encoded, collapse = ","), "]")
}

write_analysis_provenance <- function(path, workflow, inputs, settings, sources, summary, output_directory) {
  # Inventory is captured last so it records hashes/sizes of all preceding
  # tabular, graphical, RDS, manifest, and session outputs in one portable file.
  output_files <- sort(list.files(output_directory, recursive = TRUE, full.names = TRUE))
  output_files <- output_files[file.info(output_files)$isdir %in% FALSE]
  artifacts <- lapply(output_files, function(file) {
    info <- file.info(file)
    list(
      path = file,
      size_bytes = unname(info$size),
      md5 = file_md5(file)
    )
  })
  provenance <- list(
    schema_version = "1.0",
    workflow = workflow,
    run_time_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    inputs = inputs,
    settings = settings,
    sources = sources,
    summary = summary,
    artifacts = artifacts
  )
  writeLines(.json_encode(provenance), path, useBytes = TRUE)
  invisible(provenance)
}
