# NOTE: readable analysis code from the private project. It needs the fire-department
# transport records, which are not distributed, and is not sourced by the public
# pipeline (_targets.R). Kept so that the published estimates can be audited.

# Shared helpers for the paper-C command-line drivers.

parse_args <- function(
  defaults,
  args = commandArgs(trailingOnly = TRUE)
) {
  for (arg in args) {
    key <- sub("^--", "", arg)
    if (!grepl("=", key, fixed = TRUE)) {
      # A bare switch must carry the `--` prefix, as each driver's own parser
      # required before these were merged.
      if (
        startsWith(arg, "--") &&
          key %in% names(defaults) &&
          is.logical(defaults[[key]])
      ) {
        defaults[[key]] <- TRUE
        next
      }
      rlang::abort(glue::glue("Unknown argument: {arg}"))
    }
    parts <- strsplit(key, "=", fixed = TRUE)[[1L]]
    if (length(parts) != 2L || !parts[[1L]] %in% names(defaults)) {
      rlang::abort(glue::glue("Unknown argument: {arg}"))
    }
    value <- parts[[2L]]
    defaults[[parts[[1L]]]] <- if (is.integer(defaults[[parts[[1L]]]])) {
      as.integer(value)
    } else if (is.logical(defaults[[parts[[1L]]]])) {
      parsed <- switch(
        tolower(value),
        true = TRUE,
        false = FALSE,
        rlang::abort(glue::glue(
          "Logical argument {parts[[1L]]} must be true or false."
        ))
      )
      parsed
    } else {
      value
    }
  }
  defaults
}

resolve_path <- function(path) {
  if (startsWith(path, "/")) path else here::here(path)
}

# Path of the resumable partial file that sits next to a report. The report
# path must end in ".rds": otherwise sub() would return it unchanged, the
# partial results would be written over the report path, and the final
# unlink() of the partial file would delete the finished report.
partial_path_for <- function(output_path, suffix = ".partial.rds") {
  if (!grepl("\\.rds$", output_path)) {
    rlang::abort(glue::glue(
      "--output must end in .rds so that the partial file gets its own path; ",
      "got {output_path}."
    ))
  }
  partial_path <- sub("\\.rds$", suffix, output_path)
  stopifnot(!identical(partial_path, output_path))
  partial_path
}

clamp_workers <- function(workers) {
  available_cores <- parallel::detectCores()
  if (!is.na(available_cores)) min(workers, available_cores) else workers
}

saitama_departments <- function(panel) {
  panel |>
    dplyr::filter(prefecture == "埼玉県") |>
    dplyr::distinct(department) |>
    dplyr::pull(department)
}

# Data-use guard ---------------------------------------------------------------

# Saitama's data-use conditions keep its department names out of any output
# (CLAUDE.md), so a report is checked before it is written. The check walks
# the serialized bytes rather than the object tree: a report carries
# `utils::sessionInfo()`, whose loaded-namespace entries are self-referential
# once enough namespaces are attached, and a recursive walk over that never
# terminates. Serializing is bounded, and it also sees names hidden in
# attributes, factor levels and environments, which a walk over list elements
# would miss. The one thing it cannot see is a string that does not exist
# yet: an ALTREP deferred character vector is serialized in its compact
# form, so a name that would only be materialized by as.character() would be
# missed. Department names come from the data as character, so none is.
report_contains_string <- function(object, targets) {
  connection <- rawConnection(raw(0), "w")
  on.exit(close(connection))
  serialize(object, connection)
  blob <- rawConnectionValue(connection)
  matched <- purrr::map_lgl(targets, function(target) {
    length(grepRaw(charToRaw(target), blob, fixed = TRUE, all = FALSE)) > 0L
  })
  targets[matched]
}

check_no_saitama <- function(object, departments, label) {
  found <- report_contains_string(object, departments)
  if (length(found) > 0L) {
    rlang::abort(glue::glue(
      "{label} carries {length(found)} Saitama department name(s)."
    ))
  }
  invisible(object)
}

strip_fit <- function(spec) {
  spec$fit <- NULL
  spec$design_kept_rows <- NULL
  spec
}

paper_c_session <- function(packages = c("fixest", "sf", "targets")) {
  versions <- stats::setNames(
    lapply(packages, utils::packageVersion),
    packages
  )
  c(
    list(r = R.version.string),
    versions,
    # As text, not as the object: sessionInfo()'s loaded-namespace entries
    # are self-referential once enough namespaces are loaded.
    list(session_info = utils::capture.output(utils::sessionInfo()))
  )
}
