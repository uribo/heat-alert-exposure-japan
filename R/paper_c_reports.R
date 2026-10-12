# NOTE: readable analysis code from the private project. It needs the fire-department
# transport records, which are not distributed, and is not sourced by the public
# pipeline (_targets.R). Kept so that the published estimates can be audited.

# Writers for the reports the manuscript reads (Issue #50 (a)). The reports
# are targets; these write them where draft/paper_c.qmd reads them, as the
# commands of `format = "file"` targets.

# Write a report object through save_rds_atomic(), so that an interrupted
# target leaves the previous report in place rather than a truncated one, and
# return the path for the file target.
write_report_rds <- function(report, path) {
  save_rds_atomic(report, path)
  path
}

# Write each table of a named list to the CSV of the same name in `paths`
# and return the paths for the file target. Each CSV goes through a temporary
# file in the same directory and is renamed into place, as in
# save_rds_atomic().
write_report_csvs <- function(tables, paths) {
  if (!setequal(names(tables), names(paths))) {
    rlang::abort("Every table needs exactly one path.")
  }
  purrr::iwalk(tables, \(table, name) {
    path <- paths[[name]]
    temporary <- tempfile(
      pattern = paste0(basename(path), "."),
      tmpdir = dirname(path),
      fileext = ".tmp"
    )
    readr::write_csv(table, temporary)
    if (!file.rename(temporary, path)) {
      unlink(temporary)
      rlang::abort(glue::glue("Could not move {temporary} to {path}."))
    }
  })
  unname(paths[names(tables)])
}

# Eight-category place-of-onset composition (paper C, step 5g (iv)).
# Descriptive only, no test and no alert split: 7.15 fixes the eight
# categories as "reported descriptively, without tests". Rows: the stratified
# outcome panel restricted to department-days with a recorded place breakdown,
# where the eight place columns sum to transports (asserted). Counts and
# shares only; no department name leaves this function.
build_paper_c_place_composition <- function(panel) {
  place_cols <- c(
    place_home = "Home",
    place_work_industry = "Workplace (factories and similar)",
    place_work_agri = "Workplace (agriculture and livestock)",
    place_school = "Educational institution",
    place_public_indoor = "Public place, indoor",
    place_public_outdoor = "Public place, outdoor",
    place_road = "Road",
    place_other = "Other"
  )
  stopifnot(all(names(place_cols) %in% names(panel)))
  # WBGT-complete rows only, so that the composition covers the same
  # department-days as the stratified outcome fits (tbl-outcomes).
  recorded <- panel |>
    dplyr::filter(!is.na(.data$local_high), !is.na(.data$place_home))
  place_sum <- rowSums(recorded[, names(place_cols)])
  stopifnot(
    !anyNA(place_sum),
    all(place_sum == recorded$transports)
  )
  composition <- tibble::tibble(
    category = unname(place_cols),
    cases = as.integer(colSums(recorded[, names(place_cols)]))
  ) |>
    dplyr::mutate(share = cases / sum(cases))
  meta <- tibble::tibble(
    department_days = nrow(recorded),
    department_days_panel = nrow(panel),
    prefectures = dplyr::n_distinct(recorded$prefecture),
    departments = dplyr::n_distinct(recorded$department),
    date_min = as.character(min(recorded$date)),
    date_max = as.character(max(recorded$date)),
    transports = sum(recorded$transports)
  )
  stopifnot(sum(composition$cases) == meta$transports)
  list(composition = composition, meta = meta)
}

# Sample counts that the Sample section of the manuscript states (Issue #59,
# stage 6): prefectures, departments, prefecture-years and their completeness,
# blank records in the primary sample and in the stratified prefectures, and
# department-days absent within the department-years present, by prefecture
# and month. A prefecture-year is complete when every date of the 122-day
# window is present, whether or not every department reports on it. The
# stratified prefectures are those whose outcomes are defined in
# `panel_outcomes`. The distance from each department's jurisdiction centroid
# to its assigned WBGT station (Issue #43, item 12) is read from `covariates`,
# the department table of result 1 whose `dist_m` the result 1 regressions
# log, so the Sample text and the covariate cannot drift apart. Counts only;
# no department name leaves this function.
build_paper_c_sample_summary <- function(
  panel,
  panel_outcomes,
  covariates,
  window_days = 122L
) {
  stopifnot(
    nrow(panel_outcomes) == nrow(panel),
    "y_mild" %in% names(panel_outcomes)
  )
  panel_departments <- dplyr::distinct(panel, prefecture, department)
  key <- dplyr::join_by("prefecture", "department")
  if (
    anyDuplicated(covariates[c("prefecture", "department")]) > 0L ||
      nrow(dplyr::anti_join(panel_departments, covariates, by = key)) > 0L ||
      nrow(dplyr::anti_join(covariates, panel_departments, by = key)) > 0L
  ) {
    rlang::abort(glue::glue(
      "The covariate table does not hold one row for each of the ",
      "{nrow(panel_departments)} panel departments."
    ))
  }
  # `[[` with exact matching: a missing column is NULL rather than a partial
  # match, and NULL fails the length check instead of passing as numeric(0).
  dist_m <- covariates[["dist_m", exact = TRUE]]
  if (
    !is.numeric(dist_m) ||
      length(dist_m) != nrow(panel_departments) ||
      !all(is.finite(dist_m)) ||
      any(dist_m <= 0)
  ) {
    rlang::abort(glue::glue(
      "`dist_m` must hold one finite positive distance for each of the ",
      "{nrow(panel_departments)} panel departments."
    ))
  }
  dist_km <- dist_m / 1000
  dates_by_year <- panel |>
    dplyr::summarise(
      n_dates = dplyr::n_distinct(.data$date),
      .by = c("prefecture", "year")
    )
  outcome_prefectures <- unique(
    panel_outcomes$prefecture[!is.na(panel_outcomes$y_mild)]
  )
  blank <- is.na(panel$transports)
  missing <- panel |>
    dplyr::summarise(
      n_rows = dplyr::n(),
      .by = c("prefecture", "department", "year")
    ) |>
    dplyr::filter(.data$n_rows < window_days) |>
    dplyr::select("prefecture", "department", "year")
  absent <- missing |>
    dplyr::inner_join(
      panel |>
        dplyr::distinct(prefecture, year, date),
      by = dplyr::join_by("prefecture", "year"),
      relationship = "many-to-many"
    ) |>
    dplyr::anti_join(
      panel,
      by = dplyr::join_by("prefecture", "department", "date")
    ) |>
    dplyr::mutate(month = as.integer(format(.data$date, "%m"))) |>
    dplyr::count(prefecture, year, month, name = "department_days")
  summary <- tibble::tibble(
    department_days = nrow(panel),
    prefectures = dplyr::n_distinct(panel$prefecture),
    departments = dplyr::n_distinct(panel$department),
    prefecture_years = nrow(dates_by_year),
    complete_prefecture_years = sum(dates_by_year$n_dates == window_days),
    blank_records = sum(blank),
    stratified_prefectures = length(outcome_prefectures),
    blank_records_stratified = sum(
      blank & panel$prefecture %in% outcome_prefectures
    ),
    absent_department_days = sum(absent$department_days),
    station_dist_km_median = stats::median(dist_km),
    station_dist_km_max = max(dist_km)
  )
  list(summary = summary, absent = absent)
}

# A report without the fields that change on every rebuild without changing
# any estimate: wall times, the worker count and the session record (Issue
# #50, 2026-09-30 decision 1). Moved reports are checked against their
# previous copies with identical() on what this returns. List elements are
# dropped by name at any depth; `columns` are dropped from data frames.
report_without_volatile <- function(
  report,
  fields = c(
    "timings",
    "session",
    "seconds_fit",
    "workers",
    "verified_at",
    "seconds",
    "total_wall_seconds",
    "seconds_null_fit",
    "seconds_full_fit"
  ),
  columns = "seconds"
) {
  if (is.data.frame(report)) {
    return(report[, setdiff(names(report), columns), drop = FALSE])
  }
  if (!is.list(report)) {
    return(report)
  }
  kept <- if (is.null(names(report))) {
    report
  } else {
    report[!names(report) %in% fields]
  }
  kept[] <- lapply(kept, report_without_volatile, fields, columns)
  kept
}
