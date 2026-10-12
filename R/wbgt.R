# Guard for the MOE WBGT "final" hourly files (file_type "2-D").
#
# parse_moe_wbgt_csv() renames four columns with purrr::set_names(), so a file
# that reads with a different width aborts with a 4-vs-2 message naming no path.
# On 2026-09-02 that killed a 12-hour build at file 6,889 of 35,082, where
# OneDrive Files On-Demand had served a partially hydrated file. A truncated file
# of the right width would have passed silently instead, so check the row count
# too: every published file covers one whole month at hourly resolution, which
# makes the expected shape a function of the year-month in the file name.
#
# R/read_moe_wbgt.R is vendored and must stay byte-identical to upstream, so the
# guard lives here and wraps the vendored parser rather than editing it.
parse_moe_wbgt_final <- function(path) {
  header <- readr::read_lines(path, n_max = 1L)
  if (length(header) != 1L) {
    rlang::abort(glue::glue("{path}: no header line."))
  }
  n_field <- stringr::str_count(header, ",") + 1L
  if (n_field != 4L) {
    rlang::abort(glue::glue(
      "{path}: expected a 4-column header (Date, Time, WBGT, Tg), got {n_field}."
    ))
  }
  year_month <- stringr::str_extract(basename(path), "[0-9]{6}(?=\\.csv$)")
  if (is.na(year_month)) {
    rlang::abort(glue::glue("{path}: the file name carries no YYYYMM."))
  }
  n_expected <- as.integer(lubridate::days_in_month(lubridate::ym(
    year_month
  ))) *
    24L
  df <- parse_moe_wbgt_csv(path, file_type = "2-D")
  if (nrow(df) != n_expected) {
    rlang::abort(glue::glue(
      "{path}: expected {n_expected} hourly rows for {year_month}, got {nrow(df)}."
    ))
  }
  df
}
