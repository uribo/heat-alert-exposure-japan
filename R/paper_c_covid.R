# NOTE: readable analysis code from the private project. It needs the fire-department
# transport records, which are not distributed, and is not sourced by the public
# pipeline (_targets.R). Kept so that the published estimates can be audited.

# Department-day indicators of COVID-19 emergency declarations (prefecture
# level) and priority measures (municipal level, aggregated to the department
# through the municipality crosswalk), paper C Issue #43 item 11.

# Strip the district prefix from a municipal name: 石狩郡当別町 -> 当別町.
# Non-greedy so that 赤穂郡上郡町 -> 上郡町.
strip_district <- function(x) {
  ifelse(
    grepl("郡", x) & grepl("[町村]$", x),
    sub("^[^郡]*郡", "", x),
    x
  )
}

# Source corrections listed in data-raw/covid19.R.
covid_municipality_fixes <- c(
  "茂木市" = "茂木町",
  "旧唐津市" = "唐津市"
)

# Expand (from, to) periods to one row per date.
expand_periods <- function(x) {
  x |>
    dplyr::mutate(date = purrr::map2(date_from, date_to, seq, by = "day")) |>
    tidyr::unnest(date) |>
    dplyr::select(-date_from, -date_to)
}

# Returns a list: `indicators`, a tibble of department, date, emergency,
# priority (0/1) over `dates` for the departments of `crosswalk`; and
# `unmatched`, the priority-measure municipalities of the sample prefectures
# that no department in the crosswalk serves (reported, not imputed).
build_covid_indicators <- function(emergency, priority, crosswalk, dates) {
  stopifnot(
    all(c("area", "date_from", "date_to") %in% names(emergency)),
    all(c("prefecture", "area", "date_from", "date_to") %in% names(priority)),
    all(c("prefecture", "department", "administration") %in% names(crosswalk))
  )
  departments <- crosswalk |>
    dplyr::distinct(prefecture, department)
  window <- range(dates)
  emergency_days <- emergency |>
    dplyr::filter(date_to >= window[1], date_from <= window[2]) |>
    dplyr::transmute(prefecture = area, date_from, date_to) |>
    expand_periods() |>
    dplyr::distinct(prefecture, date) |>
    dplyr::mutate(emergency = 1)
  priority_named <- priority |>
    dplyr::filter(
      date_to >= window[1],
      date_from <= window[2],
      prefecture %in% departments$prefecture
    ) |>
    dplyr::mutate(
      administration = strip_district(area),
      administration = dplyr::coalesce(
        unname(covid_municipality_fixes[administration]),
        administration
      )
    )
  matched <- priority_named |>
    dplyr::inner_join(
      crosswalk |> dplyr::distinct(prefecture, administration, department),
      by = dplyr::join_by(prefecture, administration),
      relationship = "many-to-many"
    )
  unmatched <- priority_named |>
    dplyr::distinct(prefecture, area, administration) |>
    dplyr::anti_join(
      crosswalk |> dplyr::distinct(prefecture, administration),
      by = dplyr::join_by(prefecture, administration)
    )
  priority_days <- matched |>
    dplyr::select(prefecture, department, date_from, date_to) |>
    expand_periods() |>
    dplyr::distinct(prefecture, department, date) |>
    dplyr::mutate(priority = 1)
  indicators <- tidyr::expand_grid(departments, date = dates) |>
    dplyr::left_join(emergency_days, by = dplyr::join_by(prefecture, date)) |>
    dplyr::left_join(
      priority_days,
      by = dplyr::join_by(prefecture, department, date)
    ) |>
    dplyr::mutate(
      emergency = tidyr::replace_na(emergency, 0),
      priority = tidyr::replace_na(priority, 0)
    )
  stopifnot(
    nrow(indicators) == nrow(departments) * length(dates),
    !anyNA(indicators$emergency),
    !anyNA(indicators$priority)
  )
  list(indicators = indicators, unmatched = unmatched)
}

# Provenance of the two committed source files (data-raw/covid19.R documents
# the copy), stored with the sensitivity report.
paper_c_covid19_source <- function() {
  list(
    repository = "https://github.com/uribo/japan-covid19",
    commit = "c4299b6",
    compiled = as.Date("2021-09-18"),
    copied = as.Date("2026-09-29")
  )
}
