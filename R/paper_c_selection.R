# NOTE: readable analysis code from the private project. It needs the fire-department
# transport records, which are not distributed, and is not sourced by the public
# pipeline (_targets.R). Kept so that the published estimates can be audited.

# Selection discordance between the fire-department unit and the prefecture unit
# (design note §8 ②, pre-specified 2026-09-10). Neither unit is treated as
# correct: what is reported is the two-way difference between the sets the two
# rules select, never an error rate of one against the other.
#
# Two selection rules, fixed before the comparison and never mixed:
#   (a) "rate"  - transports / observed department-days / 2020 population * 10,000
#   (b) "daily" - transports / observed department-days
# Raw counts are not ranked: the panel is unbalanced (Hokkaido holds 2022 only,
# Saitama has no 2022, and 639 department-days are absent), so a department's
# total depends on how many days it reported.

# Department-level table for one period. `years = NULL` keeps the primary sample.
build_selection_table <- function(panel, population, years = NULL) {
  data <- panel |>
    dplyr::mutate(year = lubridate::year(date))
  if (!is.null(years)) {
    data <- data |> dplyr::filter(year %in% years)
  }
  out <- data |>
    dplyr::summarise(
      transports = sum(transports_analysis, na.rm = TRUE),
      days = dplyr::n(),
      .by = c(prefecture, department)
    ) |>
    dplyr::left_join(
      population |> dplyr::select(prefecture, department, population),
      by = dplyr::join_by(prefecture, department)
    )
  missing <- sum(is.na(out$population))
  if (missing > 0L) {
    rlang::abort(paste0(missing, " department(s) carry no 2020 population"))
  }
  out |>
    dplyr::mutate(
      daily = transports / days,
      rate = daily / population * 10000
    )
}

# Accumulate units in descending order of `metric` until the cumulative
# population first reaches `target`, then include every unit tied with the
# boundary unit. Returns the whole table when the target is never reached.
accumulate_to <- function(data, metric, target) {
  ordered <- data |>
    dplyr::arrange(dplyr::desc(.data[[metric]])) |>
    dplyr::mutate(cumulative = cumsum(population))
  reached <- which(ordered$cumulative >= target)
  if (length(reached) == 0L) {
    return(list(selected = ordered, boundary = NA_real_, reached = FALSE))
  }
  boundary <- ordered[[metric]][reached[1]]
  list(
    selected = ordered |> dplyr::filter(.data[[metric]] >= boundary),
    boundary = boundary,
    reached = TRUE
  )
}

# Prefecture aggregate: the daily load is the sum of its departments' daily
# loads and the rate is that sum over the summed population. A prefecture here
# is the aggregate of the departments observed in the panel, not the whole
# prefecture.
aggregate_prefectures <- function(data) {
  data |>
    dplyr::summarise(
      daily = sum(daily),
      population = sum(population),
      .by = prefecture
    ) |>
    dplyr::mutate(rate = daily / population * 10000)
}

# One comparison: departments picked directly, against the departments that sit
# inside the prefectures picked to cover the same population. The boundary
# prefecture (with its ties) is reported both included and excluded; the two
# versions bracket the arbitrariness of an indivisible unit and are not a
# confidence interval.
compare_selection_units <- function(
  data,
  metric = c("rate", "daily"),
  share = 0.2
) {
  metric <- match.arg(metric)
  total_population <- sum(data$population)
  departments <- accumulate_to(data, metric, share * total_population)$selected
  selected_population <- sum(departments$population)

  prefectures <- aggregate_prefectures(data)
  boundary <- accumulate_to(prefectures, metric, selected_population)
  versions <- list(
    include = prefectures |>
      dplyr::filter(
        is.na(boundary$boundary) | .data[[metric]] >= boundary$boundary
      ),
    exclude = prefectures |>
      dplyr::filter(
        !is.na(boundary$boundary) & .data[[metric]] > boundary$boundary
      )
  )

  key <- function(x) paste(x$prefecture, x$department)
  summaries <- purrr::imap(versions, function(picked, label) {
    inside <- data |> dplyr::filter(prefecture %in% picked$prefecture)
    tibble::tibble(
      boundary_prefecture = label,
      n_prefectures = nrow(picked),
      n_departments_in_prefectures = nrow(inside),
      only_department_rule = sum(!key(departments) %in% key(inside)),
      only_prefecture_rule = sum(!key(inside) %in% key(departments)),
      population_difference = sum(picked$population) - selected_population
    )
  }) |>
    purrr::list_rbind() |>
    dplyr::mutate(
      discordance = only_department_rule + only_prefecture_rule,
      .after = only_prefecture_rule
    )

  list(
    metric = metric,
    share = share,
    n_departments_total = nrow(data),
    n_departments_selected = nrow(departments),
    population_total = total_population,
    population_selected = selected_population,
    departments = departments,
    versions = summaries
  )
}
