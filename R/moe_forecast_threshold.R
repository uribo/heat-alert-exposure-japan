# Forecast-threshold design of paper C (secondary analysis,
# notes/paper_strategy.md 6.9.1). The running variable is the area maximum of
# the daily-maximum forecast WBGT from the bulletins behind the alert rule;
# the reported quantities are the coefficient path over its integer values
# and the excess step at the alert threshold.

# Alert area of each forecast station. Stations in the alert-area assignment
# (2023 roster) keep their area; stations the roster lacks take the area of
# their prefecture code when that code maps to one area (41171 is 栃木;
# 6.7 notes it leaves the series in 2022-08). A code split across areas
# aborts, since the station could not be placed.
moe_forecast_station_area <- function(daily, assignment) {
  unplaced <- daily |>
    dplyr::distinct(wbgt_no) |>
    dplyr::mutate(pref_cd = substr(wbgt_no, 1L, 2L)) |>
    dplyr::anti_join(assignment, by = c(wbgt_no = "station_no"))
  code_area <- assignment |>
    dplyr::distinct(pref_cd = substr(station_no, 1L, 2L), alert_area) |>
    dplyr::filter(dplyr::n() == 1L, .by = pref_cd)
  station_area <- dplyr::bind_rows(
    dplyr::select(assignment, station_no, alert_area),
    unplaced |>
      dplyr::inner_join(code_area, by = "pref_cd") |>
      dplyr::select(station_no = wbgt_no, alert_area)
  )
  if (!all(unplaced$wbgt_no %in% station_area$station_no)) {
    rlang::abort("A forecast station has no alert area.")
  }
  station_area
}

# Area-day running variable: x_am and x_pm are the area maxima of the
# station daily maxima in the morning bulletin of D and the afternoon
# bulletin of D-1, and x the larger of the two. The archive holds no
# bulletin issued 01:00-07:00 on 2022-08-14, so that day has x_pm only and x
# is left NA (6.9.1: undefined), so that no consumer can use it by accident;
# any other missing bulletin aborts.
moe_forecast_area_day <- function(
  daily,
  station_area,
  areas,
  gap_dates = as.Date("2022-08-14")
) {
  area_day <- daily |>
    dplyr::inner_join(station_area, by = c(wbgt_no = "station_no")) |>
    dplyr::filter(alert_area %in% areas) |>
    dplyr::mutate(bulletin = dplyr::if_else(issue_hour >= 12L, "pm", "am")) |>
    dplyr::summarise(
      x = max(forecast_max),
      .by = c(alert_area, target_date, bulletin)
    ) |>
    tidyr::pivot_wider(
      id_cols = c(alert_area, target_date),
      names_from = bulletin,
      values_from = x,
      names_prefix = "x_"
    ) |>
    dplyr::mutate(x = pmax(x_pm, x_am))
  gaps <- unique(area_day$target_date[is.na(area_day$x_am)])
  if (anyNA(area_day$x_pm) || !all(gaps %in% gap_dates)) {
    rlang::abort(
      "An area-day lacks a bulletin outside the known 2022-08-14 gap."
    )
  }
  area_day
}

# 0/1 columns for the values of x, with x = 32 as the reference. Values at or
# below 29 and at or above 35 are pooled (35 and 36+ hold 27 and 9 area-days
# in 2021-2022). Returns the design with the columns and, as an attribute,
# the terms that vary in it (a subsample without x <= 29 has no x_le29).
forecast_path_levels <- c("le29", "30", "31", "32", "33", "34", "ge35")

add_forecast_path_terms <- function(design, x = "x") {
  value <- dplyr::case_when(
    design[[x]] <= 29 ~ "le29",
    design[[x]] >= 35 ~ "ge35",
    .default = as.character(design[[x]])
  )
  levels <- setdiff(forecast_path_levels, "32")
  for (level in levels) {
    design[[paste0("x_", level)]] <- as.integer(value == level)
  }
  terms <- paste0("x_", levels)
  attr(design, "path_terms") <- terms[purrr::map_lgl(terms, \(term) {
    dplyr::n_distinct(design[[term]]) > 1L
  })]
  design
}

# With b32 = 0: step32 = -b31, step33 = b33, step34 = b34 - b33. The primary
# contrast is the excess step C = step33 - (step32 + step34) / 2.
forecast_path_contrasts <- function(terms) {
  weights <- list(
    "excess step: step33 - mean(step32, step34)" = c(
      x_31 = 0.5,
      x_33 = 1.5,
      x_34 = -0.5
    ),
    "step33 - step32" = c(x_31 = 1, x_33 = 1),
    "step33 - step34" = c(x_33 = 2, x_34 = -1)
  )
  unit <- purrr::map(terms, \(term) stats::setNames(1, term)) |>
    stats::setNames(terms)
  purrr::map(c(unit, weights), \(w) {
    if (!all(names(w) %in% terms)) {
      rlang::abort("A contrast uses a value of x absent from the fit.")
    }
    a <- stats::setNames(numeric(length(terms)), terms)
    a[names(w)] <- w
    a
  })
}

# Critical value of a bias-adjusted interval: the 1 - alpha quantile of
# |T + b| with T ~ t(df) and b = M / se, the worst-case bias in standard
# errors. This borrows the fixed-length construction of Kolesár and Rothe
# (2018), but their coverage guarantee rests on their conditions and a normal
# approximation; with a Poisson fixed-effects fit, the CV3J se and t(G - 1)
# (chosen so that M = 0 reproduces the CV3J interval of 7.13) the interval is
# an approximation whose coverage in this design has not been checked.
bias_adjusted_critical_value <- function(bias_ratio, df, level = 0.95) {
  purrr::map_dbl(bias_ratio, \(b) {
    stats::uniroot(
      \(c) stats::pt(c - b, df) - stats::pt(-c - b, df) - level,
      lower = 0,
      upper = b + stats::qt(1 - (1 - level) / 2, df) + 10
    )$root
  })
}

# Bias-adjusted intervals, under an assumed curvature bound M, for the
# contrasts whose worst-case bias under |second difference of the smooth
# path| <= M is M (6.9.1): the excess step and the two one-sided step
# differences.
forecast_bias_adjusted_intervals <- function(cv3, m_grid, df) {
  cv3 |>
    dplyr::filter(!startsWith(contrast, "x_")) |>
    tidyr::expand_grid(m = m_grid) |>
    dplyr::mutate(
      critical_value = bias_adjusted_critical_value(m / se_cv3, df),
      adjusted_low = estimate - critical_value * se_cv3,
      adjusted_high = estimate + critical_value * se_cv3
    )
}

# Fit the path on one design, with CR1 and (when `jackknife` is TRUE) CV3J
# inference on every coefficient and contrast.
fit_forecast_path <- function(design, label, jackknife = TRUE, workers = 1L) {
  terms <- attr(design, "path_terms")
  fit <- fit_paper_c_power(design, paper_c_result2_formula(terms, "0.1"))
  kept <- design[fit_kept_rows(fit, nrow(design)), ]
  n_clusters <- dplyr::n_distinct(kept$area_primary)
  b <- stats::coef(fit)[terms]
  v <- stats::vcov(fit)[terms, terms]
  contrasts <- forecast_path_contrasts(terms)
  critical_cr1 <- stats::qt(0.975, df = n_clusters - 1L)
  cr1 <- purrr::imap(contrasts, \(a, name) {
    estimate <- sum(a * b)
    se <- sqrt(drop(t(a) %*% v %*% a))
    tibble::tibble(
      sample = label,
      contrast = name,
      estimate = estimate,
      se_cr1 = se,
      conf_low_cr1 = estimate - critical_cr1 * se,
      conf_high_cr1 = estimate + critical_cr1 * se
    )
  }) |>
    purrr::list_rbind()
  out <- dplyr::mutate(cr1, n = stats::nobs(fit), clusters = n_clusters)
  if (jackknife) {
    jk <- paper_c_jackknife_vcov(kept, terms, "0.1", workers = workers)
    out <- dplyr::left_join(
      out,
      summarise_jackknife(jk, b, contrasts),
      by = "contrast"
    )
    attr(out, "jackknife") <- jk$estimates
  }
  out
}

# Input chain (Issue #50 (a)) ------------------------------------------------------

# One frozen bulletin set: the manifest, the bodies it lists (read_moe_forecast()
# checks each against its sha256 and its bulletin hour) and the station daily
# maxima. The bodies are gitignored; on a machine without them, copy
# OneDrive _frozen/23_Heatstroke/moe_forecast/ and check it with shasum first.
read_moe_forecast_set <- function(
  manifest_path,
  expected_files,
  dir = here::here("data-raw/frozen/moe_forecast")
) {
  manifest <- readr::read_tsv(
    manifest_path,
    col_types = readr::cols(
      from_date = readr::col_date(),
      to_date = readr::col_date(),
      records = readr::col_integer(),
      bytes = readr::col_integer(),
      .default = readr::col_character()
    )
  )
  if (nrow(manifest) != expected_files) {
    rlang::abort(glue::glue(
      "{basename(manifest_path)}: {nrow(manifest)} files, expected {expected_files}."
    ))
  }
  read_moe_forecast(manifest, dir) |>
    moe_forecast_daily_max()
}

# Area-day running variable for each bulletin set, on the released area roster,
# joined to the alert record. `daily_sets` is a named list of station
# daily maxima; the name becomes the `set` column.
build_moe_forecast_area_day <- function(
  daily_sets,
  wbgt_stations,
  heat_alert,
  alert_records,
  area_roster
) {
  assignment <- assign_wbgt_station_alert_area(
    wbgt_stations,
    heat_alert
  )$assignment
  roster_areas <- unique(area_roster$area_primary)
  alert_days <- make_moe_alert_days(alert_records, heat_alert$area) |>
    dplyr::rename(alert_area = area, target_date = date)
  purrr::imap(daily_sets, \(daily, set) {
    station_area <- moe_forecast_station_area(daily, assignment)
    area_day <- moe_forecast_area_day(daily, station_area, roster_areas) |>
      dplyr::inner_join(alert_days, by = c("alert_area", "target_date"))
    if (nrow(area_day) != length(roster_areas) * 2L * 122L) {
      rlang::abort(glue::glue(
        "{set}: expected {length(roster_areas) * 2L * 122L} area-days, ",
        "got {nrow(area_day)}."
      ))
    }
    dplyr::mutate(area_day, set = set, .before = 1L)
  }) |>
    purrr::list_rbind()
}

# Report (Issue #50 (a)) -----------------------------------------------------------

# The trial's coefficient path (tracked; copied from the gitignored trial
# report), read and checked before any fit.
read_moe_forecast_trial_path <- function(path) {
  trial_path <- readr::read_tsv(
    path,
    comment = "#",
    col_types = readr::cols(term = "c", estimate = "d")
  )
  if (
    nrow(trial_path) != 6L ||
      !setequal(trial_path$term, c("<=29", "30", "31", "33", "34", ">=35")) ||
      !all(is.finite(trial_path$estimate))
  ) {
    rlang::abort(
      "notes/moe_forecast_trial_path.tsv does not hold the six trial terms."
    )
  }
  trial_path
}

# The counts of the 2021-2022 sample fixed in 6.9.1 on the pre-specified 23
# prefectures. They are the default so that src/ drivers keep their meaning;
# the target passes the counts of the current sample (Issue #59).
moe_forecast_expected_23 <- function() {
  list(areas = 30L, both_sides = 22L, n = 81404L)
}

# The 25-prefecture sample (Issue #59 A, 2026-10-05): 富山, 沖縄本島, 宮古島
# and 八重山 add four areas, all observed on both sides of 32 | 33, and 5,900
# department-days to the fit.
moe_forecast_expected_25 <- function() {
  list(areas = 34L, both_sides = 26L, n = 87304L)
}

# The secondary analysis fixed in notes/paper_strategy.md 6.9.1 on 2026-09-30.
# Primary: the 05:00/17:00 bulletins, the coefficient path on all 2021-2022
# department-days, the excess step C with CR1, CV3J and bias-adjusted
# intervals under an assumed curvature bound M (approximate; coverage in this
# design unchecked). Robustness: the path on x in [30, 35]. Sensitivity: the
# 04:00/16:00 bulletins. By year: CR1 only. `design_full` is the full sample
# of build_paper_c_analysis_design(); the fits are deterministic, so the
# worker count does not change the report. The path on the pre-specified 23
# prefectures (CR1 only) must reproduce the 2026-09-29 trial and is kept as
# `prespecified` for the disclosure of the move to 25 (Issue #59).
build_moe_forecast_threshold_report <- function(
  design_full,
  area_day_sets,
  trial_path,
  m_grid = c(0, 0.02, 0.05, 0.10),
  expected = moe_forecast_expected_23(),
  workers = 1L
) {
  area_days <- split(
    dplyr::select(area_day_sets, !set),
    area_day_sets$set
  )[c("05/17", "04/16")]

  # Agreement of each set with the alert record, by bulletin and for the day.
  first_stage <- purrr::imap(area_days, \(area_day, set) {
    dplyr::bind_rows(
      dplyr::count(area_day, rule = x_pm >= 33, alert = alert_17pm) |>
        dplyr::mutate(level = "afternoon bulletin (D-1)"),
      dplyr::count(
        dplyr::filter(area_day, !is.na(x_am)),
        rule = x_am >= 33,
        alert = alert_5am
      ) |>
        dplyr::mutate(level = "morning bulletin (D)"),
      # The rule as reconstructed in 6.9: on 2022-08-14, where x is undefined,
      # the afternoon bulletin alone decides.
      dplyr::count(
        area_day,
        rule = pmax(x_pm, x_am, na.rm = TRUE) >= 33,
        alert = alert_day
      ) |>
        dplyr::mutate(level = "day")
    ) |>
      dplyr::mutate(set = set)
  }) |>
    purrr::list_rbind()

  # How far the 04:00/16:00 bulletins move x, value by value.
  x_crosstab <- dplyr::inner_join(
    dplyr::select(area_days[["05/17"]], alert_area, target_date, x_0517 = x),
    dplyr::select(area_days[["04/16"]], alert_area, target_date, x_0416 = x),
    by = c("alert_area", "target_date")
  ) |>
    dplyr::filter(target_date != as.Date("2022-08-14")) |>
    dplyr::count(x_0517, x_0416)

  design <- design_full |>
    dplyr::filter(lubridate::year(date) %in% 2021:2022) |>
    dplyr::left_join(
      dplyr::select(area_days[["05/17"]], alert_area, target_date, x),
      by = c(area_primary = "alert_area", date = "target_date")
    ) |>
    dplyr::left_join(
      dplyr::select(area_days[["04/16"]], alert_area, target_date, x_0416 = x),
      by = c(area_primary = "alert_area", date = "target_date")
    ) |>
    # x is undefined without the morning bulletin on 2022-08-14 (6.9.1).
    dplyr::filter(date != as.Date("2022-08-14"))
  if (anyNA(design$x) || anyNA(design$x_0416)) {
    rlang::abort(
      "A department-day in 2021-2022 has no forecast running variable."
    )
  }

  count_both_sides <- function(design) {
    design |>
      dplyr::filter(x %in% c(32, 33)) |>
      dplyr::summarise(sides = dplyr::n_distinct(x), .by = area_primary) |>
      dplyr::count(sides)
  }
  both_sides <- count_both_sides(design)
  # The fits below check n after fixest's drops.
  if (
    dplyr::n_distinct(design$area_primary) != expected$areas ||
      both_sides$n[both_sides$sides == 2L] != expected$both_sides
  ) {
    rlang::abort(glue::glue(
      "The 2021-2022 sample no longer has {expected$areas} areas, ",
      "{expected$both_sides} on both sides."
    ))
  }

  primary <- fit_forecast_path(
    add_forecast_path_terms(design),
    "primary (05/17)",
    workers = workers
  )
  if (primary$n[1] != expected$n) {
    rlang::abort(glue::glue(
      "{primary$sample[1]}: n = {primary$n[1]}, expected {expected$n}."
    ))
  }

  # The path on the pre-specified 23 prefectures must reproduce the trial's
  # (same rows, same model, dummies instead of i()); a difference means the
  # inputs moved. Its counts are those of 6.9.1 and do not follow `expected`.
  design_23 <- dplyr::filter(
    design,
    prefecture %in% paper_c_prespecified_prefectures()
  )
  both_sides_23 <- count_both_sides(design_23)
  if (
    dplyr::n_distinct(design_23$area_primary) != 30L ||
      both_sides_23$n[both_sides_23$sides == 2L] != 22L
  ) {
    rlang::abort(
      "The pre-specified 23 prefectures no longer give 30 areas, 22 on both sides."
    )
  }
  prespecified <- fit_forecast_path(
    add_forecast_path_terms(design_23),
    "pre-specified 23 prefectures (05/17)",
    jackknife = FALSE
  )
  if (prespecified$n[1] != 81404L) {
    rlang::abort(glue::glue(
      "{prespecified$sample[1]}: n = {prespecified$n[1]}, expected 81404."
    ))
  }
  reproduced <- prespecified |>
    dplyr::filter(startsWith(contrast, "x_")) |>
    dplyr::mutate(
      term = dplyr::case_when(
        contrast == "x_le29" ~ "<=29",
        contrast == "x_ge35" ~ ">=35",
        .default = sub("^x_", "", contrast)
      )
    ) |>
    dplyr::inner_join(trial_path, by = "term", suffix = c("", "_trial"))
  if (
    nrow(reproduced) != 6L ||
      max(abs(reproduced$estimate - reproduced$estimate_trial)) > 1e-8
  ) {
    rlang::abort(
      "The pre-specified path does not reproduce the 2026-09-29 trial."
    )
  }

  window <- design |>
    dplyr::filter(dplyr::between(x, 30, 35)) |>
    add_forecast_path_terms() |>
    fit_forecast_path("robustness: x in [30, 35]", workers = workers)

  sensitivity <- design |>
    dplyr::mutate(x = x_0416) |>
    add_forecast_path_terms() |>
    fit_forecast_path("sensitivity: 04/16 bulletins", workers = workers)

  by_year <- purrr::map(2021:2022, \(y) {
    design |>
      dplyr::filter(lubridate::year(date) == y) |>
      add_forecast_path_terms() |>
      fit_forecast_path(paste("year", y), jackknife = FALSE)
  }) |>
    purrr::list_rbind()

  bias_adjusted <- purrr::map(list(primary, window, sensitivity), \(r) {
    forecast_bias_adjusted_intervals(r, m_grid, df = r$clusters[1] - 1L)
  }) |>
    purrr::list_rbind()

  list(
    specification = "notes/paper_strategy.md 6.9.1 (fixed 2026-09-30)",
    expected = expected,
    first_stage = first_stage,
    x_crosstab = x_crosstab,
    both_sides = both_sides,
    results = dplyr::bind_rows(primary, window, sensitivity, by_year),
    prespecified = prespecified,
    bias_adjusted = bias_adjusted,
    jackknife = list(
      primary = attr(primary, "jackknife"),
      window = attr(window, "jackknife"),
      sensitivity = attr(sensitivity, "jackknife")
    ),
    m_grid = m_grid,
    fixest_version = as.character(utils::packageVersion("fixest"))
  )
}
