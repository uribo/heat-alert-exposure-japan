# NOTE: readable analysis code from the private project. It needs the fire-department
# transport records, which are not distributed, and is not sourced by the public
# pipeline (_targets.R). Kept so that the published estimates can be audited.

# Reproduce the independent step-3 panel without reading or rebuilding targets.
build_paper_c_panel <- function(
  df_hs_japan,
  df_fd_alert_area,
  df_moe_wbgt_daily_agg,
  df_moe_alert_historical,
  sf_heat_alert,
  expected
) {
  alert_days <- make_moe_alert_days(
    df_moe_alert_historical,
    sf::st_drop_geometry(sf_heat_alert)$area
  )
  stationless <- df_fd_alert_area |>
    dplyr::filter(is.na(near_wbgt_station)) |>
    dplyr::select(prefecture, department)

  attach_exposure <- function(
    data,
    assignment = df_fd_alert_area,
    allow_stationless = stationless
  ) {
    joined <- data |>
      dplyr::mutate(
        year = lubridate::year(date),
        month = lubridate::month(date)
      ) |>
      dplyr::left_join(
        assignment,
        by = dplyr::join_by(prefecture, department),
        relationship = "many-to-one"
      )
    # Count after the join, not before: a department with no assignment would
    # otherwise lose its area, its station, its WBGT and finally every row.
    unmatched <- joined |>
      dplyr::filter(is.na(area_primary) | is.na(near_wbgt_station)) |>
      dplyr::distinct(prefecture, department) |>
      dplyr::anti_join(
        allow_stationless,
        by = dplyr::join_by(prefecture, department)
      )
    if (nrow(unmatched) > 0L) {
      rlang::abort(paste0(
        "No alert area or station for ",
        nrow(unmatched),
        " undeclared department(s): ",
        paste(unmatched$department, collapse = ", ")
      ))
    }
    joined |>
      dplyr::left_join(
        df_moe_wbgt_daily_agg |>
          dplyr::transmute(
            station_no,
            date,
            wbgt_max = dplyr::if_else(
              is.finite(wbgt_max),
              wbgt_max,
              NA_real_
            )
          ),
        by = dplyr::join_by(near_wbgt_station == station_no, date),
        relationship = "many-to-one",
        # No department is stationless since the crosswalk repair, but an NA
        # key must never match if one appears again.
        na_matches = "never"
      ) |>
      dplyr::left_join(
        alert_days,
        by = dplyr::join_by(area_primary == area, date),
        relationship = "many-to-one",
        na_matches = "never"
      ) |>
      dplyr::mutate(
        alert_day = tidyr::replace_na(alert_day, FALSE),
        alert_5am = tidyr::replace_na(alert_5am, FALSE),
        alert_17pm = tidyr::replace_na(alert_17pm, FALSE),
        local_high = wbgt_max >= 33,
        # Declared imputation. All NA transports in this window are rows where the
        # whole record is blank (dispatches is NA too), which is what the export
        # produces for a zero-transport day rather than for an absent unit. That
        # reading is the working assumption of paper_strategy 1.3 and is still
        # pending confirmation against one raw file per prefecture (TODO.md), so
        # the count is reported rather than buried.
        transports_analysis = tidyr::replace_na(transports, 0)
      )
  }

  # The 2020 trial region (1都8県) is taken from the areas with a 2020 record
  # in the issuance history. Pin it, so a refetched history that carries 2020
  # rows for other areas cannot silently widen the trial.
  trial_areas <- df_moe_alert_historical |>
    dplyr::filter(lubridate::year(datetime) == 2020L) |>
    dplyr::distinct(area) |>
    dplyr::pull(area)
  expected_trial <- c(
    "茨城",
    "栃木",
    "群馬",
    "埼玉",
    "千葉",
    "東京",
    "神奈川",
    "山梨",
    "長野"
  )
  if (!setequal(trial_areas, expected_trial)) {
    rlang::abort(paste(
      "2020 trial areas changed:",
      paste(sort(trial_areas), collapse = ", ")
    ))
  }

  panel <- df_hs_japan |>
    dplyr::filter(
      lubridate::year(date) %in% 2018:2022,
      lubridate::month(date) %in% 6:9,
      prefecture != "千葉県"
    ) |>
    attach_exposure() |>
    dplyr::mutate(
      system_active = year >= 2021L |
        (year == 2020L &
          date >= as.Date("2020-07-01") &
          area_primary %in% trial_areas),
      alert_status = dplyr::case_when(
        alert_day ~ "alert",
        system_active ~ "active, no alert",
        .default = "structural no-alert"
      )
    )

  stopifnot(
    nrow(stationless) == 0L,
    !anyNA(panel$area_primary),
    !any(panel$alert_day & !panel$system_active),
    # The imputation above assumes every NA transports row is a fully blank
    # record; outside this window df_hs_japan has 23,372 rows with dispatches
    # present and transports NA, so the assumption is asserted, not trusted.
    !any(is.na(panel$transports) & !is.na(panel$dispatches))
  )
  # The 2020 status split is what separates the trial from the structural
  # no-alert rows; the four-cell guard keys on alert_day only and would not
  # notice a wrong trial classification.
  status_2020 <- panel |>
    dplyr::filter(year == 2020L) |>
    dplyr::count(alert_status) |>
    dplyr::arrange(alert_status)
  cells <- panel |>
    dplyr::filter(!is.na(local_high)) |>
    dplyr::count(alert_day, local_high) |>
    dplyr::arrange(alert_day, local_high)
  # The counts are fixed by the caller (`expected`), so a change of sample is
  # a change of the target definition. History of the 23-prefecture values:
  # the step-3 cells were 179,410 / 5,424 / 8,132 / 1,588 with the
  # distance-only station assignment; the coverage-aware assignment (issue
  # #22) moved two departments' stations and added their 488 previously
  # station-less days; repairing the Wikipedia 管轄区域 cells (2026-09-08) gave
  # the last ten departments a geometry and brought their 5,119 FD-days into
  # the cells, WBGT coverage growing from 195,042 to 200,161 of 200,173 rows.
  observed <- list(
    rows = nrow(panel),
    departments = dplyr::n_distinct(panel$department),
    areas = dplyr::n_distinct(panel$area_primary),
    blank_transports = sum(is.na(panel$transports)),
    transports = sum(panel$transports_analysis),
    status_2020 = stats::setNames(status_2020$n, status_2020$alert_status),
    cells = cells$n,
    alert_2020 = sum(panel$year == 2020L & panel$alert_day)
  )
  check_paper_c_panel_counts(observed, expected)
  message(
    "Panel guards passed: ",
    observed$rows,
    " rows; four cells ",
    paste(observed$cells, collapse = "/"),
    "; 2020 alert FD-days: ",
    observed$alert_2020
  )
  panel
}

# Compare the panel's counts with the values the target fixes. Every field
# is compared as numbers with names kept, and a mismatch reports all observed
# values at once so that a change of sample can be checked in one run.
check_paper_c_panel_counts <- function(observed, expected) {
  fields <- c(
    "rows",
    "departments",
    "areas",
    "blank_transports",
    "transports",
    "status_2020",
    "cells",
    "alert_2020"
  )
  if (!setequal(names(expected), fields)) {
    rlang::abort(paste(
      "`expected` must name exactly:",
      paste(fields, collapse = ", ")
    ))
  }
  differs <- fields[
    !vapply(
      fields,
      \(field) {
        identical(
          stats::setNames(
            as.numeric(observed[[field]]),
            names(observed[[field]])
          ),
          stats::setNames(
            as.numeric(expected[[field]]),
            names(expected[[field]])
          )
        )
      },
      logical(1)
    )
  ]
  if (length(differs) > 0L) {
    shown <- vapply(
      fields,
      \(field) {
        value <- observed[[field]]
        paste0(
          field,
          " = ",
          paste(
            if (is.null(names(value))) {
              value
            } else {
              paste0(names(value), ":", value)
            },
            collapse = "/"
          )
        )
      },
      character(1)
    )
    rlang::abort(c(
      paste("Paper-C panel counts differ in:", paste(differs, collapse = ", ")),
      i = paste("Observed:", paste(shown, collapse = "; "))
    ))
  }
  invisible(observed)
}
