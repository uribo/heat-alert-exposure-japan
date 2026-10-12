# NOTE: readable analysis code from the private project. It needs the fire-department
# transport records, which are not distributed, and is not sourced by the public
# pipeline (_targets.R). Kept so that the published estimates can be audited.

# Result 1 of paper C (notes/paper_strategy.md 7; steps 5a-5b): concordance of
# the alert with the realized local WBGT, department onset and slope profiles,
# their dispersion within and between areas, and the descriptive covariate
# regressions. Moved from the chunks of notes/paper_c_result1.qmd (Issue #50
# (a)); the note now only displays the objects returned here. The code is the
# chunks' code, unchanged apart from the function wrapper.

# The complementary GLMM, fitted in a child process by callr::r() so that the
# 30-minute bound holds even during native optimization. Defined at top level
# so that the child receives the function without its caller's environment.
fit_paper_c_result1_glmm <- function(data) {
  conditions <- character()
  fit <- withCallingHandlers(
    lme4::glmer(
      transports_analysis ~ hinge +
        (1 + hinge | department) +
        (1 + hinge | area_primary) +
        factor(year),
      family = stats::poisson(),
      data = data
    ),
    warning = function(w) {
      conditions <<- c(conditions, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      conditions <<- c(conditions, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  list(
    variance = as.data.frame(lme4::VarCorr(fit)),
    fixed_slope = lme4::fixef(fit)[["hinge"]],
    fixed_slope_se = sqrt(diag(as.matrix(vcov(fit))))[["hinge"]],
    singular = lme4::isSingular(fit),
    optimizer_code = fit@optinfo$conv$opt,
    convergence_messages = fit@optinfo$conv$lme4$messages,
    conditions = unique(conditions)
  )
}

# Everything notes/paper_c_result1.qmd displays, plus `report`, the aggregate
# list draft/paper_c.qmd reads (no raw series, no Saitama names).
build_paper_c_result1 <- function(
  df_paper_c_panel,
  df_fd_alert_area,
  df_fd_merge_pop2020_sum,
  sf_fd_merge_all
) {
  # Assignments, exposure, and cell shares --------------------------------
  assignments <- df_fd_alert_area |>
    dplyr::arrange(prefecture, department) |>
    dplyr::mutate(local_id = dplyr::row_number(), .by = prefecture) |>
    dplyr::mutate(
      department_display = dplyr::if_else(
        prefecture == "埼玉県",
        sprintf("Saitama FD %02d", local_id),
        department
      )
    )
  stationed <- assignments |> dplyr::filter(!is.na(near_wbgt_station))
  stopifnot(nrow(stationed) == 417L)
  # Keep every joined department and expose failed matches before any analysis.
  join_audit <- tibble::tibble(
    join = character(),
    departments = integer(),
    unmatched = integer()
  )
  audit_join <- function(x, y, by, label) {
    y$.joined_marker <- TRUE
    out <- dplyr::left_join(x, y, by = by, relationship = "many-to-one")
    failures <- out |>
      dplyr::filter(is.na(.joined_marker)) |>
      dplyr::distinct(prefecture, department)
    join_audit <<- dplyr::bind_rows(
      join_audit,
      tibble::tibble(
        join = label,
        departments = nrow(dplyr::distinct(x, prefecture, department)),
        unmatched = nrow(failures)
      )
    )
    message(label, ": unmatched departments = ", nrow(failures))
    dplyr::select(out, -.joined_marker)
  }
  panel <- df_paper_c_panel |> dplyr::filter(!is.na(local_high))
  stopifnot(nrow(panel) == 215613L)
  # A department with a station can still have no WBGT-complete day: its station
  # may carry no data in its only season. Such departments are named here and
  # excluded from every department-level table, so the per-department joins
  # below count them as declared rather than as silent NA rows.
  no_exposure <- stationed |>
    dplyr::anti_join(
      dplyr::distinct(panel, prefecture, department),
      by = dplyr::join_by(prefecture, department)
    ) |>
    dplyr::left_join(
      df_paper_c_panel |>
        dplyr::summarise(
          fd_days = dplyr::n(),
          wbgt_missing = sum(is.na(wbgt_max)),
          years = paste(sort(unique(year)), collapse = ", "),
          .by = c(prefecture, department)
        ),
      by = dplyr::join_by(prefecture, department)
    )
  exposed <- stationed |>
    dplyr::anti_join(no_exposure, by = dplyr::join_by(prefecture, department))
  stopifnot(nrow(exposed) + nrow(no_exposure) == 417L)
  statuses <- c("alert", "active, no alert", "structural no-alert")
  safe_share <- function(n, d) {
    n / dplyr::if_else(d > 0, as.numeric(d), NA_real_)
  }
  cell_summary <- function(data, groups = character()) {
    data |>
      dplyr::summarise(
        fd_days = dplyr::n(),
        total_transports = sum(transports_analysis),
        .by = dplyr::all_of(c(groups, "alert_status", "local_high"))
      )
  }
  add_shares <- function(data, groups = character()) {
    data |>
      dplyr::mutate(
        day_row_share = safe_share(fd_days, sum(fd_days)),
        transport_row_share = safe_share(
          total_transports,
          sum(total_transports)
        ),
        .by = dplyr::all_of(c(groups, "alert_status"))
      ) |>
      dplyr::mutate(
        day_column_share = safe_share(fd_days, sum(fd_days)),
        transport_column_share = safe_share(
          total_transports,
          sum(total_transports)
        ),
        .by = dplyr::all_of(c(groups, "local_high"))
      ) |>
      dplyr::mutate(
        day_share = safe_share(fd_days, sum(fd_days)),
        transport_share = safe_share(total_transports, sum(total_transports)),
        transports_per_day = safe_share(total_transports, fd_days),
        .by = dplyr::all_of(groups)
      )
  }
  cells <- cell_summary(panel) |>
    tidyr::complete(
      alert_status = statuses,
      local_high = c(FALSE, TRUE),
      fill = list(fd_days = 0L, total_transports = 0)
    ) |>
    add_shares()
  area_year <- cell_summary(panel, c("area_primary", "year")) |>
    tidyr::complete(
      area_primary = unique(assignments$area_primary),
      year = 2018:2022,
      alert_status = statuses,
      local_high = c(FALSE, TRUE),
      fill = list(fd_days = 0L, total_transports = 0)
    ) |>
    add_shares(c("area_primary", "year"))
  year_cells <- cell_summary(panel, "year") |>
    tidyr::complete(
      year = 2018:2022,
      alert_status = statuses,
      local_high = c(FALSE, TRUE),
      fill = list(fd_days = 0L, total_transports = 0)
    ) |>
    add_shares("year")
  active_shares <- panel |>
    dplyr::summarise(
      active_days = sum(system_active),
      joint_days = sum(system_active & alert_day & local_high),
      alert_days = sum(system_active & alert_day),
      high_days = sum(system_active & local_high),
      .by = c(prefecture, department)
    )
  shares <- audit_join(
    exposed,
    active_shares,
    dplyr::join_by(prefecture, department),
    "Active-day denominators"
  ) |>
    dplyr::mutate(
      share_high_given_alert = safe_share(joint_days, alert_days),
      share_alert_given_high = safe_share(joint_days, high_days)
    )

  # Covariates ------------------------------------------------------------
  # Two earlier gaps are both closed. Seven departments once carried population
  # 0 with prop_age_class3 NaN, because the upstream sum ran na.rm = TRUE over an
  # all-NA join; the 2026-09-08 crosswalk repair ended that. Eight more lacked
  # population entirely - the six sharing their only municipality with another
  # department, plus the two abolished 旧) departments - until the apportionment
  # of 2026-09-10 (issue #17) gave every one of the 391 a weighted population.
  # The zero-as-missing handling below is kept even so: it costs nothing, and it
  # is what would surface a regression in the crosswalk rather than letting a
  # department enter the covariates with a population of zero.
  population <- df_fd_merge_pop2020_sum |>
    dplyr::transmute(
      prefecture,
      department,
      population_zero = !is.na(population) & population == 0,
      population = dplyr::if_else(population_zero, NA_real_, population),
      prop_age_class3 = dplyr::if_else(
        population_zero,
        NA_real_,
        prop_age_class3
      )
    )
  all_population <- audit_join(
    assignments,
    population,
    dplyr::join_by(prefecture, department),
    "Population on 417 departments"
  )
  covariates <- audit_join(
    exposed,
    population,
    dplyr::join_by(prefecture, department),
    "Population on 417 departments"
  )
  geometry <- sf_fd_merge_all |> sf::st_make_valid() |> sf::st_transform(6933)
  message("Constructing nationwide union boundary")
  coast <- sf::st_boundary(sf::st_union(sf::st_geometry(geometry)))
  message("Measuring shared boundary lengths")
  shared <- vapply(
    seq_len(nrow(geometry)),
    function(i) {
      intersection <- sf::st_intersection(
        sf::st_boundary(sf::st_geometry(geometry)[i]),
        coast
      )
      sum(as.numeric(sf::st_length(intersection)))
    },
    numeric(1)
  )
  geometry_data <- geometry |>
    dplyr::mutate(
      area_km2 = as.numeric(sf::st_area(geometry)) / 1e6,
      shared_boundary_m = shared,
      coastal = shared > 0
    ) |>
    sf::st_drop_geometry() |>
    dplyr::select(prefecture, department, area_km2, shared_boundary_m, coastal)
  covariates <- audit_join(
    covariates,
    geometry_data,
    dplyr::join_by(prefecture, department),
    "Jurisdiction geometry on 417 departments"
  )
  # jmastats::stations repeats 37 station numbers (relocated or multi-type
  # stations, e.g. 帯広 at 38 m and 45 m), so collapse to one row per station
  # and keep a flag for the ones whose rows disagree on elevation.
  station_elevations <- jmastats::stations |>
    sf::st_drop_geometry() |>
    dplyr::summarise(
      station_elevation_m = as.numeric(dplyr::first(elevation)),
      station_elevation_conflict = dplyr::n_distinct(elevation) > 1L,
      .by = station_no
    ) |>
    dplyr::mutate(
      near_wbgt_station = as.character(station_no),
      .keep = "unused"
    )
  covariates <- audit_join(
    covariates,
    station_elevations,
    dplyr::join_by(near_wbgt_station),
    "Station elevation on 417 departments"
  ) |>
    dplyr::mutate(
      density = population / area_km2,
      log_density = log(density),
      log_distance = log(dist_m),
      jurisdiction_elevation_m = NA_real_
    )
  inventory_columns <- c(
    "dist_m",
    "prop_age_class3",
    "density",
    "coastal",
    "station_elevation_m",
    "jurisdiction_elevation_m"
  )
  inventory <- purrr::map(inventory_columns, function(column) {
    present <- !is.na(covariates[[column]]) & is.finite(covariates[[column]])
    tibble::tibble(
      covariate = column,
      available = sum(present),
      total = nrow(covariates),
      lacking = paste(covariates$department_display[!present], collapse = "; ")
    )
  }) |>
    purrr::list_rbind()
  covariate_join_audit <- join_audit

  # Stage one: department onset and slope profiles ------------------------
  profile_department <- function(data, grid = seq(24, 32, 0.5)) {
    warnings <- character()
    errors <- character()
    # factor(year) needs two or more observed years; 57 Hokkaido departments
    # have 2022 only, so the year term is dropped for single-year departments
    # rather than letting glm() fail on a one-level factor.
    multi_year <- dplyr::n_distinct(data$year) > 1L
    fits <- lapply(grid, function(k) {
      tryCatch(
        withCallingHandlers(
          # Build the formula here so that its environment holds k.
          stats::glm(
            if (multi_year) {
              transports_analysis ~ pmax(wbgt_max - k, 0) + factor(year)
            } else {
              transports_analysis ~ pmax(wbgt_max - k, 0)
            },
            family = stats::poisson(),
            data = data
          ),
          warning = function(w) {
            warnings <<- c(warnings, paste0("k=", k, ": ", conditionMessage(w)))
            invokeRestart("muffleWarning")
          }
        ),
        error = function(e) {
          errors <<- c(errors, paste0("k=", k, ": ", conditionMessage(e)))
          NULL
        }
      )
    })
    # A fit whose hinge coefficient is aliased (every day at or below k) has no
    # slope; summary() would silently drop the row and [2, 1] would then read
    # the first year dummy. Treat such fits as unavailable at that k.
    fits <- lapply(fits, function(x) {
      if (!is.null(x) && anyNA(stats::coef(x))) NULL else x
    })
    deviances <- vapply(
      fits,
      function(x) if (is.null(x)) Inf else x$deviance,
      numeric(1)
    )
    order_fit <- order(deviances, grid)
    best <- order_fit[1]
    valid <- is.finite(deviances[best])
    coefficient <- if (valid) summary(fits[[best]])$coefficients else NULL
    tibble::tibble(
      onset = if (valid) grid[best] else NA_real_,
      slope = if (valid) coefficient[2, 1] else NA_real_,
      se = if (valid) coefficient[2, 2] else NA_real_,
      deviance_gap = if (valid) {
        deviances[order_fit[2]] - deviances[best]
      } else {
        NA_real_
      },
      total_transports = sum(data$transports_analysis),
      days = nrow(data),
      years = dplyr::n_distinct(data$year),
      year_term = dplyr::n_distinct(data$year) > 1L,
      flat_profile = if (valid) {
        deviances[order_fit[2]] - deviances[best] < 1
      } else {
        NA
      },
      onset_at_bound = if (valid) grid[best] %in% range(grid) else NA,
      under_30 = sum(data$transports_analysis) < 30,
      glm_warnings = paste(unique(warnings), collapse = " | "),
      glm_errors = paste(unique(errors), collapse = " | "),
      converged = if (valid) fits[[best]]$converged else FALSE
    )
  }
  stage_one <- panel |>
    dplyr::group_by(prefecture, department) |>
    dplyr::group_modify(~ profile_department(.x)) |>
    dplyr::ungroup()
  # Sensitivity: most onsets sit on the 24 degree floor of the specified grid, so
  # the profile is repeated over 18-32 to see whether the floor, not the data,
  # is choosing the onset. Only the onset-related columns are kept.
  stage_one_wide <- panel |>
    dplyr::group_by(prefecture, department) |>
    dplyr::group_modify(~ profile_department(.x, grid = seq(18, 32, 0.5))) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      prefecture,
      department,
      onset_wide = onset,
      slope_wide = slope,
      onset_wide_at_bound = onset_at_bound,
      deviance_gap_wide = deviance_gap
    )
  stage_one <- stage_one |>
    dplyr::left_join(
      stage_one_wide,
      by = dplyr::join_by(prefecture, department),
      relationship = "one-to-one"
    )
  stage_one <- audit_join(
    exposed |>
      dplyr::select(prefecture, department, department_display, area_primary),
    stage_one,
    dplyr::join_by(prefecture, department),
    "Stage-one estimates"
  )
  flag_counts <- stage_one |>
    dplyr::summarise(
      attempted = dplyr::n(),
      estimated = sum(is.finite(onset) & is.finite(slope)),
      flat = sum(flat_profile, na.rm = TRUE),
      at_bound = sum(onset_at_bound, na.rm = TRUE),
      at_bound_wide = sum(onset_wide_at_bound, na.rm = TRUE),
      under_30 = sum(under_30),
      warnings = sum(nzchar(glm_warnings)),
      errors = sum(nzchar(glm_errors)),
      nonconverged = sum(!converged)
    )
  onset_table <- stage_one |>
    dplyr::count(onset, onset_wide, name = "departments") |>
    dplyr::arrange(onset, onset_wide)

  # Diagnostic: two-slope change point ------------------------------------
  profile_kink <- function(data, grid = seq(22, 32, 0.5)) {
    multi_year <- dplyr::n_distinct(data$year) > 1L
    # A fit with an aliased coefficient (e.g. every day above k, so the hinge
    # term equals wbgt_max - k) is treated as unavailable at that k.
    fit_quiet <- function(formula) {
      fit <- tryCatch(
        suppressWarnings(stats::glm(
          formula,
          family = stats::poisson(),
          data = data
        )),
        error = function(e) NULL
      )
      if (!is.null(fit) && anyNA(stats::coef(fit))) NULL else fit
    }
    base <- fit_quiet(
      if (multi_year) {
        transports_analysis ~ wbgt_max + factor(year)
      } else {
        transports_analysis ~ wbgt_max
      }
    )
    fits <- lapply(grid, function(k) {
      fit_quiet(
        if (multi_year) {
          transports_analysis ~ I(wbgt_max - k) +
            pmax(wbgt_max - k, 0) +
            factor(year)
        } else {
          transports_analysis ~ I(wbgt_max - k) + pmax(wbgt_max - k, 0)
        }
      )
    })
    deviances <- vapply(
      fits,
      function(x) if (is.null(x)) Inf else x$deviance,
      numeric(1)
    )
    best <- order(deviances, grid)[1]
    valid <- is.finite(deviances[best]) && !is.null(base)
    coefficient <- if (valid) summary(fits[[best]])$coefficients else NULL
    tibble::tibble(
      kink_k = if (valid) grid[best] else NA_real_,
      slope_below = if (valid) coefficient[2, 1] else NA_real_,
      kink = if (valid) coefficient[3, 1] else NA_real_,
      kink_se = if (valid) coefficient[3, 2] else NA_real_,
      improvement = if (valid) base$deviance - deviances[best] else NA_real_,
      kink_at_bound = if (valid) grid[best] %in% range(grid) else NA,
      loglinear_slope = if (!is.null(base)) {
        stats::coef(base)[["wbgt_max"]]
      } else {
        NA_real_
      }
    )
  }
  kink <- panel |>
    dplyr::group_by(prefecture, department) |>
    dplyr::group_modify(~ profile_kink(.x)) |>
    dplyr::ungroup()
  kink <- audit_join(
    exposed |>
      dplyr::select(prefecture, department, department_display, area_primary),
    kink,
    dplyr::join_by(prefecture, department),
    "Kink diagnostic"
  )
  kink_counts <- kink |>
    dplyr::summarise(
      departments = dplyr::n(),
      estimated = sum(is.finite(kink)),
      at_bound = sum(kink_at_bound, na.rm = TRUE),
      improvement_over_3.84 = sum(improvement > 3.84, na.rm = TRUE),
      improvement_over_10 = sum(improvement > 10, na.rm = TRUE),
      positive_kink = sum(kink > 0, na.rm = TRUE),
      median_improvement = stats::median(improvement, na.rm = TRUE),
      median_kink = stats::median(kink, na.rm = TRUE),
      median_loglinear_slope = stats::median(loglinear_slope, na.rm = TRUE)
    )
  kink_location <- kink |>
    dplyr::filter(improvement > 3.84) |>
    dplyr::count(kink_k, name = "departments_over_3.84")

  # Stage two: within- and between-area dispersion ------------------------
  capture_fit <- function(expr) {
    conditions <- character()
    fit <- tryCatch(
      withCallingHandlers(
        expr,
        warning = function(w) {
          conditions <<- c(conditions, conditionMessage(w))
          invokeRestart("muffleWarning")
        },
        message = function(m) {
          conditions <<- c(conditions, conditionMessage(m))
          invokeRestart("muffleMessage")
        }
      ),
      error = function(e) {
        conditions <<- c(conditions, conditionMessage(e))
        NULL
      }
    )
    list(fit = fit, conditions = unique(conditions))
  }
  onset_input <- stage_one |> dplyr::filter(is.finite(onset))
  slope_input <- stage_one |>
    dplyr::filter(is.finite(slope), is.finite(se), se > 0)
  onset_fit <- capture_fit(lme4::lmer(
    onset ~ 1 + (1 | area_primary),
    data = onset_input
  ))
  slope_fit <- capture_fit(lme4::lmer(
    slope ~ 1 + (1 | area_primary),
    weights = 1 / se^2,
    data = slope_input
  ))
  # Unweighted companion so the within-area SD is on the slope's own scale.
  slope_fit_unweighted <- capture_fit(lme4::lmer(
    slope ~ 1 + (1 | area_primary),
    data = slope_input
  ))
  dispersion_row <- function(result, variable, data) {
    fit <- result$fit
    vc <- if (!is.null(fit)) as.data.frame(lme4::VarCorr(fit)) else NULL
    between <- if (!is.null(vc)) {
      vc$sdcor[vc$grp == "area_primary"]
    } else {
      NA_real_
    }
    within <- if (!is.null(vc)) vc$sdcor[vc$grp == "Residual"] else NA_real_
    tibble::tibble(
      variable = variable,
      departments = nrow(data),
      areas = dplyr::n_distinct(data$area_primary),
      between_sd = between,
      within_sd = within,
      between_within_ratio = between / within,
      singular = if (!is.null(fit)) lme4::isSingular(fit) else NA,
      conditions = paste(
        c(result$conditions, if (!is.null(fit)) fit@optinfo$conv$lme4$messages),
        collapse = " | "
      )
    )
  }
  stage_two <- dplyr::bind_rows(
    dispersion_row(onset_fit, "onset", onset_input),
    dispersion_row(slope_fit, "slope (weighted 1/se^2)", slope_input),
    dispersion_row(slope_fit_unweighted, "slope (unweighted)", slope_input)
  )

  # Complementary GLMM ----------------------------------------------------
  k_common <- median(stage_one$onset, na.rm = TRUE)
  glmm_data <- panel |> dplyr::mutate(hinge = pmax(wbgt_max - k_common, 0))
  glmm_start <- proc.time()[["elapsed"]]
  # A child process enforces the wall-clock bound even during native optimization.
  glmm_result <- tryCatch(
    callr::r(
      fit_paper_c_result1_glmm,
      args = list(data = glmm_data),
      timeout = 1800
    ),
    error = function(e) list(error = conditionMessage(e))
  )
  glmm_elapsed <- proc.time()[["elapsed"]] - glmm_start
  if (is.null(glmm_result$error)) {
    glmm_summary <- tibble::tibble(
      k_common = k_common,
      seconds = glmm_elapsed,
      converged = all(glmm_result$optimizer_code == 0) &&
        is.null(glmm_result$convergence_messages),
      singular = glmm_result$singular,
      fixed_slope = glmm_result$fixed_slope,
      fixed_slope_se = glmm_result$fixed_slope_se,
      department_slope_sd = glmm_result$variance$sdcor[
        glmm_result$variance$grp == "department" &
          glmm_result$variance$var1 == "hinge" &
          is.na(glmm_result$variance$var2)
      ],
      area_slope_sd = glmm_result$variance$sdcor[
        glmm_result$variance$grp == "area_primary" &
          glmm_result$variance$var1 == "hinge" &
          is.na(glmm_result$variance$var2)
      ],
      conditions = paste(
        c(glmm_result$conditions, glmm_result$convergence_messages),
        collapse = " | "
      )
    )
  } else {
    glmm_summary <- tibble::tibble(
      k_common = k_common,
      seconds = glmm_elapsed,
      converged = FALSE,
      singular = NA,
      fixed_slope = NA_real_,
      fixed_slope_se = NA_real_,
      department_slope_sd = NA_real_,
      area_slope_sd = NA_real_,
      conditions = glmm_result$error
    )
  }

  # Descriptive regressions of the concordance shares ---------------------
  regression_data <- audit_join(
    shares,
    covariates |>
      dplyr::select(
        prefecture,
        department,
        prop_age_class3,
        log_density,
        coastal,
        log_distance
      ),
    dplyr::join_by(prefecture, department),
    "Regression covariates"
  )
  regression_tables <- list()
  regression_counts <- list()
  for (outcome in c("share_high_given_alert", "share_alert_given_high")) {
    denominator <- if (outcome == "share_high_given_alert") {
      "alert_days"
    } else {
      "high_days"
    }
    trial_data <- regression_data |>
      dplyr::mutate(trials = .data[[denominator]]) |>
      dplyr::filter(trials > 0)
    complete <- with(
      trial_data,
      is.finite(prop_age_class3) &
        is.finite(log_density) &
        !is.na(coastal) &
        is.finite(log_distance)
    )
    used <- trial_data[complete, ]
    for (area_effect in c(TRUE, FALSE)) {
      model_name <- paste(outcome, if (area_effect) "area FE" else "no area FE")
      formula <- stats::as.formula(paste(
        "cbind(joint_days, trials - joint_days) ~ prop_age_class3 + log_density + coastal + log_distance",
        if (area_effect) "+ factor(area_primary)" else ""
      ))
      result <- capture_fit(stats::glm(
        formula,
        family = stats::binomial(),
        data = used
      ))
      covariance <- if (!is.null(result$fit)) {
        capture_fit(sandwich::vcovCL(
          result$fit,
          cluster = used$area_primary
        ))
      } else {
        list(fit = NULL, conditions = character())
      }
      coefficients <- if (!is.null(result$fit)) {
        stats::coef(result$fit)
      } else {
        numeric()
      }
      se <- rep(NA_real_, length(coefficients))
      names(se) <- names(coefficients)
      if (!is.null(covariance$fit)) {
        se[names(diag(covariance$fit))] <- sqrt(diag(covariance$fit))
      }
      regression_tables[[model_name]] <- tibble::tibble(
        model = model_name,
        term = names(coefficients),
        estimate = as.numeric(coefficients),
        clustered_se = as.numeric(se)
      )
      regression_counts[[model_name]] <- tibble::tibble(
        model = model_name,
        departments_used = nrow(used),
        clusters = dplyr::n_distinct(used$area_primary),
        zero_denominator = sum(regression_data[[denominator]] == 0),
        covariate_dropped = sum(!complete),
        covariate_dropped_names = paste(
          trial_data$department_display[!complete],
          collapse = "; "
        ),
        converged = if (!is.null(result$fit)) result$fit$converged else FALSE,
        conditions = paste(
          c(result$conditions, covariance$conditions),
          collapse = " | "
        )
      )
    }
  }
  regression_coefficients <- purrr::list_rbind(regression_tables)
  regression_usage <- purrr::list_rbind(regression_counts)
  # Ratio of |coefficient| to clustered SE, so the text cannot drift from the table.
  covariate_ratios <- regression_coefficients |>
    dplyr::filter(
      term %in%
        c("prop_age_class3", "log_density", "coastalTRUE", "log_distance")
    ) |>
    dplyr::mutate(ratio = abs(estimate) / clustered_se) |>
    dplyr::arrange(dplyr::desc(ratio))
  primary_ratios <- covariate_ratios |>
    dplyr::filter(!grepl("no area FE", model, fixed = TRUE))

  report <- list(
    cells = cells,
    inventory = inventory,
    population_gaps = all_population |>
      dplyr::filter(is.na(population)) |>
      dplyr::select(prefecture, department_display),
    stage_one_counts = flag_counts,
    onset_table = onset_table,
    kink_counts = kink_counts,
    kink_location = kink_location,
    no_exposure = no_exposure |>
      dplyr::select(prefecture, department_display, years, wbgt_missing),
    slope_quantiles = stats::quantile(
      stage_one$slope,
      c(0.05, 0.25, 0.5, 0.75, 0.95),
      na.rm = TRUE
    ),
    stage_two = stage_two,
    glmm = glmm_summary,
    regression_usage = regression_usage,
    regression_coefficients = regression_coefficients,
    zero_denominators = c(
      high_given_alert = sum(shares$alert_days == 0),
      alert_given_high = sum(shares$high_days == 0)
    )
  )
  list(
    assignments = assignments,
    panel = panel,
    no_exposure = no_exposure,
    exposed = exposed,
    cells = cells,
    area_year = area_year,
    year_cells = year_cells,
    shares = shares,
    all_population = all_population,
    covariates = covariates,
    covariate_join_audit = covariate_join_audit,
    inventory = inventory,
    stage_one = stage_one,
    flag_counts = flag_counts,
    onset_table = onset_table,
    kink = kink,
    kink_counts = kink_counts,
    kink_location = kink_location,
    stage_two = stage_two,
    glmm_summary = glmm_summary,
    regression_usage = regression_usage,
    regression_coefficients = regression_coefficients,
    covariate_ratios = covariate_ratios,
    primary_ratios = primary_ratios,
    join_audit = join_audit,
    report = report
  )
}

# Alert-area-year concordance for the manuscript (paper C, step 5g (viii)),
# moved from src/paper_c_area_year_shares.R (Issue #50 (a)). Reproduces the
# area_year cells of result 1, including combinations without observations,
# and checks the totals against `result1_cells` (the report's cells). Shares
# are weighted by department-days; a zero denominator is undefined, never a
# zero share. Only area-level aggregates are returned; no department name.
build_paper_c_area_year_shares <- function(
  df_paper_c_panel,
  df_fd_alert_area,
  result1_cells
) {
  panel <- df_paper_c_panel |>
    dplyr::filter(!is.na(local_high)) |>
    dplyr::select(
      area_primary,
      year,
      alert_status,
      local_high,
      transports_analysis
    )
  assignments <- df_fd_alert_area |>
    dplyr::distinct(area_primary)
  reference <- result1_cells
  statuses <- c("alert", "active, no alert", "structural no-alert")
  # English names of the 30 alert forecast areas in the sample (Hepburn
  # romanization; the eight Hokkaido areas are sub-prefectural and marked).
  area_english <- c(
    "三重" = "Mie",
    "佐賀" = "Saga",
    "兵庫" = "Hyogo",
    "和歌山" = "Wakayama",
    "埼玉" = "Saitama",
    "山口" = "Yamaguchi",
    "山形" = "Yamagata",
    "山梨" = "Yamanashi",
    "岩手" = "Iwate",
    "島根" = "Shimane",
    "徳島" = "Tokushima",
    "愛媛" = "Ehime",
    "愛知" = "Aichi",
    "新潟" = "Niigata",
    "栃木" = "Tochigi",
    "神奈川" = "Kanagawa",
    "福井" = "Fukui",
    "福島" = "Fukushima",
    "群馬" = "Gunma",
    "長野" = "Nagano",
    "静岡" = "Shizuoka",
    "香川" = "Kagawa",
    "上川・留萌" = "Kamikawa-Rumoi (Hokkaido)",
    "十勝" = "Tokachi (Hokkaido)",
    "宗谷" = "Soya (Hokkaido)",
    "渡島・檜山" = "Oshima-Hiyama (Hokkaido)",
    "石狩・空知・後志" = "Ishikari-Sorachi-Shiribeshi (Hokkaido)",
    "網走・北見・紋別" = "Abashiri-Kitami-Monbetsu (Hokkaido)",
    "胆振・日高" = "Iburi-Hidaka (Hokkaido)",
    "釧路・根室" = "Kushiro-Nemuro (Hokkaido)",
    # Issue #59
    "富山" = "Toyama",
    "沖縄本島" = "Okinawa Main Island (Okinawa)",
    "宮古島" = "Miyakojima (Okinawa)",
    "八重山" = "Yaeyama (Okinawa)"
  )
  safe_share <- function(n, d) {
    n / dplyr::if_else(d > 0, as.numeric(d), NA_real_)
  }
  cell_summary <- function(data, groups = character()) {
    data |>
      dplyr::summarise(
        fd_days = dplyr::n(),
        total_transports = sum(transports_analysis),
        .by = dplyr::all_of(c(groups, "alert_status", "local_high"))
      )
  }
  stopifnot(
    nrow(panel) == 215613L,
    !anyNA(panel),
    !anyNA(assignments),
    all(panel$alert_status %in% statuses),
    setequal(panel$area_primary, assignments$area_primary),
    setequal(assignments$area_primary, names(area_english)),
    !anyDuplicated(area_english),
    setequal(panel$year, 2018:2022)
  )
  area_year_cells <- cell_summary(panel, c("area_primary", "year")) |>
    tidyr::complete(
      area_primary = unique(assignments$area_primary),
      year = 2018:2022,
      alert_status = statuses,
      local_high = c(FALSE, TRUE),
      fill = list(fd_days = 0L, total_transports = 0)
    )
  area_year <- area_year_cells |>
    dplyr::summarise(
      n_all = sum(fd_days),
      alert_days = sum(fd_days[alert_status == "alert"]),
      alert_high_share = safe_share(
        sum(fd_days[alert_status == "alert" & local_high]),
        alert_days
      ),
      active_no_alert_days = sum(fd_days[alert_status == "active, no alert"]),
      active_no_alert_high_share = safe_share(
        sum(fd_days[alert_status == "active, no alert" & local_high]),
        active_no_alert_days
      ),
      structural_days = sum(fd_days[alert_status == "structural no-alert"]),
      .by = c(area_primary, year)
    ) |>
    dplyr::rename(area = area_primary, fd_days = n_all) |>
    dplyr::mutate(area_en = unname(area_english[area]), .after = area) |>
    dplyr::arrange(area_en, year)
  meta <- tibble::tibble(
    rows = nrow(area_year),
    areas = dplyr::n_distinct(area_year$area),
    years = dplyr::n_distinct(area_year$year),
    fd_days_total = sum(area_year$fd_days),
    alert_days_total = sum(area_year$alert_days),
    transports_total = sum(area_year_cells$total_transports)
  )
  stopifnot(
    meta$rows == 170L,
    meta$areas == 34L,
    meta$years == 5L,
    meta$fd_days_total == 215613L,
    meta$fd_days_total == sum(reference$fd_days),
    meta$alert_days_total ==
      sum(reference$fd_days[reference$alert_status == "alert"]),
    sum(area_year$active_no_alert_days) ==
      sum(reference$fd_days[reference$alert_status == "active, no alert"]),
    meta$transports_total == sum(reference$total_transports),
    !anyNA(area_year[c(
      "fd_days",
      "alert_days",
      "active_no_alert_days",
      "structural_days"
    )]),
    all(
      area_year$fd_days ==
        area_year$alert_days +
          area_year$active_no_alert_days +
          area_year$structural_days
    ),
    all(is.na(area_year$alert_high_share) == (area_year$alert_days == 0)),
    all(
      is.na(area_year$active_no_alert_high_share) ==
        (area_year$active_no_alert_days == 0)
    )
  )
  list(area_year = area_year, meta = meta)
}
