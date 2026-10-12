# NOTE: readable analysis code from the private project. It needs the fire-department
# transport records, which are not distributed, and is not sourced by the public
# pipeline (_targets.R). Kept so that the published estimates can be audited.

# Power simulation for paper C under the pre-specified model (paper_strategy 7.2).
#
# The estimating model is a Poisson regression of fire-department x day
# transports on the alert indicators, with department, date, area-maximum WBGT
# bins at t and t-1, and local WBGT bins as fixed effects, clustered on the
# alert area. The data-generating process is anchored on the null fit
# (the model without `alert_high`): a parametric negative binomial draw with
# the effect applied to the `alert_high` rows, plus an analytic noncentral-t
# curve from the real-data clustered SE.
#
# Every function is pure. The driver in `src/paper_c_power.R` reads the store
# and calls them in order.

# Prefectures whose raw target carries a severity column (中等症 or 重症).
# Same rule as notes/paper_c_alert_join.qmd, lifted into a function.
paper_c_severity_prefectures <- function(
  store,
  sample_prefectures,
  expected_count = 14L
) {
  codes <- zipangu::jpnprefs |>
    dplyr::filter(prefecture_kanji %in% sample_prefectures) |>
    dplyr::pull(jis_code)
  pref_data <- purrr::map(codes, function(code) {
    targets::tar_read_raw(paste0("df_hs_pref", code), store = store)
  }) |>
    rlang::set_names(codes)
  paper_c_severity_prefectures_from_data(
    pref_data,
    sample_prefectures,
    expected_count
  )
}

# The same rule on the raw prefecture targets passed in, named by JIS code, so
# that a target can declare them as dependencies instead of reading the store.
# `expected_count` defaults to the 14 of the pre-specified 23-prefecture
# sample, which the scripts in src/ read; the targets pass the count of the
# current sample (Issue #59).
paper_c_severity_prefectures_from_data <- function(
  pref_data,
  sample_prefectures,
  expected_count = 14L
) {
  candidates <- zipangu::jpnprefs |>
    dplyr::filter(prefecture_kanji %in% sample_prefectures)
  missing <- setdiff(candidates$jis_code, names(pref_data))
  if (length(missing) > 0L) {
    rlang::abort(glue::glue(
      "No raw target was passed for sample prefecture codes: ",
      "{paste(missing, collapse = ', ')}."
    ))
  }
  has_severity <- purrr::map_lgl(candidates$jis_code, function(code) {
    any(c("中等症", "重症") %in% names(pref_data[[code]]))
  })
  severity <- candidates$prefecture_kanji[has_severity]
  if (length(severity) != expected_count) {
    rlang::abort(glue::glue(
      "Expected {expected_count} severity prefectures, got {length(severity)}: ",
      "{paste(severity, collapse = ', ')}."
    ))
  }
  severity
}

# Assign WBGT stations to the 58 alert-area polygons. Same discipline as
# assign_fd_alert_area(): make both layers valid and project to EPSG:6933
# before the spatial predicate. Stations in no polygon (coastal or island
# points) take the nearest polygon, and the count of these fallbacks is
# returned with the audit.
assign_wbgt_station_alert_area <- function(
  wbgt_station_geometry,
  alert_geometry
) {
  if (
    nrow(alert_geometry) != 58L || dplyr::n_distinct(alert_geometry$area) != 58L
  ) {
    rlang::abort("Expected 58 unique heat-alert areas.")
  }
  if (
    is.na(sf::st_crs(wbgt_station_geometry)) ||
      is.na(sf::st_crs(alert_geometry))
  ) {
    rlang::abort("Station and alert-area geometries need a CRS.")
  }

  stations <- wbgt_station_geometry |>
    dplyr::mutate(station_no = as.character(station_no)) |>
    sf::st_make_valid() |>
    sf::st_transform(6933)
  areas <- alert_geometry |>
    dplyr::select(area) |>
    sf::st_make_valid() |>
    sf::st_transform(6933)

  within <- sf::st_within(stations, areas)
  if (any(lengths(within) > 1L)) {
    rlang::abort("A station falls within more than one alert-area polygon.")
  }
  area_index <- vapply(
    within,
    function(x) if (length(x) == 1L) x[[1L]] else NA_integer_,
    integer(1)
  )
  fallback <- is.na(area_index)
  if (any(fallback)) {
    area_index[fallback] <- sf::st_nearest_feature(
      stations[fallback, ],
      areas
    )
  }

  assignment <- stations |>
    sf::st_drop_geometry() |>
    dplyr::transmute(
      station_no,
      station_name,
      station_label = area,
      alert_area = areas$area[area_index],
      nearest_fallback = fallback
    )

  # Audit 1: a station label that is itself an alert-area name must match.
  labelled <- assignment |>
    dplyr::filter(station_label %in% alert_geometry$area)
  label_mismatch <- labelled |>
    dplyr::filter(station_label != alert_area)
  if (nrow(label_mismatch) > 0L) {
    rlang::abort(paste0(
      "Station label disagrees with the spatial alert area for ",
      nrow(label_mismatch),
      " station(s): ",
      paste(
        label_mismatch$station_no,
        label_mismatch$station_label,
        label_mismatch$alert_area,
        sep = "/",
        collapse = ", "
      )
    ))
  }

  # Audit 2: Hokkaido labels finer than the alert areas must be a
  # `・`-separated component of the spatial area name.
  hokkaido_labels <- c(
    "上川",
    "留萌",
    "石狩",
    "空知",
    "後志",
    "根室",
    "釧路",
    "胆振",
    "日高",
    "渡島",
    "檜山"
  )
  hokkaido <- assignment |>
    dplyr::filter(station_label %in% hokkaido_labels) |>
    dplyr::mutate(
      component = purrr::map2_lgl(
        station_label,
        alert_area,
        function(label, area) label %in% strsplit(area, "・")[[1L]]
      )
    )
  if (!all(hokkaido$component)) {
    bad <- hokkaido |>
      dplyr::filter(!component)
    rlang::abort(paste0(
      "Hokkaido station label is not a component of its spatial area for ",
      nrow(bad),
      " station(s): ",
      paste(bad$station_label, bad$alert_area, sep = "/", collapse = ", ")
    ))
  }

  # Audit 3: labels coarser than the alert areas (鹿児島, 沖縄) are only
  # counted; the split is reported, not asserted.
  coarse_split <- assignment |>
    dplyr::filter(station_label %in% c("鹿児島", "沖縄")) |>
    dplyr::count(station_label, alert_area, name = "stations")

  unaudited <- assignment |>
    dplyr::filter(
      !station_label %in% alert_geometry$area,
      !station_label %in% hokkaido_labels,
      !station_label %in% c("鹿児島", "沖縄")
    )
  if (nrow(unaudited) > 0L) {
    rlang::abort(paste0(
      "Station labels outside the three audit rules: ",
      paste(unique(unaudited$station_label), collapse = ", ")
    ))
  }

  list(
    assignment = assignment,
    audit = list(
      stations = nrow(assignment),
      within_polygon = sum(!fallback),
      nearest_fallback = sum(fallback),
      fallback_stations = assignment |>
        dplyr::filter(nearest_fallback) |>
        dplyr::select(station_no, station_name, station_label, alert_area),
      label_matches = nrow(labelled),
      hokkaido_component_matches = nrow(hokkaido),
      coarse_split = coarse_split
    )
  )
}

# Area-maximum WBGT per alert area x date over the full station series, so the
# t-1 lag exists on June 1. The lag is a self-join on date - 1 rather than
# dplyr::lag(): the series has between-season gaps (October to April) and a
# positional lag would carry October into the next April.
build_area_wbgt_max <- function(wbgt_daily, station_assignment) {
  station_daily <- wbgt_daily |>
    dplyr::mutate(station_no = as.character(station_no)) |>
    dplyr::filter(is.finite(wbgt_max)) |>
    dplyr::inner_join(
      station_assignment |>
        dplyr::select(station_no, alert_area),
      by = dplyr::join_by(station_no),
      relationship = "many-to-one"
    )
  area_daily <- station_daily |>
    dplyr::summarise(
      area_wbgt_max = max(wbgt_max),
      stations_reporting = dplyr::n(),
      .by = c(alert_area, date)
    )
  area_daily |>
    dplyr::left_join(
      area_daily |>
        dplyr::transmute(
          alert_area,
          date = date + 1L,
          area_wbgt_max_lag1 = area_wbgt_max
        ),
      by = dplyr::join_by(alert_area, date),
      relationship = "one-to-one"
    ) |>
    pointblank::rows_distinct(columns = c(alert_area, date))
}

# Attach the area-maximum WBGT and the integer bin codes to the analysis rows.
# Bins are integer codes (round(x * 10) at 0.1 degC, round(x) at 1 degC), never
# factors, so the session locale cannot reorder them.
build_paper_c_power_design <- function(
  panel,
  area_wbgt_max,
  taxable_income = NULL,
  income_fiscal_year = 2020L,
  expected_departments = 391L
) {
  required <- c(
    "prefecture",
    "date",
    "department",
    "year",
    "area_primary",
    "near_wbgt_station",
    "wbgt_max",
    "alert_day",
    "local_high",
    "transports_analysis"
  )
  missing <- setdiff(required, names(panel))
  if (length(missing) > 0L) {
    rlang::abort(paste(
      "Panel lacks columns:",
      paste(missing, collapse = ", ")
    ))
  }

  design <- panel |>
    dplyr::filter(!is.na(local_high)) |>
    dplyr::select(dplyr::all_of(required)) |>
    dplyr::left_join(
      area_wbgt_max |>
        dplyr::select(
          alert_area,
          date,
          area_wbgt_max,
          area_wbgt_max_lag1,
          stations_reporting
        ),
      by = dplyr::join_by(area_primary == alert_area, date),
      relationship = "many-to-one"
    )

  if (!is.null(taxable_income)) {
    panel_departments <- panel |>
      dplyr::distinct(prefecture, department)
    income <- taxable_income |>
      dplyr::filter(fiscal_year == income_fiscal_year) |>
      dplyr::select(
        prefecture,
        department,
        income_per_taxpayer,
        income_per_resident
      )
    duplicate_income <- income |>
      dplyr::count(prefecture, department) |>
      dplyr::filter(n != 1L)
    missing_income <- dplyr::anti_join(
      panel_departments,
      income,
      by = dplyr::join_by(prefecture, department)
    )
    extra_income <- dplyr::anti_join(
      income,
      panel_departments,
      by = dplyr::join_by(prefecture, department)
    )
    if (
      nrow(income) != expected_departments ||
        nrow(panel_departments) != expected_departments ||
        nrow(duplicate_income) > 0L ||
        nrow(missing_income) > 0L ||
        nrow(extra_income) > 0L
    ) {
      rlang::abort(glue::glue(
        "Expected an exact {expected_departments}-department income join ",
        "for fiscal year ",
        "{income_fiscal_year}; got {nrow(income)} income rows, ",
        "{nrow(duplicate_income)} duplicate keys, {nrow(missing_income)} ",
        "missing keys, and {nrow(extra_income)} extra keys."
      ))
    }
    if (
      any(!is.finite(income$income_per_taxpayer)) ||
        any(!is.finite(income$income_per_resident))
    ) {
      rlang::abort("The income modifiers must be finite for all departments.")
    }
    taxpayer_mean <- mean(income$income_per_taxpayer)
    taxpayer_sd <- stats::sd(income$income_per_taxpayer)
    resident_mean <- mean(income$income_per_resident)
    resident_sd <- stats::sd(income$income_per_resident)
    if (
      !is.finite(taxpayer_sd) ||
        taxpayer_sd <= 0 ||
        !is.finite(resident_sd) ||
        resident_sd <= 0
    ) {
      rlang::abort(
        "The income modifiers need positive finite standard deviations."
      )
    }
    income <- income |>
      dplyr::mutate(
        income_z = (income_per_taxpayer - taxpayer_mean) / taxpayer_sd,
        income_resident_z = (income_per_resident - resident_mean) / resident_sd
      )
    rows_before_income <- nrow(design)
    design <- design |>
      dplyr::left_join(
        income,
        by = dplyr::join_by(prefecture, department),
        relationship = "many-to-one"
      )
    if (
      nrow(design) != rows_before_income ||
        anyNA(design$income_z) ||
        anyNA(design$income_resident_z)
    ) {
      rlang::abort(glue::glue(
        "Income join changed rows from {rows_before_income} to {nrow(design)} ",
        "or left {sum(is.na(design$income_z))} taxpayer / ",
        "{sum(is.na(design$income_resident_z))} resident z values missing."
      ))
    }
  }

  if (anyNA(design$area_wbgt_max) || anyNA(design$area_wbgt_max_lag1)) {
    rlang::abort(glue::glue(
      "Area-maximum WBGT missing on {sum(is.na(design$area_wbgt_max))} rows ",
      "and its lag on {sum(is.na(design$area_wbgt_max_lag1))} rows."
    ))
  }
  alert_rows <- design |>
    dplyr::filter(alert_day)
  if (anyNA(alert_rows$area_wbgt_max) || anyNA(alert_rows$area_wbgt_max_lag1)) {
    rlang::abort("An alert row lacks the area-maximum WBGT or its lag.")
  }

  design <- design |>
    dplyr::mutate(
      area_bin_t = as.integer(round(area_wbgt_max * 10)),
      area_bin_lag = as.integer(round(area_wbgt_max_lag1 * 10)),
      local_bin = as.integer(round(wbgt_max * 10)),
      area_bin_t_1c = as.integer(round(area_wbgt_max)),
      area_bin_lag_1c = as.integer(round(area_wbgt_max_lag1)),
      local_bin_1c = as.integer(round(wbgt_max)),
      alert_day = as.numeric(alert_day)
    )
  alert_terms <- c("alert_trial", "alert_high")
  if (!is.null(taxable_income)) {
    alert_terms <- c(alert_terms, "alert_income", "alert_income_resident")
  }
  rebuild_alert_terms(design, alert_terms)
}

# Guard the two samples against the counts fixed in the brief.
check_paper_c_power_sample <- function(
  design,
  expected_rows,
  expected_cells,
  expected_areas,
  expected_alerted_areas,
  label
) {
  cells <- design |>
    dplyr::count(alert_day, local_high) |>
    dplyr::arrange(alert_day, local_high) |>
    dplyr::pull(n)
  areas <- dplyr::n_distinct(design$area_primary)
  alerted <- design |>
    dplyr::filter(alert_day == 1) |>
    dplyr::pull(area_primary) |>
    dplyr::n_distinct()
  observed <- c(nrow(design), cells, areas, alerted)
  expected <- c(
    expected_rows,
    expected_cells,
    expected_areas,
    expected_alerted_areas
  )
  if (!identical(as.integer(observed), as.integer(expected))) {
    rlang::abort(glue::glue(
      "{label} sample guard failed. Expected rows/cells/areas/alerted = ",
      "{paste(expected, collapse = '/')}, observed {paste(observed, collapse = '/')}."
    ))
  }
  invisible(design)
}

# Model formula from explicit regressors at a given bin resolution.
paper_c_power_formula_terms <- function(
  terms,
  bin_resolution = c("0.1", "1")
) {
  bin_resolution <- match.arg(bin_resolution)
  suffix <- if (bin_resolution == "1") "_1c" else ""
  rhs <- if (length(terms) == 0L) "1" else paste(terms, collapse = " + ")
  stats::as.formula(paste0(
    "transports_analysis ~ ",
    rhs,
    " | department + date + area_bin_t",
    suffix,
    " + area_bin_lag",
    suffix,
    " + local_bin",
    suffix
  ))
}

# Model formulas for the local-high power analysis.
paper_c_power_formula <- function(bin_resolution = c("0.1", "1"), full = TRUE) {
  terms <- c("alert_day", "alert_trial", if (full) "alert_high")
  paper_c_power_formula_terms(terms, match.arg(bin_resolution))
}

paper_c_power_ssc <- function() {
  fixest::ssc(adj = TRUE, cluster.adj = TRUE)
}

# Fit the Poisson model. `outcome` replaces `transports_analysis` when given.
fit_paper_c_power <- function(design, formula, outcome = NULL) {
  if (!is.null(outcome)) {
    design$transports_analysis <- outcome
  }
  fixest::fepois(
    formula,
    data = design,
    cluster = ~area_primary,
    ssc = paper_c_power_ssc(),
    nthreads = 1L,
    notes = FALSE
  )
}

# Row indices of `design` the fit kept (fixest drops singleton / all-zero
# fixed-effect levels).
fit_kept_rows <- function(fit, n_rows) {
  removed <- fit$obs_selection$obsRemoved
  if (is.null(removed)) {
    return(seq_len(n_rows))
  }
  seq_len(n_rows)[removed]
}

# Cluster-robust inference on one coefficient: t = est / se, p from the t
# distribution with G - 1 degrees of freedom, G = clusters in the fitted rows.
coefficient_test <- function(fit, term, n_clusters) {
  estimate <- stats::coef(fit)[[term]]
  se <- fixest::se(
    fit,
    cluster = ~area_primary,
    ssc = paper_c_power_ssc()
  )[[term]]
  statistic <- estimate / se
  tibble::tibble(
    term = term,
    estimate = estimate,
    se = se,
    statistic = statistic,
    clusters = n_clusters,
    p_value = 2 * stats::pt(-abs(statistic), df = n_clusters - 1L)
  )
}

# Coefficient table with the G - 1 t-based p-value next to fixest's z-based one.
coefficient_table <- function(fit, n_clusters) {
  z_table <- fixest::coeftable(
    fit,
    cluster = ~area_primary,
    ssc = paper_c_power_ssc()
  )
  purrr::map(
    rownames(z_table),
    coefficient_test,
    fit = fit,
    n_clusters = n_clusters
  ) |>
    purrr::list_rbind() |>
    dplyr::mutate(p_value_z = z_table[, "Pr(>|z|)"])
}

# Removed fixed-effect levels and observations, from the fit object.
fit_removals <- function(fit) {
  list(
    nobs = fit$nobs,
    nobs_origin = fit$nobs_origin,
    obs_removed = fit$nobs_origin - fit$nobs,
    fixef_sizes = fit$fixef_sizes,
    fixef_removed = purrr::map_int(fit$fixef_removed, length),
    collinear = fit$collin.var
  )
}

# Moment estimate of the negative-binomial size from the Poisson null fit:
# 1 / theta = mean((y - mu)^2 - mu) / mean(mu^2).
theta_moment <- function(y, mu) {
  inverse <- mean((y - mu)^2 - mu) / mean(mu^2)
  if (!is.finite(inverse) || inverse <= 0) {
    rlang::abort(glue::glue(
      "Moment estimate of 1/theta is {inverse}; the null fit shows no overdispersion."
    ))
  }
  1 / inverse
}

# Null fit, real-data full fit, fitted means, and theta for one sample.
fit_paper_c_power_anchors <- function(
  design,
  bin_resolution = c("0.1", "1"),
  theta_method = c("moment", "fenegbin"),
  label,
  effect_term = "alert_high",
  terms = c("alert_day", "alert_trial", effect_term),
  real_inference = NULL
) {
  bin_resolution <- match.arg(bin_resolution)
  theta_method <- match.arg(theta_method)
  if (!effect_term %in% terms || !effect_term %in% names(design)) {
    rlang::abort(glue::glue(
      "{label}: `{effect_term}` must occur in both `terms` and the design."
    ))
  }
  null_terms <- setdiff(terms, effect_term)
  null_formula <- paper_c_power_formula_terms(null_terms, bin_resolution)
  full_formula <- paper_c_power_formula_terms(terms, bin_resolution)

  timing <- system.time(
    null_fit <- fit_paper_c_power(design, null_formula)
  )[["elapsed"]]
  timing_full <- system.time(
    full_fit <- fit_paper_c_power(design, full_formula)
  )[["elapsed"]]

  for (fit in list(null_fit, full_fit)) {
    if (anyNA(stats::coef(fit)) || !is.null(fit$collin.var)) {
      rlang::abort(glue::glue(
        "{label}: fixest dropped a regressor (NA coefficient or collinearity)."
      ))
    }
  }
  if (
    length(stats::coef(null_fit)) != length(null_terms) ||
      length(stats::coef(full_fit)) != length(terms)
  ) {
    rlang::abort(glue::glue("{label}: unexpected number of coefficients."))
  }

  kept <- fit_kept_rows(null_fit, nrow(design))
  kept_full <- fit_kept_rows(full_fit, nrow(design))
  if (!identical(kept, kept_full)) {
    rlang::abort(glue::glue(
      "{label}: the null and full fits keep different rows."
    ))
  }
  kept_design <- design[kept, ]
  mu <- unname(stats::fitted(null_fit))
  if (length(mu) != nrow(kept_design) || anyNA(mu)) {
    rlang::abort(glue::glue(
      "{label}: fitted means do not align with kept rows."
    ))
  }
  if (anyNA(kept_design[[effect_term]])) {
    rlang::abort(glue::glue(
      "{label}: `{effect_term}` is missing on a fitted row."
    ))
  }
  dropped_effect_rows <- sum(design[[effect_term]][-kept] != 0)
  n_clusters <- dplyr::n_distinct(kept_design$area_primary)
  real_fit_modifier <- coefficient_test(full_fit, effect_term, n_clusters)
  if (!is.null(real_inference)) {
    if (
      nrow(real_inference) != 1L ||
        abs(real_inference$estimate - real_fit_modifier$estimate) > 1e-10 ||
        abs(real_inference$se - real_fit_modifier$se) > 1e-6 ||
        real_inference$clusters != n_clusters
    ) {
      rlang::abort(glue::glue(
        "{label}: the refit of `{effect_term}` does not reproduce the ",
        "reviewed modifier report."
      ))
    }
    real_modifier <- real_inference
  } else {
    real_modifier <- real_fit_modifier
  }

  theta_estimate <- theta_moment(kept_design$transports_analysis, mu)
  theta <- list(
    method = "moment",
    value = theta_estimate,
    fenegbin = NULL
  )
  if (theta_method == "fenegbin") {
    timing_nb <- system.time(
      nb_fit <- fixest::fenegbin(
        null_formula,
        data = design,
        cluster = ~area_primary,
        nthreads = 1L,
        notes = FALSE
      )
    )[["elapsed"]]
    theta$method <- "fenegbin"
    theta$moment <- theta_estimate
    theta$value <- nb_fit$theta
    theta$fenegbin_seconds <- timing_nb
  }

  anchor <- list(
    label = label,
    effect_term = effect_term,
    bin_resolution = bin_resolution,
    null_formula = null_formula,
    full_formula = full_formula,
    kept_design = kept_design,
    mu = mu,
    n_clusters = n_clusters,
    theta = theta,
    null_coefficients = coefficient_table(null_fit, n_clusters),
    full_coefficients = coefficient_table(full_fit, n_clusters),
    real_modifier = real_modifier,
    real_modifier_refit = real_fit_modifier,
    null_removals = fit_removals(null_fit),
    full_removals = fit_removals(full_fit),
    dropped_effect_rows = dropped_effect_rows,
    seconds_null_fit = timing,
    seconds_full_fit = timing_full
  )
  if (effect_term == "alert_high") {
    anchor$real_alert_high <- real_modifier
    anchor$dropped_alert_high_rows <- dropped_effect_rows
  }
  anchor
}

# Draw one simulated outcome from the parametric negative binomial DGP:
# y ~ NegBin(mean = mu_hat * exp(beta * modifier), size = theta_hat).
simulate_outcome <- function(anchor, beta) {
  mu <- anchor$mu
  effect_term <- if (is.null(anchor$effect_term)) {
    "alert_high"
  } else {
    anchor$effect_term
  }
  modifier <- anchor$kept_design[[effect_term]]
  stats::rnbinom(
    length(mu),
    mu = mu * exp(beta * modifier),
    size = anchor$theta$value
  )
}

# One replication: simulate, refit the full model, and test the modifier.
run_paper_c_power_replication <- function(anchor, beta, seed) {
  assign(".Random.seed", seed, envir = globalenv())
  outcome <- simulate_outcome(anchor, beta)
  tryCatch(
    {
      fit <- fit_paper_c_power(
        anchor$kept_design,
        anchor$full_formula,
        outcome = outcome
      )
      kept <- fit_kept_rows(fit, nrow(anchor$kept_design))
      n_clusters <- dplyr::n_distinct(anchor$kept_design$area_primary[kept])
      effect_term <- if (is.null(anchor$effect_term)) {
        "alert_high"
      } else {
        anchor$effect_term
      }
      test <- coefficient_test(fit, effect_term, n_clusters)
      tibble::tibble(
        beta_hat = test$estimate,
        se = test$se,
        p_value = test$p_value,
        clusters = n_clusters,
        nobs = fit$nobs,
        failed = is.na(test$estimate) || is.na(test$p_value)
      )
    },
    error = function(e) {
      tibble::tibble(
        beta_hat = NA_real_,
        se = NA_real_,
        p_value = NA_real_,
        clusters = NA_integer_,
        nobs = NA_integer_,
        failed = TRUE
      )
    }
  )
}

# Task grid with one L'Ecuyer-CMRG stream per task, so the result depends on
# the seed only and not on the number of workers, the chunking, or a resume.
build_paper_c_power_tasks <- function(
  samples,
  relative_effects,
  replications,
  seed
) {
  grid <- tidyr::expand_grid(
    sample = samples,
    dgp = "nb",
    relative_effect = relative_effects,
    replication = seq_len(replications)
  )
  # The caller's RNG state and kind are restored afterwards, and the global
  # seed is read with get() so that targets does not record `.Random.seed` as
  # a dependency of every target that calls this (Issue #50 self-review).
  grid$stream <- withr::with_preserve_seed({
    set.seed(seed, kind = "L'Ecuyer-CMRG")
    current <- get(".Random.seed", envir = globalenv())
    streams <- vector("list", nrow(grid))
    for (i in seq_len(nrow(grid))) {
      current <- parallel::nextRNGStream(current)
      streams[[i]] <- current
    }
    streams
  })
  grid$chunk <- dplyr::consecutive_id(grid$sample, grid$relative_effect)
  grid
}

# Character key of an L'Ecuyer-CMRG stream, stored with every result row so a
# resumed run can prove that each completed task used the stream the current
# grid assigns to it.
stream_key <- function(stream) {
  paste(stream, collapse = "_")
}

# Write an RDS through a temporary file in the same directory and rename it
# into place, so a reader never sees a half-written file and a failure while
# writing leaves the previous version intact.
save_rds_atomic <- function(object, path) {
  temporary <- tempfile(
    pattern = paste0(basename(path), "."),
    tmpdir = dirname(path),
    fileext = ".tmp"
  )
  saveRDS(object, temporary)
  if (!file.rename(temporary, path)) {
    unlink(temporary)
    rlang::abort(glue::glue("Could not move {temporary} to {path}."))
  }
  invisible(path)
}

# The settings a partial file must share with the run that resumes it.
paper_c_power_partial_settings <- function(settings, theta = NULL) {
  list(
    seed = settings$seed,
    replications = settings$replications,
    bin = settings$bin,
    theta_method = settings$theta,
    theta_value = theta,
    relative_effects = settings$relative_effects,
    effect_term = if (is.null(settings$effect_term)) {
      "alert_high"
    } else {
      settings$effect_term
    },
    fixest_version = settings$fixest_version,
    # Content fingerprints, set by the driver once the designs and anchors
    # exist (NULL for drivers that do not compute them).
    design_hash = settings$design_hash,
    anchor_hash = settings$anchor_hash
  )
}

# Abort unless the partial file belongs to exactly this simulation: same
# settings, and for every completed key the stored stream equals the one the
# current task grid assigns. `theta` is checked only when supplied (it is
# known after the anchor fits).
validate_paper_c_power_partial <- function(
  partial,
  tasks,
  settings,
  partial_path,
  theta = NULL
) {
  if (
    !is.list(partial) ||
      !all(c("settings", "results") %in% names(partial))
  ) {
    rlang::abort(glue::glue(
      "{partial_path} is not a partial file (needs `settings` and `results`)."
    ))
  }
  expected <- paper_c_power_partial_settings(settings, theta)
  stored <- partial$settings
  if (is.null(theta)) {
    # Before the anchors are fitted only the grid can be checked; theta and
    # the anchor fingerprint are compared on the second call.
    expected$theta_value <- NULL
    stored$theta_value <- NULL
    expected$anchor_hash <- NULL
    stored$anchor_hash <- NULL
  }
  for (name in names(expected)) {
    if (!identical(stored[[name]], expected[[name]])) {
      rlang::abort(glue::glue(
        "--resume refused: `{name}` in {partial_path} is ",
        "{paste(format(stored[[name]]), collapse = ', ')} but the current run ",
        "has {paste(format(expected[[name]]), collapse = ', ')}. ",
        "Delete the partial file to start over."
      ))
    }
  }
  keys <- c("sample", "dgp", "relative_effect", "replication")
  results <- partial$results
  if (!all(c(keys, "stream_key") %in% names(results))) {
    rlang::abort(glue::glue(
      "{partial_path} results lack the task keys or the stream key."
    ))
  }
  current <- tasks |>
    dplyr::mutate(current_stream_key = purrr::map_chr(stream, stream_key)) |>
    dplyr::select(dplyr::all_of(keys), current_stream_key)
  matched <- results |>
    dplyr::select(dplyr::all_of(keys), stream_key) |>
    dplyr::left_join(current, by = keys)
  unknown <- matched |>
    dplyr::filter(is.na(current_stream_key))
  if (nrow(unknown) > 0L) {
    rlang::abort(glue::glue(
      "--resume refused: {nrow(unknown)} completed task(s) in {partial_path} ",
      "are not in the current task grid."
    ))
  }
  mismatched <- matched |>
    dplyr::filter(stream_key != current_stream_key)
  if (nrow(mismatched) > 0L) {
    rlang::abort(glue::glue(
      "--resume refused: {nrow(mismatched)} completed task(s) in ",
      "{partial_path} used a different random stream than the current grid ",
      "assigns (first: {mismatched$sample[[1L]]}, effect ",
      "{mismatched$relative_effect[[1L]]}, replication ",
      "{mismatched$replication[[1L]]})."
    ))
  }
  invisible(partial)
}

# One chunk of the task grid with forked workers: the task keys, the key of
# each task's stream, and the fit. Every task carries its own L'Ecuyer-CMRG
# stream, so the rows depend on the tasks alone, not on the worker count or
# the chunking.
run_paper_c_power_chunk <- function(
  anchors,
  chunk,
  workers,
  failure_note = ""
) {
  keys <- c("sample", "dgp", "relative_effect", "replication")
  worker <- function(i) {
    fixest::setFixest_nthreads(1L)
    task <- chunk[i, ]
    run_paper_c_power_replication(
      anchor = anchors[[task$sample]],
      beta = log1p(task$relative_effect),
      seed = task$stream[[1L]]
    )
  }
  chunk_results <- parallel::mclapply(
    seq_len(nrow(chunk)),
    worker,
    mc.cores = workers,
    mc.preschedule = FALSE,
    mc.set.seed = FALSE
  )
  broken <- purrr::map_lgl(chunk_results, function(x) {
    is.null(x) || inherits(x, "try-error")
  })
  if (any(broken)) {
    first_error <- purrr::detect(chunk_results, inherits, "try-error")
    rlang::abort(paste0(
      failure_note,
      if (nzchar(failure_note)) ": " else "",
      sum(broken),
      " NULL or error result(s)",
      if (!is.null(first_error)) {
        paste0(": ", conditionMessage(attr(first_error, "condition")))
      } else {
        " (NULL: a forked worker died)"
      }
    ))
  }
  chunk |>
    dplyr::mutate(stream_key = purrr::map_chr(stream, stream_key)) |>
    dplyr::select(dplyr::all_of(keys), stream_key) |>
    dplyr::bind_cols(purrr::list_rbind(chunk_results))
}

# Run the tasks chunk by chunk (one chunk per sample x relative effect) with
# forked workers. After every chunk the accumulated results go to
# `partial_path` as list(settings, results), so an interrupted run can resume
# from it; a resumed run skips the tasks whose keys are already in
# `completed` (validated beforehand by validate_paper_c_power_partial()).
run_paper_c_power <- function(
  anchors,
  tasks,
  workers,
  partial_path,
  partial_settings,
  completed = NULL
) {
  keys <- c("sample", "dgp", "relative_effect", "replication")
  if (!is.null(completed)) {
    done <- tasks |>
      dplyr::semi_join(completed, by = keys)
    tasks <- tasks |>
      dplyr::anti_join(completed, by = keys)
    message(
      "Resuming: ",
      nrow(done),
      " fits already in ",
      partial_path,
      "; ",
      nrow(tasks),
      " to run"
    )
  }
  results <- completed
  chunk_ids <- unique(tasks$chunk)
  started <- Sys.time()
  for (chunk_index in seq_along(chunk_ids)) {
    chunk <- tasks |>
      dplyr::filter(chunk == chunk_ids[[chunk_index]])
    results <- dplyr::bind_rows(
      results,
      run_paper_c_power_chunk(
        anchors,
        chunk,
        workers,
        failure_note = paste0(
          "Chunk ",
          chunk_index,
          " of ",
          length(chunk_ids),
          "; partial results are in ",
          partial_path
        )
      )
    )
    save_rds_atomic(
      list(settings = partial_settings, results = results),
      partial_path
    )
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    message(
      "Chunk ",
      chunk_index,
      "/",
      length(chunk_ids),
      " (",
      chunk$sample[[1L]],
      ", effect ",
      chunk$relative_effect[[1L]],
      "): ",
      round(elapsed),
      " s elapsed, ",
      nrow(results),
      " fits saved"
    )
  }
  results |>
    dplyr::arrange(sample, dgp, relative_effect, replication)
}

# Per-cell power table.
summarise_paper_c_power <- function(replication_results) {
  replication_results |>
    dplyr::summarise(
      replications = dplyr::n(),
      power = mean(p_value < 0.05, na.rm = TRUE),
      mean_beta_hat = mean(beta_hat, na.rm = TRUE),
      sd_beta_hat = stats::sd(beta_hat, na.rm = TRUE),
      mean_se = mean(se, na.rm = TRUE),
      min_clusters = min(clusters, na.rm = TRUE),
      mean_nobs = mean(nobs, na.rm = TRUE),
      min_nobs = min(nobs, na.rm = TRUE),
      failed = mean(failed),
      .by = c(sample, dgp, relative_effect)
    ) |>
    dplyr::mutate(beta = log1p(relative_effect), .after = relative_effect) |>
    dplyr::arrange(sample, dgp, relative_effect)
}

# Analytic power from the real-data clustered modifier SE: a noncentral t test
# with G - 1 degrees of freedom and noncentrality beta / SE_real.
analytic_paper_c_power <- function(anchors, relative_effects) {
  purrr::imap(anchors, function(anchor, label) {
    real_modifier <- if (is.null(anchor$real_modifier)) {
      anchor$real_alert_high
    } else {
      anchor$real_modifier
    }
    se_real <- real_modifier$se
    clusters <- anchor$n_clusters
    df <- clusters - 1L
    t_crit <- stats::qt(0.975, df)
    beta <- log1p(relative_effects)
    ncp <- beta / se_real
    tibble::tibble(
      sample = label,
      relative_effect = relative_effects,
      beta = beta,
      se_real = se_real,
      clusters = clusters,
      power = stats::pt(-t_crit, df, ncp) + 1 - stats::pt(t_crit, df, ncp)
    )
  }) |>
    purrr::list_rbind()
}

# Real-data clustered SE next to the mean simulated SE, per sample.
compare_paper_c_se <- function(anchors, power_table) {
  simulated <- power_table |>
    dplyr::summarise(mean_simulated_se = mean(mean_se), .by = sample)
  purrr::imap(anchors, function(anchor, label) {
    real_modifier <- if (is.null(anchor$real_modifier)) {
      anchor$real_alert_high
    } else {
      anchor$real_modifier
    }
    tibble::tibble(sample = label, se_real = real_modifier$se)
  }) |>
    purrr::list_rbind() |>
    dplyr::left_join(simulated, by = dplyr::join_by(sample)) |>
    dplyr::mutate(ratio = mean_simulated_se / se_real)
}
