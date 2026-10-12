# NOTE: readable analysis code from the private project. It needs the fire-department
# transport records, which are not distributed, and is not sourced by the public
# pipeline (_targets.R). Kept so that the published estimates can be audited.

# Result 2's baseline for paper C (paper_strategy 7.2) as targets (Issue #50
# (a)): the computation of src/paper_c_result2.R, split into the fits, the
# leave-one-area-out jackknife, the randomization inference in blocks of
# seeds, and the report. The modifiers of 7.2 are deliberately absent: this
# estimates the alert main effect alone, split by regime (nationwide 2021- and
# the 2020 trial).

# Worker count for the forked fits, from PAPER_C_WORKERS. Read inside the
# commands rather than stored in a target, so that changing it never
# invalidates anything; no estimate depends on it.
paper_c_workers <- function(default = 20L) {
  clamp_workers(as.integer(Sys.getenv("PAPER_C_WORKERS", default)))
}

# The run's settings, in the order src/paper_c_result2.R stored them.
# `workers` is filled in when the report is assembled (paper_c_workers()).
paper_c_result2_settings <- function(
  ri_replications = 4999L,
  ri_chunk_size = 200L
) {
  settings <- list(
    bootstrap_replications = 0L,
    ri_replications = ri_replications,
    ri_chunk_size = ri_chunk_size,
    workers = NA_integer_,
    bin = "0.1",
    output = "notes/paper_c_result2_report.rds",
    reuse_partial = FALSE
  )
  settings$seed <- 20260910L
  settings$ri_strata <- c("area_primary", "year", "area_bin_t_1c")
  settings$fixest_version <- as.character(utils::packageVersion("fixest"))
  settings
}

# Randomization inference runs on the two main samples only.
paper_c_result2_ri_labels <- function() c("primary", "stratified")

# The cue of the randomization blocks. A block is never rebuilt once it
# exists, since a small upstream change would otherwise restart a day-long
# run; the report target instead compares every block's stored settings,
# design fingerprint included, with the current inputs and aborts on a
# difference. Rebuild deliberately with PAPER_C_RI_REBUILD=true. A change to
# the permutation code itself is not detected: rebuild after one.
paper_c_ri_cue <- function() {
  targets::tar_cue(
    mode = if (identical(Sys.getenv("PAPER_C_RI_REBUILD"), "true")) {
      "thorough"
    } else {
      "never"
    }
  )
}

# The seed grid of a randomization run in blocks of `chunk_size`, one dynamic
# branch per block. The grid depends on `seed` alone (randomization_seeds()),
# so the blocks reproduce the chunks of the driver exactly.
randomization_seed_blocks <- function(seed, replications, chunk_size) {
  seeds <- randomization_seeds(seed, replications)
  unname(split(seeds, ceiling(seq_along(seeds) / chunk_size)))
}

# One branch: the replications of one block, with the settings they were
# computed under so that the report can refuse a block whose inputs moved.
run_randomization_branch <- function(inputs, seeds, contrasts, workers) {
  list(
    settings = inputs$settings,
    replications = run_randomization_block(inputs, seeds, contrasts, workers)
  )
}

# Bind the blocks in grid order into the object run_randomization_inference()
# returns, after checking that every block belongs to the current inputs and
# that together they hold exactly the seed grid, in order.
assemble_randomization_inference <- function(
  inputs,
  blocks,
  seed,
  replications
) {
  label <- inputs$settings$label
  stale <- !purrr::map_lgl(blocks, \(block) {
    identical(block$settings, inputs$settings)
  })
  if (any(stale)) {
    rlang::abort(glue::glue(
      "{sum(stale)} of {length(blocks)} randomization blocks for {label} were ",
      "computed from different inputs (the blocks are frozen by ",
      "paper_c_ri_cue()). Rebuild them deliberately with ",
      "PAPER_C_RI_REBUILD=true."
    ))
  }
  results <- dplyr::bind_rows(purrr::map(blocks, "replications"))
  # One row per seed and contrast, matched against the grid the inputs imply.
  # The unique seeds alone would pass a block bound twice, whose draws
  # summarise_randomization_inference() would then count twice, and the
  # contrasts must come from the terms rather than from the draws, or a
  # contrast missing from every seed would pass as well.
  seeds <- randomization_seeds(seed, replications)
  expected <- tidyr::expand_grid(
    seed = seeds,
    contrast = names(paper_c_result2_contrasts(inputs$settings$terms))
  )
  draws <- dplyr::select(results, seed, contrast)
  duplicated_draws <- nrow(draws) - nrow(dplyr::distinct(draws))
  missing_draws <- nrow(dplyr::anti_join(
    expected,
    draws,
    by = c("seed", "contrast")
  ))
  unknown_draws <- nrow(dplyr::anti_join(
    draws,
    expected,
    by = c("seed", "contrast")
  ))
  if (
    !identical(unique(results$seed), seeds) ||
      duplicated_draws + missing_draws + unknown_draws > 0L
  ) {
    rlang::abort(glue::glue(
      "The randomization blocks for {label} do not hold the seed grid of ",
      "{replications} replications in order, once per contrast ",
      "({duplicated_draws} duplicated, {missing_draws} missing, ",
      "{unknown_draws} not in the grid)."
    ))
  }
  list(
    strata = summarise_permutation_strata(inputs$area_days),
    strata_vars = inputs$settings$strata_vars,
    frozen_years = inputs$settings$frozen_years,
    replications = results
  )
}

# Gate 13: the two treatment definitions, measured rather than estimated.
# 7.2 lists "either issuance" against "the 5 a.m. bulletin only" as a
# pre-specified sensitivity. In this sample every alerted department-day also
# carries a 5 a.m. bulletin, so the two definitions select the same rows and
# re-estimating under the second would reproduce the first exactly. The
# output names differ from the column names on purpose: summarise() evaluates
# its expressions in order, so reusing `alert_day` as an output name would
# make the comparisons below compare two scalars.
paper_c_result2_gate13 <- function(panel) {
  gate13 <- panel |>
    dplyr::summarise(
      n_alert_days = sum(alert_day),
      n_with_5am = sum(alert_5am),
      n_with_17pm = sum(alert_17pm),
      n_alert_without_5am = sum(alert_day & !alert_5am),
      n_17pm_without_5am = sum(alert_17pm & !alert_5am)
    )
  if (gate13$n_alert_days != gate13$n_with_5am) {
    rlang::abort(glue::glue(
      "The two treatment definitions differ on ",
      "{gate13$n_alert_without_5am} department-days; gate 13 needs a ",
      "sensitivity analysis rather than a count."
    ))
  }
  gate13
}

# The six specifications of the baseline, each a design, its terms and its bin
# resolution.
paper_c_result2_specification_designs <- function(designs, panel, bin) {
  straddling_departments <- panel |>
    dplyr::filter(straddles) |>
    dplyr::distinct(department) |>
    dplyr::pull(department)
  list(
    primary = list(
      design = designs$full,
      terms = c("alert_day", "alert_trial"),
      bin = bin
    ),
    stratified = list(
      design = designs$stratified,
      terms = c("alert_day", "alert_trial"),
      bin = bin
    ),
    exclude_2020 = list(
      design = designs$full |> dplyr::filter(year != 2020L),
      terms = "alert_day",
      bin = bin
    ),
    bins_1c = list(
      design = designs$full,
      terms = c("alert_day", "alert_trial"),
      bin = "1"
    ),
    # Exploratory, not pre-specified: the alert area that carries most of the
    # jackknife variance was picked after seeing the leave-one-out estimates
    # of 2026-09-11. It is written in by name so that a re-run cannot quietly
    # choose a different area and present it as the same specification.
    # 神奈川 is also a 2020 trial area, so dropping it moves the trial
    # coefficients as well.
    exclude_high_leverage = list(
      design = designs$full |> dplyr::filter(area_primary != "神奈川"),
      terms = c("alert_day", "alert_trial"),
      bin = bin
    ),
    exclude_straddling = list(
      design = designs$full |>
        dplyr::filter(!department %in% straddling_departments),
      terms = c("alert_day", "alert_trial"),
      bin = bin
    )
  )
}

# The fits of the six specifications. The fit objects carry the estimation
# data, Saitama department names included, so they are dropped; the kept-row
# index stays for the jackknife.
build_paper_c_result2_fits <- function(analysis_design, panel, settings) {
  gate13 <- paper_c_result2_gate13(panel)
  specification_designs <- paper_c_result2_specification_designs(
    analysis_design$designs,
    panel,
    settings$bin
  )
  specifications <- purrr::imap(specification_designs, \(spec, label) {
    fit <- run_paper_c_result2_spec(
      spec$design,
      label = label,
      terms = spec$terms,
      bin_resolution = spec$bin,
      bootstrap_replications = settings$bootstrap_replications,
      seed = settings$seed
    )
    fit$fit <- NULL
    fit
  })
  list(
    gate13 = gate13,
    straddling_departments = panel |>
      dplyr::filter(straddles) |>
      dplyr::distinct(department) |>
      dplyr::pull(department) |>
      length(),
    specifications = specifications,
    coefficients = purrr::map(specifications, "coefficients") |>
      purrr::list_rbind()
  )
}

# Leave-one-area-out jackknife (CV3) of the two main samples. The interval
# comes from here: randomization inference tests the sharp null but does not
# invert into an interval for a count outcome, and CR1's own interval rests
# on the variance estimator 7.12 showed to be biased down.
build_paper_c_result2_jackknife <- function(
  analysis_design,
  panel,
  fits,
  settings
) {
  specification_designs <- paper_c_result2_specification_designs(
    analysis_design$designs,
    panel,
    settings$bin
  )
  labels <- paper_c_result2_ri_labels()
  purrr::map(rlang::set_names(labels), \(label) {
    spec <- fits$specifications[[label]]
    design_kept <- specification_designs[[label]]$design[
      spec$design_kept_rows,
    ]
    fit <- paper_c_jackknife_vcov(
      design_kept,
      spec$terms,
      spec$bin_resolution,
      workers = paper_c_workers()
    )
    fit$summary <- summarise_jackknife(
      fit,
      spec$estimates,
      paper_c_result2_contrasts(spec$terms)
    )
    fit
  })
}

# Inputs of the randomization runs of the two main samples.
build_paper_c_result2_ri_inputs <- function(
  analysis_design,
  panel,
  fits,
  settings
) {
  specification_designs <- paper_c_result2_specification_designs(
    analysis_design$designs,
    panel,
    settings$bin
  )
  labels <- paper_c_result2_ri_labels()
  purrr::map(rlang::set_names(labels), \(label) {
    spec <- fits$specifications[[label]]
    prepare_randomization_inputs(
      specification_designs[[label]]$design,
      terms = spec$terms,
      bin_resolution = spec$bin_resolution,
      replications = settings$ri_replications,
      seed = settings$seed,
      strata_vars = settings$ri_strata,
      label = label
    )
  })
}

# The report src/paper_c_result2.R wrote, from the targets above. `timings`
# stays empty: the wall times are in tar_meta().
build_paper_c_result2_report <- function(
  analysis_design,
  panel,
  fits,
  jackknife,
  ri_inputs,
  ri_blocks,
  settings
) {
  randomization <- purrr::imap(ri_inputs, \(inputs, label) {
    draws <- assemble_randomization_inference(
      inputs,
      ri_blocks[[label]],
      settings$seed,
      settings$ri_replications
    )
    draws$summary <- summarise_randomization_inference(
      draws$replications,
      fits$specifications[[label]]$coefficients
    )
    draws
  })
  report <- list(
    design_summary = analysis_design$summary,
    gate13 = fits$gate13,
    straddling_departments = fits$straddling_departments,
    specifications = purrr::map(fits$specifications, strip_fit),
    jackknife = jackknife,
    coefficients = fits$coefficients,
    randomization = randomization,
    settings = c(
      utils::modifyList(settings, list(workers = paper_c_workers())),
      list(timings = list())
    ),
    session = paper_c_session()
  )
  check_no_saitama(report, saitama_departments(panel), "The report")
  report
}

# Nationwide-only randomization (paper_strategy 7.13) --------------------------

# The randomization of the baseline permutes every alert label, so the null
# it tests is the global one. Holding 2020's assignment at its observed value
# and permuting only the nationwide years leaves the trial coefficient
# estimated from real data, so what is perturbed under the null is the
# nationwide coefficient by itself. Settings in the order
# src/paper_c_result2_nationwide.R stored them; `workers` is filled in when
# the result is assembled.
paper_c_result2_nationwide_settings <- function(
  ri_replications = 4999L,
  ri_chunk_size = 200L
) {
  list(
    ri_replications = ri_replications,
    ri_chunk_size = ri_chunk_size,
    workers = NA_integer_,
    frozen_year = 2020L,
    report = "notes/paper_c_result2_report.rds",
    output = "notes/paper_c_result2_nationwide.rds",
    seed = 20260912L,
    ri_strata = c("area_primary", "year", "area_bin_t_1c")
  )
}

# The baseline's two samples under the names of the specifications that use
# them.
paper_c_result2_ri_designs <- function(analysis_design) {
  list(
    primary = analysis_design$designs$full,
    stratified = analysis_design$designs$stratified
  )
}

build_paper_c_result2_nationwide_ri_inputs <- function(
  analysis_design,
  fits,
  settings
) {
  designs <- paper_c_result2_ri_designs(analysis_design)
  purrr::imap(designs, \(design, label) {
    spec <- fits$specifications[[label]]
    prepare_randomization_inputs(
      design,
      terms = spec$terms,
      bin_resolution = spec$bin_resolution,
      replications = settings$ri_replications,
      seed = settings$seed,
      strata_vars = settings$ri_strata,
      frozen_years = settings$frozen_year,
      label = label
    )
  })
}

# What the p-values are computed against: the baseline's observed statistics,
# confirmed by refitting the model on the design (two fits, about a minute),
# so that the result cannot be paired with a different run's statistics.
build_paper_c_result2_nationwide_provenance <- function(analysis_design, fits) {
  designs <- paper_c_result2_ri_designs(analysis_design)
  purrr::imap(designs, \(design, label) {
    spec <- fits$specifications[[label]]
    stored_hash <- spec$design_hash
    design_hash <- paper_c_design_hash(design)
    hash_checked <- !is.null(stored_hash)
    if (hash_checked && !identical(stored_hash, design_hash)) {
      rlang::abort(glue::glue(
        "The design for {label} differs from the one the fits were built on."
      ))
    }
    refit <- run_paper_c_result2_spec(
      design,
      label = label,
      terms = spec$terms,
      bin_resolution = spec$bin_resolution,
      bootstrap_replications = 0L
    )
    stored <- fits$coefficients |>
      dplyr::filter(specification == label) |>
      dplyr::arrange(contrast)
    fresh <- refit$coefficients |> dplyr::arrange(contrast)
    if (!identical(stored$contrast, fresh$contrast)) {
      rlang::abort(glue::glue("Contrasts differ for {label}."))
    }
    gap <- max(abs(stored$statistic - fresh$statistic))
    if (gap > 1e-4) {
      rlang::abort(glue::glue(
        "Refitting {label} reproduces the fits' statistics only to {gap}."
      ))
    }
    list(
      label = label,
      design_hash = design_hash,
      hash_checked = hash_checked,
      statistic_gap = gap,
      observed = stored |> dplyr::select(contrast, estimate, se, statistic)
    )
  })
}

# The result src/paper_c_result2_nationwide.R wrote. There is no `seconds`
# or `verified_at`: the wall time is in tar_meta(), and the provenance is
# computed with the result rather than stamped on afterwards.
build_paper_c_result2_nationwide <- function(
  panel,
  fits,
  ri_inputs,
  ri_blocks,
  provenance,
  settings
) {
  randomization <- purrr::imap(ri_inputs, \(inputs, label) {
    draws <- assemble_randomization_inference(
      inputs,
      ri_blocks[[label]],
      settings$seed,
      settings$ri_replications
    )
    draws$summary <- summarise_randomization_inference(
      draws$replications,
      fits$coefficients |> dplyr::filter(specification == label)
    )
    draws
  })
  result <- list(
    randomization = randomization,
    settings = utils::modifyList(settings, list(workers = paper_c_workers())),
    session = paper_c_session("fixest"),
    provenance = provenance
  )
  check_no_saitama(result, saitama_departments(panel), "The result")
  result
}

# CR1-only sensitivity specifications (Issue #43 items 1, 4, 6, 8, 11) ------------

# The specifications of src/paper_c_result2_sensitivity.R. Before any of them,
# the primary specification is refitted and required to reproduce the
# baseline fits exactly, so every coefficient is read against the baseline the
# manuscript prints. No randomization inference and no jackknife.
build_paper_c_sensitivity_report <- function(
  analysis_design,
  panel,
  fits,
  result2_settings,
  fd_merge_pop2020,
  covid_emergency_file,
  covid_priority_file
) {
  settings <- list(
    bin = "0.1",
    output = "notes/paper_c_sensitivity_report.rds",
    baseline = "notes/paper_c_result2_report.rds"
  )
  settings$seed <- 20260910L
  settings$fixest_version <- as.character(utils::packageVersion("fixest"))
  designs <- analysis_design$designs

  # Baseline equivalence.
  if (!identical(result2_settings$fixest_version, settings$fixest_version)) {
    rlang::abort(glue::glue(
      "The baseline was fitted with fixest ",
      "{result2_settings$fixest_version}; this session has ",
      "{settings$fixest_version}."
    ))
  }
  baseline_stored <- as.data.frame(fits$coefficients) |>
    dplyr::filter(specification == "primary") |>
    dplyr::select(contrast, estimate_stored = estimate, se_stored = se) |>
    dplyr::arrange(contrast)
  primary <- run_paper_c_result2_spec(
    designs$full,
    label = "primary",
    terms = c("alert_day", "alert_trial"),
    bin_resolution = settings$bin,
    seed = settings$seed
  )
  baseline_current <- as.data.frame(primary$coefficients) |>
    dplyr::select(contrast, estimate_current = estimate, se_current = se) |>
    dplyr::arrange(contrast)
  baseline_equivalence <- baseline_stored |>
    dplyr::inner_join(baseline_current, by = dplyr::join_by(contrast)) |>
    dplyr::mutate(
      estimate_difference = estimate_current - estimate_stored,
      se_difference = se_current - se_stored
    )
  if (
    nrow(baseline_equivalence) != nrow(baseline_stored) ||
      nrow(baseline_equivalence) != nrow(baseline_current) ||
      any(baseline_equivalence$estimate_difference != 0) ||
      any(baseline_equivalence$se_difference != 0) ||
      primary$rows_estimated != fits$specifications$primary$rows_estimated
  ) {
    rlang::abort(
      "The refitted primary specification does not reproduce the baseline."
    )
  }

  # Blank records: keys of the panel rows whose transport count is NA (the
  # panel builder asserts that these are fully blank records: 378 of them in
  # the 25-prefecture sample, 338 in the pre-specified 23; Issue #59).
  blank_keys <- panel |>
    dplyr::filter(is.na(transports)) |>
    dplyr::distinct(prefecture, date, department)
  stopifnot(nrow(blank_keys) == 378L)
  # Lead and lag of the area-day alert. The panel covers June-September, so
  # the lead of 30 September and the lag of 1 June fall outside it; those
  # edges are set to 0 and counted.
  area_alert <- designs$full |>
    dplyr::distinct(area_primary, date, alert_day)
  stopifnot(
    nrow(area_alert) == nrow(dplyr::distinct(area_alert, area_primary, date))
  )
  lead_lag <- area_alert |>
    dplyr::left_join(
      area_alert |>
        dplyr::transmute(
          area_primary,
          date = date - 1L,
          alert_lead1 = alert_day
        ),
      by = dplyr::join_by(area_primary, date)
    ) |>
    dplyr::left_join(
      area_alert |>
        dplyr::transmute(
          area_primary,
          date = date + 1L,
          alert_lag1 = alert_day
        ),
      by = dplyr::join_by(area_primary, date)
    )
  edge_counts <- c(
    lead_missing = sum(is.na(lead_lag$alert_lead1)),
    lag_missing = sum(is.na(lead_lag$alert_lag1))
  )
  lead_lag <- lead_lag |>
    dplyr::mutate(
      alert_lead1 = tidyr::replace_na(alert_lead1, 0),
      alert_lag1 = tidyr::replace_na(alert_lag1, 0)
    ) |>
    dplyr::select(area_primary, date, alert_lead1, alert_lag1)
  # Issue #43 item 11: COVID-19 emergency declarations (prefecture-day) and
  # priority measures (department-day through the municipality crosswalk),
  # June-September 2021 only; 2020 and 2022 have none inside the window.
  crosswalk <- fd_merge_pop2020 |>
    dplyr::distinct(prefecture, department, administration) |>
    dplyr::semi_join(
      panel |> dplyr::distinct(department),
      by = dplyr::join_by(department)
    )
  stopifnot(
    dplyr::n_distinct(crosswalk$department) ==
      dplyr::n_distinct(panel$department)
  )
  covid <- build_covid_indicators(
    emergency = readr::read_csv(covid_emergency_file, show_col_types = FALSE),
    priority = readr::read_csv(covid_priority_file, show_col_types = FALSE),
    crosswalk = crosswalk,
    dates = sort(unique(designs$full$date))
  )
  specification_designs <- list(
    no_local_bins = list(
      design = designs$full,
      terms = c("alert_day", "alert_trial"),
      drop_local = TRUE
    ),
    # The primary and no-local-bins fits keep different rows (the local bins
    # create singleton levels the other model lacks), so the comparison the
    # manuscript reads is made on the rows the primary specification
    # estimates.
    no_local_bins_common = list(
      design = designs$full[primary$design_kept_rows, ],
      terms = c("alert_day", "alert_trial"),
      drop_local = TRUE
    ),
    years_2021_2022 = list(
      design = designs$full |> dplyr::filter(year >= 2021L),
      terms = "alert_day",
      drop_local = FALSE
    ),
    # Issue #43 item 8: the 338 department-days whose record is blank in every
    # column are read as zero-transport days in the primary sample (a declared
    # imputation, Sample section). Dropping them instead.
    drop_blank_records = list(
      design = designs$full |>
        dplyr::anti_join(
          blank_keys,
          by = dplyr::join_by(prefecture, date, department)
        ),
      terms = c("alert_day", "alert_trial"),
      drop_local = FALSE
    ),
    # Issue #43 item 6: the alert on the previous day (alert_lag1; calls
    # brought forward would show as a deficit on the day after an alert) and
    # on the next day (alert_lead1; a placebo, with the caveat that the 17:00
    # bulletin announcing the next day's alert acts on the evening before).
    # Both enter one model with the primary terms.
    alert_lead_lag = list(
      design = designs$full |>
        dplyr::left_join(lead_lag, by = dplyr::join_by(area_primary, date)),
      terms = c("alert_day", "alert_trial", "alert_lead1", "alert_lag1"),
      drop_local = FALSE
    ),
    covid_controls = list(
      design = designs$full |>
        dplyr::left_join(
          covid$indicators,
          by = dplyr::join_by(prefecture, department, date)
        ),
      terms = c("alert_day", "alert_trial", "emergency", "priority"),
      drop_local = FALSE
    ),
    # Issue #59: the pre-specified 23-prefecture sample, refitted from the
    # current design. Its rows, area WBGT and fixed effects are those of the
    # sample as disclosed, so this fit is the current code's reading of the
    # results frozen in data-raw/frozen/paper_c_p23; the jackknife below
    # gives its CV3J interval.
    prespecified_23 = list(
      design = designs$full |>
        dplyr::filter(prefecture %in% paper_c_prespecified_prefectures()),
      terms = c("alert_day", "alert_trial"),
      drop_local = FALSE
    )
  )
  covid_design <- specification_designs$covid_controls$design
  stopifnot(
    nrow(covid_design) == nrow(designs$full),
    !anyNA(covid_design$emergency),
    !anyNA(covid_design$priority),
    all(
      covid_design$year[
        covid_design$emergency == 1 | covid_design$priority == 1
      ] ==
        2021L
    ),
    nrow(specification_designs$drop_blank_records$design) ==
      nrow(designs$full) - nrow(blank_keys),
    nrow(specification_designs$alert_lead_lag$design) == nrow(designs$full),
    !anyNA(specification_designs$alert_lead_lag$design$alert_lead1),
    !anyNA(specification_designs$alert_lead_lag$design$alert_lag1),
    nrow(specification_designs$no_local_bins_common$design) ==
      primary$rows_estimated,
    setequal(specification_designs$years_2021_2022$design$year, c(2021, 2022)),
    identical(
      sort(unique(designs$full$alert_day)),
      sort(unique(specification_designs$alert_lead_lag$design$alert_lead1))
    ),
    all(specification_designs$years_2021_2022$design$alert_trial == 0)
  )
  specifications <- purrr::imap(specification_designs, \(spec, label) {
    run_paper_c_result2_spec(
      spec$design,
      label = label,
      terms = spec$terms,
      bin_resolution = settings$bin,
      seed = settings$seed,
      drop_local = spec$drop_local
    )
  })
  specifications$primary <- primary
  coefficients <- purrr::map(specifications, "coefficients") |>
    purrr::list_rbind()
  prespecified <- specifications$prespecified_23
  prespecified_jackknife <- paper_c_jackknife_vcov(
    specification_designs$prespecified_23$design[
      prespecified$design_kept_rows,
    ],
    prespecified$terms,
    prespecified$bin_resolution,
    workers = paper_c_workers()
  )
  prespecified_jackknife$summary <- summarise_jackknife(
    prespecified_jackknife,
    prespecified$estimates,
    paper_c_result2_contrasts(prespecified$terms)
  )
  # The no-local-bins fit must have exactly the local bins fewer fixed
  # effects and the same rows going in; the 2021-2022 fit must keep only
  # those years.
  stopifnot(
    identical(specifications$no_local_bins$rows, specifications$primary$rows),
    !"local_bin" %in% names(specifications$no_local_bins$removals$fixef_sizes),
    "local_bin" %in% names(specifications$primary$removals$fixef_sizes),
    specifications$years_2021_2022$rows ==
      nrow(specification_designs$years_2021_2022$design)
  )

  report <- list(
    design_summary = analysis_design$summary,
    baseline_equivalence = baseline_equivalence,
    blank_records = nrow(blank_keys),
    lead_lag_edges = edge_counts,
    covid = list(
      source = paper_c_covid19_source(),
      emergency_department_days = sum(covid$indicators$emergency),
      priority_department_days = sum(covid$indicators$priority),
      unmatched = covid$unmatched
    ),
    specifications = purrr::map(specifications, strip_fit),
    coefficients = coefficients,
    prespecified_23_jackknife = prespecified_jackknife,
    settings = c(settings, list(timings = list())),
    session = paper_c_session()
  )
  check_no_saitama(report, saitama_departments(panel), "The sensitivity report")
  report
}

# Alert-area-day concordance (paper C, Issue #45 item A), moved from
# src/paper_c_area_day_concordance.R (Issue #50 (a)). For system-active
# area-days: how often an alert day is an area-day whose realized maximum
# reached 33 degC and how often such an area-day carried an alert. It also
# reproduces the permutation-strata summary of result 2's primary
# randomization (from its inputs, which the report's summary is computed
# from) and counts the unalerted days the permutation exchanges with alerted
# ones by their alert status. Only area-level aggregates are returned.
build_paper_c_area_day_concordance <- function(
  analysis_design,
  panel,
  ri_inputs_primary
) {
  design <- analysis_design$designs$full
  statuses <- c("alert", "active, no alert", "structural no-alert")

  # Area-days of the estimation design ------------------------------------------
  # alert_day and the area maximum are area-level quantities joined onto
  # department rows; both must be constant within an area-day.
  area_days <- design |>
    dplyr::distinct(
      area_primary,
      date,
      year,
      alert_day,
      area_wbgt_max,
      area_bin_t_1c
    )
  stopifnot(
    nrow(area_days) == nrow(dplyr::distinct(area_days, area_primary, date)),
    all(is.finite(area_days$area_wbgt_max))
  )
  status <- panel |>
    dplyr::distinct(area_primary, date, system_active, alert_status)
  stopifnot(nrow(status) == nrow(dplyr::distinct(status, area_primary, date)))
  area_days <- area_days |>
    dplyr::inner_join(status, by = dplyr::join_by(area_primary, date)) |>
    dplyr::mutate(
      area_high = area_wbgt_max >= 33,
      alert_status = factor(alert_status, levels = statuses)
    )
  stopifnot(
    nrow(area_days) == nrow(dplyr::distinct(design, area_primary, date)),
    all((area_days$alert_status == "alert") == (area_days$alert_day == 1))
  )
  # Area-days present in the panel but absent from the design: every department
  # of the area lacked local WBGT that day (the design drops NA local_high).
  panel_area_days <- panel |>
    dplyr::distinct(area_primary, date)
  absent <- dplyr::anti_join(
    panel_area_days,
    area_days,
    by = dplyr::join_by(area_primary, date)
  )

  # Cells over system-active area-days, overall and by year ---------------------
  cells <- function(x, ...) {
    x |>
      dplyr::filter(system_active) |>
      dplyr::summarise(
        area_days = dplyr::n(),
        alert_days = sum(alert_day),
        high_days = sum(area_high),
        alert_and_high = sum(alert_day == 1 & area_high),
        alert_not_high = sum(alert_day == 1 & !area_high),
        no_alert_high = sum(alert_day == 0 & area_high),
        no_alert_not_high = sum(alert_day == 0 & !area_high),
        .by = c(...)
      ) |>
      dplyr::mutate(
        share_high_given_alert = alert_and_high / alert_days,
        share_alert_given_high = alert_and_high / high_days,
        share_alert_given_not_high = alert_not_high / (area_days - high_days)
      )
  }
  overall <- cells(area_days) |> dplyr::mutate(year = NA_integer_, .before = 1)
  by_year <- cells(area_days, year) |> dplyr::arrange(year)
  stopifnot(
    sum(by_year$area_days) == overall$area_days,
    overall$alert_and_high +
      overall$alert_not_high +
      overall$no_alert_high +
      overall$no_alert_not_high ==
      overall$area_days,
    overall$alert_days == sum(area_days$alert_day)
  )

  # Permutation strata: reproduce the stored summary and classify the exchanged
  # unalerted days by status ----------------------------------------------------
  stored_strata <- as.data.frame(
    summarise_permutation_strata(ri_inputs_primary$area_days)
  )
  frozen_years <- ri_inputs_primary$settings$frozen_years
  strata <- build_alert_permutation_strata(
    design,
    strata_vars = ri_inputs_primary$settings$strata_vars,
    frozen_years = frozen_years
  )
  recomputed <- summarise_permutation_strata(strata)
  required <- c(
    "strata",
    "mixed_strata",
    "alerted_area_days",
    "alerted_in_mixed_strata"
  )
  missing <- setdiff(required, names(stored_strata))
  if (length(missing) > 0L) {
    rlang::abort(glue::glue(
      "The stored strata summary lacks {paste(missing, collapse = ', ')}."
    ))
  }
  for (column in intersect(names(recomputed), names(stored_strata))) {
    if (
      !identical(
        as.numeric(recomputed[[column]]),
        as.numeric(stored_strata[[column]])
      )
    ) {
      rlang::abort(glue::glue(
        "Permutation strata differ from the stored report in `{column}`: ",
        "{recomputed[[column]]} versus {stored_strata[[column]]}."
      ))
    }
  }
  mixed_ids <- strata |>
    dplyr::filter(!frozen) |>
    dplyr::summarise(alerted = sum(alert_day), n = dplyr::n(), .by = stratum) |>
    dplyr::filter(alerted > 0, alerted < n) |>
    dplyr::pull(stratum)
  exchanged_unalerted <- strata |>
    dplyr::filter(!frozen, stratum %in% mixed_ids, alert_day == 0) |>
    dplyr::inner_join(status, by = dplyr::join_by(area_primary, date)) |>
    dplyr::count(alert_status, name = "area_days") |>
    dplyr::mutate(alert_status = factor(alert_status, levels = statuses)) |>
    dplyr::arrange(alert_status)
  stopifnot(
    !"alert" %in% as.character(exchanged_unalerted$alert_status),
    sum(exchanged_unalerted$area_days) ==
      sum(
        strata$stratum %in% mixed_ids & !strata$frozen & strata$alert_day == 0
      )
  )

  # Output ----------------------------------------------------------------------
  concordance <- dplyr::bind_rows(overall, by_year)
  meta <- tibble::tibble(
    design_area_days = nrow(area_days),
    panel_area_days = nrow(panel_area_days),
    absent_area_days = nrow(absent),
    system_active_area_days = overall$area_days,
    areas = dplyr::n_distinct(area_days$area_primary),
    strata = recomputed$strata,
    mixed_strata = recomputed$mixed_strata,
    alerted_area_days = recomputed$alerted_area_days,
    alerted_in_mixed_strata = recomputed$alerted_in_mixed_strata,
    frozen_years = paste(frozen_years, collapse = ";"),
    frozen_alerted_area_days = recomputed$frozen_alerted_area_days,
    exchanged_unalerted_active = exchanged_unalerted$area_days[
      exchanged_unalerted$alert_status == "active, no alert"
    ],
    exchanged_unalerted_structural = sum(exchanged_unalerted$area_days[
      exchanged_unalerted$alert_status == "structural no-alert"
    ])
  )
  stopifnot(nrow(meta) == 1L, !anyNA(meta))
  list(concordance = concordance, meta = meta)
}
