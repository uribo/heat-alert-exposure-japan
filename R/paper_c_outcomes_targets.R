# NOTE: readable analysis code from the private project. It needs the fire-department
# transport records, which are not distributed, and is not sourced by the public
# pipeline (_targets.R). Kept so that the published estimates can be audited.

# Result 2's five pre-specified severity and place outcomes (paper_strategy
# 7.15) as targets (Issue #50 (a)): the computation of src/paper_c_outcomes.R.
# The driver ran a smoke stage (fits, jackknife, B = 20) and a production stage
# (B = 999) that reused the smoke fits after a design-hash check; here both are
# one pipeline and the production randomization runs in blocks of seeds.

paper_c_outcome_names <- function() {
  c("y_moderate_plus", "y_mild", "y_home", "y_indoor_nonhome", "y_outdoor")
}

# Settings of the production run, in the driver's order. The paths are the
# repository-relative ones the stored report carried; `workers` is filled in
# when the report is assembled.
paper_c_outcomes_settings <- function(
  ri_replications = 999L,
  ri_chunk_size = 100L
) {
  settings <- list(
    ri_replications = ri_replications,
    ri_chunk_size = ri_chunk_size,
    workers = NA_integer_,
    bin = "0.1",
    output = "notes/paper_c_outcomes_report.rds",
    reuse_report = "notes/paper_c_outcomes_report.rds",
    stage = "production",
    jackknife = TRUE
  )
  settings$seed <- 20260914L
  settings$ri_strata <- c("area_primary", "year", "area_bin_t_1c")
  settings$fixest_version <- as.character(utils::packageVersion("fixest"))
  settings$terms <- c("alert_day", "alert_trial")
  settings
}

# The stratified design with the five outcome columns joined on, and the join
# checks.
build_paper_c_outcomes_design <- function(
  analysis_design,
  panel,
  panel_outcomes
) {
  outcomes <- paper_c_outcome_names()
  design_base <- analysis_design$designs$stratified
  outcome_columns <- panel |>
    dplyr::select(
      prefecture,
      date,
      department,
      transports_panel = transports
    ) |>
    dplyr::left_join(
      panel_outcomes |>
        dplyr::select(
          prefecture,
          date,
          department,
          transports_outcomes = transports,
          dplyr::all_of(outcomes)
        ),
      by = dplyr::join_by(prefecture, date, department),
      relationship = "one-to-one"
    )
  design <- design_base |>
    dplyr::left_join(
      outcome_columns,
      by = dplyr::join_by(prefecture, date, department),
      relationship = "one-to-one"
    )
  join_checks <- tibble::tibble(
    check = c(
      "joined rows unchanged",
      "expected stratified rows",
      "outcomes complete",
      "transports agree"
    ),
    observed = c(
      nrow(design),
      nrow(design),
      sum(is.na(design[, outcomes, drop = FALSE])),
      sum(
        !((is.na(design$transports_panel) &
          is.na(design$transports_outcomes)) |
          (!is.na(design$transports_panel) &
            !is.na(design$transports_outcomes) &
            design$transports_panel == design$transports_outcomes))
      )
    ),
    expected = c(nrow(design_base), 142996L, 0L, 0L)
  ) |>
    dplyr::mutate(passed = observed == expected)
  if (!all(join_checks$passed)) {
    rlang::abort(
      "The outcome join failed its row, completeness, or transport check."
    )
  }
  list(design = design, join_checks = join_checks)
}

# The checks that precede the fits: the synthetic multiplicity and finiteness
# checks, the all-transport baseline against the tracked outcome scan and the
# baseline fits, the shared permutation fingerprints, and the shared labels.
build_paper_c_outcomes_checks <- function(
  outcomes_design,
  result2_fits,
  scan_file,
  settings
) {
  outcomes <- paper_c_outcome_names()
  terms <- settings$terms
  design <- outcomes_design$design
  wy_checks <- check_westfall_young()
  finiteness_checks <- check_randomization_finiteness()
  baseline_explicit <- run_paper_c_result2_spec(
    design,
    "stratified",
    terms = terms,
    bin_resolution = settings$bin,
    outcome = "transports_analysis"
  )

  # The tracked scan ran on the pre-specified 23-prefecture sample, so it is
  # compared with a refit on those prefectures' rows; since Issue #59 the
  # design holds 富山 and 沖縄 as well. Their rows, area WBGT and fixed effects
  # are unchanged, so the refit reproduces the scan to 1e-10 or the design of
  # the pre-specified rows has moved.
  scan <- readRDS(scan_file)$outcomes
  scan_selected <- scan |>
    dplyr::filter(outcome %in% c("y_all", outcomes))
  scan_baseline <- scan_selected |> dplyr::filter(outcome == "y_all")
  baseline_scan <- run_paper_c_result2_spec(
    design |>
      dplyr::filter(prefecture %in% paper_c_prespecified_prefectures()),
    "stratified",
    terms = terms,
    bin_resolution = settings$bin,
    outcome = "transports_analysis"
  )
  baseline_fixest_se <- fixest::se(baseline_scan$fit)
  baseline_counts <- tibble::tibble(
    check = c("rows", "departments", "dates", "clusters"),
    observed = c(
      baseline_scan$rows_estimated,
      baseline_scan$removals$fixef_sizes[["department"]],
      baseline_scan$removals$fixef_sizes[["date"]],
      baseline_scan$clusters
    ),
    expected = c(125248L, 256L, 610L, 21L)
  )
  baseline_se <- baseline_scan$coefficients |>
    dplyr::filter(contrast %in% c("alert_day", "alert_trial")) |>
    dplyr::select(contrast, reported = se) |>
    dplyr::mutate(observed = baseline_fixest_se[contrast]) |>
    dplyr::left_join(
      tibble::tibble(
        contrast = c("alert_day", "alert_trial"),
        expected = c(scan_baseline$se_alert, scan_baseline$se_trial),
        tolerance = 1e-10
      ),
      by = dplyr::join_by(contrast),
      relationship = "one-to-one"
    ) |>
    dplyr::mutate(
      check = paste0("baseline SE ", contrast),
      difference = observed - expected,
      reported_difference = reported - expected,
      passed = abs(difference) <= tolerance & abs(reported_difference) <= 1e-8
    )
  baseline_counts <- baseline_counts |>
    dplyr::mutate(difference = observed - expected, passed = difference == 0)
  if (!all(baseline_counts$passed) || !all(baseline_se$passed)) {
    rlang::abort(glue::glue(
      "The all-transport baseline differs from the tracked scan: count gaps ",
      "{paste(baseline_counts$difference, collapse = '/')}; SE gaps ",
      "{paste(signif(baseline_se$difference, 8), collapse = '/')}; reported ",
      "CR1 gaps {paste(signif(baseline_se$reported_difference, 8), collapse = '/')} ."
    ))
  }

  # The baseline fits are the result-2 report's source, so they are always
  # available here.
  stored <- result2_fits$coefficients |>
    dplyr::filter(specification == "stratified") |>
    dplyr::arrange(contrast)
  current <- baseline_explicit$coefficients |>
    dplyr::arrange(contrast)
  if (
    !identical(stored$estimate, current$estimate) ||
      !identical(stored$se, current$se)
  ) {
    rlang::abort("The result-2 baseline fits differ from the baseline refit.")
  }
  stored_baseline <- tibble::tibble(
    check = "stored result-2 report",
    status = "identical",
    passed = TRUE
  )

  # Common design fingerprints; direct label equality below is the
  # cross-outcome permutation check. Hash equality is true by construction.
  hash_columns <- setdiff(names(design), c(outcomes, "transports_analysis"))
  common_hashes <- list(
    design_hash = paper_c_design_hash(design, columns = hash_columns),
    strata_hash = rlang::hash(build_alert_permutation_strata(
      design,
      settings$ri_strata
    )),
    seeds_hash = rlang::hash(randomization_seeds(
      settings$seed,
      settings$ri_replications
    ))
  )
  outcome_hashes <- stats::setNames(
    lapply(outcomes, function(outcome) {
      outcome_design <- design
      outcome_design$transports_analysis <- outcome_design[[outcome]]
      list(
        design_hash = paper_c_design_hash(
          outcome_design,
          columns = hash_columns
        ),
        strata_hash = rlang::hash(build_alert_permutation_strata(
          outcome_design,
          settings$ri_strata
        )),
        seeds_hash = rlang::hash(randomization_seeds(
          settings$seed,
          settings$ri_replications
        ))
      )
    }),
    outcomes
  )
  shared_permutation_checks <- check_shared_permutations(outcome_hashes)
  if (!identical(common_hashes, outcome_hashes[[1L]])) {
    rlang::abort(
      "The common permutation hashes changed inside an outcome branch."
    )
  }

  label_seeds <- randomization_seeds(settings$seed, settings$ri_replications)[
    seq_len(min(3L, settings$ri_replications))
  ]
  label_checks <- purrr::map_dfr(label_seeds, function(replication_seed) {
    labels <- lapply(outcomes, function(outcome) {
      outcome_design <- design
      outcome_design$transports_analysis <- outcome_design[[outcome]]
      area_days <- build_alert_permutation_strata(
        outcome_design,
        settings$ri_strata
      )
      permute_alert_design(
        outcome_design,
        area_days,
        replication_seed,
        terms
      )$alert_day
    })
    tibble::tibble(
      seed = replication_seed,
      passed = all(vapply(labels[-1L], identical, logical(1), labels[[1L]]))
    )
  })
  if (!all(label_checks$passed)) {
    rlang::abort("Permuted alert labels differ across outcomes.")
  }

  list(
    baseline_equivalence = list(
      counts = baseline_counts,
      standard_errors = baseline_se,
      stored_report = stored_baseline
    ),
    shared_permutation_checks = shared_permutation_checks,
    shared_permutation_hashes = common_hashes,
    outcome_columns_hash = rlang::hash(design[, outcomes, drop = FALSE]),
    permutation_label_checks = label_checks,
    westfall_young_checks = wy_checks,
    finiteness_checks = finiteness_checks,
    scan_selected = scan_selected
  )
}

# The five outcome fits and their equivalence with the tracked scan. The fit
# objects are dropped; the kept-row index stays for the jackknife.
build_paper_c_outcomes_fits <- function(outcomes_design, checks, settings) {
  outcomes <- paper_c_outcome_names()
  design <- outcomes_design$design
  specifications <- purrr::map(outcomes, function(outcome) {
    result <- run_paper_c_result2_spec(
      design,
      label = outcome,
      terms = settings$terms,
      bin_resolution = settings$bin,
      outcome = outcome
    )
    result$fit <- NULL
    result
  }) |>
    stats::setNames(outcomes)
  coefficients <- purrr::imap_dfr(specifications, function(x, outcome) {
    x$coefficients |> dplyr::mutate(outcome = outcome, .before = 1)
  })
  # The tracked scan ran on the pre-specified 23 prefectures, so each outcome
  # is compared on a refit of those prefectures' rows (Issue #59); `total` is
  # the current sample's.
  scan_selected <- checks$scan_selected
  design_prespecified <- design |>
    dplyr::filter(prefecture %in% paper_c_prespecified_prefectures())
  outcome_equivalence <- purrr::map_dfr(outcomes, function(outcome) {
    spec <- run_paper_c_result2_spec(
      design_prespecified,
      label = outcome,
      terms = settings$terms,
      bin_resolution = settings$bin,
      outcome = outcome
    )
    scan_row <- scan_selected |>
      dplyr::filter(.data$outcome == .env$outcome)
    expected_total <- c(
      y_moderate_plus = 30931,
      y_mild = 60880,
      y_home = 37352,
      y_indoor_nonhome = 23099,
      y_outdoor = 26530
    )[[outcome]]
    observed_se <- spec$coefficients |>
      dplyr::filter(contrast == "alert_day") |>
      dplyr::pull(se)
    prespecified_total <- sum(design_prespecified[[outcome]])
    tibble::tibble(
      outcome = outcome,
      total = sum(design[[outcome]]),
      prespecified_total = prespecified_total,
      expected_total = expected_total,
      scan_total = scan_row$total,
      se_alert = observed_se,
      scan_se_alert = scan_row$se_alert,
      total_passed = prespecified_total == expected_total &&
        prespecified_total == scan_total,
      se_passed = abs(se_alert - scan_se_alert) <= 1e-8
    )
  })
  if (!all(outcome_equivalence$total_passed & outcome_equivalence$se_passed)) {
    rlang::abort(
      "An outcome total or alert-day SE differs from the tracked scan."
    )
  }
  list(
    specifications = specifications,
    coefficients = coefficients,
    outcome_equivalence = outcome_equivalence
  )
}

build_paper_c_outcomes_jackknife <- function(outcomes_design, fits) {
  outcomes <- paper_c_outcome_names()
  design <- outcomes_design$design
  purrr::map(outcomes, function(outcome) {
    spec <- fits$specifications[[outcome]]
    result <- paper_c_jackknife_vcov(
      design[spec$design_kept_rows, ],
      spec$terms,
      spec$bin_resolution,
      workers = paper_c_workers(),
      outcome = outcome
    )
    result$summary <- summarise_jackknife(
      result,
      spec$estimates,
      paper_c_result2_contrasts(spec$terms)
    )
    result
  }) |>
    stats::setNames(outcomes)
}

build_paper_c_outcomes_ri_inputs <- function(outcomes_design, fits, settings) {
  outcomes <- paper_c_outcome_names()
  purrr::map(rlang::set_names(outcomes), \(outcome) {
    spec <- fits$specifications[[outcome]]
    prepare_randomization_inputs(
      outcomes_design$design,
      terms = spec$terms,
      bin_resolution = spec$bin_resolution,
      replications = settings$ri_replications,
      seed = settings$seed,
      strata_vars = settings$ri_strata,
      label = outcome,
      outcome = outcome
    )
  })
}

# The report the production run wrote: the randomization per outcome, the
# Holm and Westfall-Young adjustments of the alert-day p-values, and the
# checks.
build_paper_c_outcomes_report <- function(
  analysis_design,
  panel,
  outcomes_design,
  checks,
  fits,
  jackknife,
  ri_inputs,
  ri_blocks,
  settings
) {
  outcomes <- paper_c_outcome_names()
  specifications <- fits$specifications
  randomization <- purrr::map(outcomes, function(outcome) {
    result <- assemble_randomization_inference(
      ri_inputs[[outcome]],
      ri_blocks[[outcome]],
      settings$seed,
      settings$ri_replications
    )
    result$summary <- summarise_randomization_inference(
      result$replications,
      specifications[[outcome]]$coefficients
    )
    result
  }) |>
    stats::setNames(outcomes)

  observed_alert <- fits$coefficients |>
    dplyr::filter(contrast == "alert_day") |>
    dplyr::transmute(outcome, statistic = abs(statistic))
  permuted_alert <- purrr::imap_dfr(randomization, function(x, outcome) {
    x$replications |>
      dplyr::filter(contrast == "alert_day") |>
      dplyr::transmute(outcome = outcome, seed, statistic = abs(statistic))
  })
  ri_alert <- purrr::imap_dfr(randomization, function(x, outcome) {
    x$summary |>
      dplyr::filter(contrast == "alert_day") |>
      dplyr::transmute(outcome = outcome, p_raw = p_value)
  })
  holm <- holm_adjust(stats::setNames(ri_alert$p_raw, ri_alert$outcome))
  wy <- westfall_young_stepdown(observed_alert, permuted_alert)
  multiplicity <- holm |>
    dplyr::left_join(
      wy |> dplyr::select(outcome, p_raw_wy = p_raw, p_wy),
      by = dplyr::join_by(outcome),
      relationship = "one-to-one"
    ) |>
    dplyr::left_join(
      ri_alert |> dplyr::rename(p_raw_ri = p_raw),
      by = dplyr::join_by(outcome),
      relationship = "one-to-one"
    )
  ri_match_passed <- identical(multiplicity$p_raw_wy, multiplicity$p_raw_ri) &&
    identical(multiplicity$p_raw, multiplicity$p_raw_ri)
  if (!ri_match_passed) {
    rlang::abort("Westfall-Young raw p-values differ from outcome-wise RI.")
  }

  report <- list(
    design_summary = analysis_design$summary,
    join_checks = outcomes_design$join_checks,
    baseline_equivalence = checks$baseline_equivalence,
    outcome_equivalence = fits$outcome_equivalence,
    shared_permutation_checks = checks$shared_permutation_checks,
    shared_permutation_hashes = checks$shared_permutation_hashes,
    outcome_columns_hash = checks$outcome_columns_hash,
    permutation_label_checks = checks$permutation_label_checks,
    westfall_young_checks = checks$westfall_young_checks,
    finiteness_checks = checks$finiteness_checks,
    ri_multiplicity_match = tibble::tibble(passed = ri_match_passed),
    specifications = purrr::map(specifications, strip_fit),
    coefficients = fits$coefficients,
    jackknife = jackknife,
    randomization = randomization,
    randomization_smoke = NULL,
    multiplicity = multiplicity,
    multiplicity_smoke = NULL,
    settings = c(
      utils::modifyList(settings, list(workers = paper_c_workers())),
      list(timings = list())
    ),
    session = paper_c_session()
  )
  check_no_saitama(report, saitama_departments(panel), "The outcome report")
  report
}
