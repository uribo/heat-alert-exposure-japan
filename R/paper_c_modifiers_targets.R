# NOTE: readable analysis code from the private project. It needs the fire-department
# transport records, which are not distributed, and is not sourced by the public
# pipeline (_targets.R). Kept so that the published estimates can be audited.

# Result 2's two pre-specified effect modifiers (paper_strategy 7.14) as
# targets (Issue #50 (a)): the computation of src/paper_c_modifiers.R. The
# driver ran in two stages, a smoke run (B = 20) that fitted the thirteen
# specifications and the jackknife, and a production run (B = 999) that reused
# them after exact design-hash checks. Both stages are one pipeline here. The
# smoke draws are the first 20 seeds of the production grid, since the grid
# depends on the seed alone (randomization_seeds()), so they are cut from the
# production blocks instead of being drawn again.

# Settings in the order the production run stored them. `workers` is filled
# in when the report is assembled; `ri_p_floor` is added there too.
paper_c_modifiers_settings <- function(
  ri_replications = 999L,
  ri_chunk_size = 20L,
  smoke_replications = 20L
) {
  settings <- list(
    ri_replications = ri_replications,
    ri_chunk_size = ri_chunk_size,
    workers = NA_integer_,
    bin = "0.1",
    output = "notes/paper_c_modifiers_report.rds",
    reuse_report = "notes/paper_c_modifiers_report.rds",
    stage = "production"
  )
  settings$seed <- 20260914L
  settings$ri_strata <- c("area_primary", "year", "area_bin_t_1c")
  settings$income_fiscal_year <- 2020L
  settings$fixest_version <- as.character(utils::packageVersion("fixest"))
  attr(settings, "smoke_replications") <- smoke_replications
  settings
}

# The four specifications that carry the jackknife and the randomization.
paper_c_modifiers_primary_labels <- function() {
  c("high_primary", "high_stratified", "income_primary", "income_stratified")
}

# The thirteen specifications: the local-heat modifier on the six baseline
# samples and the income modifier on the same six plus the resident
# denominator.
paper_c_modifiers_specification_designs <- function(designs, panel, bin) {
  straddling_departments <- panel |>
    dplyr::filter(straddles) |>
    dplyr::distinct(department) |>
    dplyr::pull(department)
  income_resident_design <- function(design) {
    design |>
      dplyr::mutate(
        income_z = income_resident_z,
        alert_income = alert_income_resident
      )
  }
  make_modifier_spec <- function(
    design,
    modifier,
    bin,
    include_trial = TRUE,
    denominator = NA_character_,
    analysis
  ) {
    modifier_term <- switch(
      modifier,
      high = "alert_high",
      income = "alert_income",
      rlang::abort(glue::glue("Unknown modifier: {modifier}"))
    )
    list(
      design = design,
      terms = c(
        "alert_day",
        if (include_trial) "alert_trial",
        modifier_term
      ),
      bin = bin,
      modifier = modifier,
      denominator = denominator,
      analysis = analysis
    )
  }
  high_specs <- list(
    high_primary = make_modifier_spec(
      designs$full,
      "high",
      bin,
      analysis = "primary"
    ),
    high_stratified = make_modifier_spec(
      designs$stratified,
      "high",
      bin,
      analysis = "stratified"
    ),
    high_exclude_2020 = make_modifier_spec(
      designs$full |> dplyr::filter(year != 2020L),
      "high",
      bin,
      include_trial = FALSE,
      analysis = "exclude_2020"
    ),
    high_bins_1c = make_modifier_spec(
      designs$full,
      "high",
      "1",
      analysis = "bins_1c"
    ),
    high_exclude_straddling = make_modifier_spec(
      designs$full |>
        dplyr::filter(!department %in% straddling_departments),
      "high",
      bin,
      analysis = "exclude_straddling"
    ),
    high_exclude_high_leverage = make_modifier_spec(
      designs$full |> dplyr::filter(area_primary != "神奈川"),
      "high",
      bin,
      analysis = "exclude_high_leverage"
    )
  )
  income_specs <- list(
    income_primary = make_modifier_spec(
      designs$full,
      "income",
      bin,
      denominator = "taxpayer",
      analysis = "primary"
    ),
    income_stratified = make_modifier_spec(
      designs$stratified,
      "income",
      bin,
      denominator = "taxpayer",
      analysis = "stratified"
    ),
    income_exclude_2020 = make_modifier_spec(
      designs$full |> dplyr::filter(year != 2020L),
      "income",
      bin,
      include_trial = FALSE,
      denominator = "taxpayer",
      analysis = "exclude_2020"
    ),
    income_bins_1c = make_modifier_spec(
      designs$full,
      "income",
      "1",
      denominator = "taxpayer",
      analysis = "bins_1c"
    ),
    income_exclude_straddling = make_modifier_spec(
      designs$full |>
        dplyr::filter(!department %in% straddling_departments),
      "income",
      bin,
      denominator = "taxpayer",
      analysis = "exclude_straddling"
    ),
    income_exclude_high_leverage = make_modifier_spec(
      designs$full |> dplyr::filter(area_primary != "神奈川"),
      "income",
      bin,
      denominator = "taxpayer",
      analysis = "exclude_high_leverage"
    ),
    income_resident = make_modifier_spec(
      income_resident_design(designs$full),
      "income",
      bin,
      denominator = "resident",
      analysis = "resident_denominator"
    )
  )
  c(high_specs, income_specs)
}

# The stage-1 fits: the thirteen specifications, their metadata, and the
# baseline equivalence (the primary baseline refitted on the income design
# must reproduce the baseline fits exactly). The fit objects are dropped; the
# kept-row index stays for the jackknife.
build_paper_c_modifiers_fits <- function(
  income_design,
  panel,
  result2_fits,
  settings
) {
  designs <- income_design$designs
  specification_designs <- paper_c_modifiers_specification_designs(
    designs,
    panel,
    settings$bin
  )
  specification_metadata <- purrr::imap(specification_designs, \(spec, label) {
    tibble::tibble(
      specification = label,
      modifier = spec$modifier,
      analysis = spec$analysis,
      denominator = spec$denominator,
      bin_resolution = spec$bin,
      rows = nrow(spec$design)
    )
  }) |>
    purrr::list_rbind()
  specifications <- purrr::imap(specification_designs, \(spec, label) {
    fit <- run_paper_c_result2_spec(
      spec$design,
      label = label,
      terms = spec$terms,
      bin_resolution = spec$bin,
      bootstrap_replications = 0L,
      seed = settings$seed
    )
    fit$fit <- NULL
    fit
  })

  baseline_refit <- run_paper_c_result2_spec(
    designs$full,
    label = "primary",
    terms = c("alert_day", "alert_trial"),
    bin_resolution = settings$bin,
    bootstrap_replications = 0L,
    seed = settings$seed
  )
  baseline_stored <- result2_fits$coefficients |>
    dplyr::filter(specification == "primary") |>
    dplyr::select(contrast, estimate_stored = estimate, se_stored = se) |>
    dplyr::arrange(contrast)
  baseline_current <- baseline_refit$coefficients |>
    dplyr::select(contrast, estimate_current = estimate, se_current = se) |>
    dplyr::arrange(contrast)
  baseline_equivalence <- baseline_stored |>
    dplyr::inner_join(
      baseline_current,
      by = dplyr::join_by(contrast),
      relationship = "one-to-one"
    ) |>
    dplyr::mutate(
      estimate_difference = estimate_current - estimate_stored,
      se_difference = se_current - se_stored
    )
  if (
    nrow(baseline_equivalence) != nrow(baseline_stored) ||
      nrow(baseline_equivalence) != nrow(baseline_current) ||
      any(baseline_equivalence$estimate_difference != 0) ||
      any(baseline_equivalence$se_difference != 0)
  ) {
    rlang::abort("The baseline primary coefficients or CR1 SEs changed.")
  }
  list(
    specification_metadata = specification_metadata,
    specifications = specifications,
    coefficients = purrr::map(specifications, "coefficients") |>
      purrr::list_rbind(),
    baseline_equivalence = baseline_equivalence,
    straddling_departments = panel |>
      dplyr::filter(straddles) |>
      dplyr::distinct(department) |>
      dplyr::pull(department) |>
      length()
  )
}

# Leave-one-area-out jackknife of the four primary specifications.
build_paper_c_modifiers_jackknife <- function(
  income_design,
  panel,
  fits,
  settings
) {
  specification_designs <- paper_c_modifiers_specification_designs(
    income_design$designs,
    panel,
    settings$bin
  )
  labels <- paper_c_modifiers_primary_labels()
  purrr::map(rlang::set_names(labels), \(label) {
    spec <- fits$specifications[[label]]
    design_kept <- specification_designs[[label]]$design[
      spec$design_kept_rows,
    ]
    result <- paper_c_jackknife_vcov(
      design_kept,
      spec$terms,
      spec$bin_resolution,
      workers = paper_c_workers()
    )
    result$summary <- summarise_jackknife(
      result,
      spec$estimates,
      paper_c_result2_contrasts(spec$terms)
    )
    result
  })
}

build_paper_c_modifiers_ri_inputs <- function(
  income_design,
  panel,
  fits,
  settings
) {
  specification_designs <- paper_c_modifiers_specification_designs(
    income_design$designs,
    panel,
    settings$bin
  )
  labels <- paper_c_modifiers_primary_labels()
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

# The acceptance check of the permuted modifier design, on the smoke grid.
build_paper_c_modifiers_permutation_acceptance <- function(
  income_design,
  settings
) {
  check_permuted_modifier_design(
    income_design$designs$full,
    replications = min(settings$ri_replications, 20L),
    seed = settings$seed,
    strata_vars = settings$ri_strata
  )
}

# The report the production run wrote.
build_paper_c_modifiers_report <- function(
  income_design,
  panel,
  fits,
  jackknife,
  permutation_acceptance,
  ri_inputs,
  ri_blocks,
  settings
) {
  specifications <- fits$specifications
  smoke_replications <- attr(settings, "smoke_replications")
  attr(settings, "smoke_replications") <- NULL
  smoke_seeds <- randomization_seeds(settings$seed, smoke_replications)

  randomization <- purrr::imap(ri_inputs, \(inputs, label) {
    draws <- assemble_randomization_inference(
      inputs,
      ri_blocks[[label]],
      settings$seed,
      settings$ri_replications
    )
    draws$summary <- summarise_randomization_inference(
      draws$replications,
      specifications[[label]]$coefficients
    )
    draws
  })
  # The smoke run's draws: the first seeds of the same grid.
  randomization_smoke <- purrr::imap(randomization, \(draws, label) {
    smoke <- draws[c("strata", "strata_vars", "frozen_years")]
    smoke$replications <- draws$replications |>
      dplyr::filter(seed %in% smoke_seeds)
    smoke$summary <- summarise_randomization_inference(
      smoke$replications,
      specifications[[label]]$coefficients
    )
    smoke
  })

  jackknife_summary <- purrr::imap(jackknife, \(x, label) {
    x$summary |>
      dplyr::mutate(specification = label, .before = 1)
  }) |>
    purrr::list_rbind()
  randomization_summary <- purrr::imap(randomization, \(result, label) {
    result$summary |>
      dplyr::transmute(specification = label, contrast, p_value_ri = p_value)
  }) |>
    purrr::list_rbind()
  coefficients <- fits$coefficients |>
    dplyr::left_join(
      jackknife_summary,
      by = dplyr::join_by(specification, contrast)
    ) |>
    dplyr::left_join(
      randomization_summary,
      by = dplyr::join_by(specification, contrast),
      relationship = "one-to-one"
    )
  settings$ri_p_floor <- 1 / (settings$ri_replications + 1)

  modifier_results <- paper_c_modifier_results(
    coefficients,
    fits$specification_metadata
  )

  report <- list(
    design_summary = income_design$summary,
    specification_metadata = fits$specification_metadata,
    specifications = purrr::map(specifications, strip_fit),
    coefficients = coefficients,
    modifier_results = modifier_results,
    jackknife = jackknife,
    randomization = randomization,
    randomization_smoke = randomization_smoke,
    permutation_acceptance = permutation_acceptance,
    baseline_equivalence = fits$baseline_equivalence,
    straddling_departments = fits$straddling_departments,
    settings = c(
      utils::modifyList(settings, list(workers = paper_c_workers())),
      list(timings = list())
    ),
    session = paper_c_session()
  )
  check_no_saitama(report, saitama_departments(panel), "The modifier report")
  report
}

# The modifier's own contrast in each specification, with its metadata.
paper_c_modifier_results <- function(coefficients, specification_metadata) {
  modifier_contrasts <- c(high = "alert_high", income = "alert_income")
  coefficients |>
    dplyr::left_join(
      specification_metadata |>
        dplyr::select(specification, modifier, analysis, denominator),
      by = dplyr::join_by(specification),
      relationship = "many-to-one"
    ) |>
    dplyr::filter(contrast == modifier_contrasts[modifier])
}

# Power of the continuous income modifier (paper_strategy 7.14) ----------------

# Settings of src/paper_c_income_power.R in its order. `workers` is filled in
# when the report is assembled; `theta` becomes the anchors' estimates there.
paper_c_income_power_settings <- function(replications = 500L) {
  settings <- list(
    replications = replications,
    workers = NA_integer_,
    bin = "0.1",
    theta = "moment",
    relative_effects = c(0, 0.05, 0.10, 0.15, 0.20, 0.30),
    output = "notes/paper_c_income_power_report.rds",
    modifier_report = "notes/paper_c_modifiers_report.rds",
    resume = FALSE
  )
  settings$seed <- 20260907L
  settings$effect_term <- "alert_income"
  settings$income_fiscal_year <- 2020L
  settings$dgp <- "nb"
  settings$ssc <- "fixest::ssc(adj = TRUE, cluster.adj = TRUE); p from t with G - 1 df"
  settings$fixef_tol <- formals(fixest::fepois)$fixef.tol
  settings$fixest_version <- as.character(utils::packageVersion("fixest"))
  settings
}

# The anchors of the simulation: the null and full fits on both samples,
# required to reproduce the reviewed CR1 estimates of the income modifier.
# Returns the anchors, the reviewed estimates, and the settings completed with
# the design and anchor fingerprints a stored chunk must share.
build_paper_c_income_power_anchors <- function(
  income_design,
  modifiers_fits,
  settings
) {
  designs <- income_design$designs
  settings$design_hash <- purrr::map_chr(designs, function(design) {
    paper_c_design_hash(design)$hash
  })
  reviewed_inference <- paper_c_modifier_results(
    modifiers_fits$coefficients,
    modifiers_fits$specification_metadata
  ) |>
    dplyr::filter(
      modifier == "income",
      analysis %in% c("primary", "stratified")
    ) |>
    dplyr::transmute(
      sample = dplyr::case_match(
        analysis,
        "primary" ~ "full",
        "stratified" ~ "stratified"
      ),
      term = contrast,
      estimate,
      se,
      statistic,
      clusters,
      p_value = p_value_cr1
    )
  if (
    nrow(reviewed_inference) != 2L ||
      !setequal(reviewed_inference$sample, names(designs))
  ) {
    rlang::abort("The modifier fits lack the two reviewed income estimates.")
  }
  anchors <- purrr::imap(designs, function(sample_design, label) {
    real_inference <- reviewed_inference |>
      dplyr::filter(sample == label) |>
      dplyr::select(-sample)
    fit_paper_c_power_anchors(
      sample_design,
      bin_resolution = settings$bin,
      theta_method = settings$theta,
      label = label,
      effect_term = settings$effect_term,
      terms = c("alert_day", "alert_trial", settings$effect_term),
      real_inference = real_inference
    )
  })
  # What a chunk must have been drawn from: the null fit's fitted means and
  # coefficients and the modifier column. Settings, seeds and theta alone do
  # not see a corrected income table, since the null model carries no income
  # term.
  settings$anchor_hash <- purrr::map(anchors, function(anchor) {
    rlang::hash(list(
      mu = anchor$mu,
      theta = anchor$theta$value,
      null_coefficients = anchor$null_coefficients,
      effect = anchor$kept_design[[anchor$effect_term]]
    ))
  }) |>
    unlist()
  theta_values <- purrr::map(anchors, function(anchor) anchor$theta$value)
  list(
    anchors = anchors,
    reviewed_inference = reviewed_inference,
    settings = settings,
    partial_settings = paper_c_power_partial_settings(settings, theta_values)
  )
}

# The task grid split into its chunks (one per sample x relative effect), one
# dynamic branch per chunk.
build_paper_c_income_power_chunks <- function(income_design, settings) {
  tasks <- build_paper_c_power_tasks(
    samples = names(income_design$designs),
    relative_effects = settings$relative_effects,
    replications = settings$replications,
    seed = settings$seed
  )
  unname(split(tasks, tasks$chunk))
}

# One branch: the fits of one chunk, with the settings they were drawn under.
run_paper_c_income_power_branch <- function(anchors, chunk, workers) {
  list(
    settings = anchors$partial_settings,
    results = run_paper_c_power_chunk(anchors$anchors, chunk, workers)
  )
}

# The report src/paper_c_income_power.R wrote, after checking that every
# chunk belongs to the current anchors and that the chunks cover the grid.
build_paper_c_income_power_report <- function(
  income_design,
  panel,
  anchors,
  chunks,
  simulation
) {
  stale <- !purrr::map_lgl(simulation, \(branch) {
    identical(branch$settings, anchors$partial_settings)
  })
  if (any(stale)) {
    rlang::abort(glue::glue(
      "{sum(stale)} of {length(simulation)} income-power chunks were computed ",
      "from different anchors or settings (the chunks are frozen by ",
      "paper_c_ri_cue()). Rebuild them deliberately with ",
      "PAPER_C_RI_REBUILD=true."
    ))
  }
  replication_results <- dplyr::bind_rows(purrr::map(simulation, "results")) |>
    dplyr::arrange(sample, dgp, relative_effect, replication)
  # Every task of the current grid exactly once, drawn from the stream the
  # grid assigns it. A row count alone passes a duplicated task that hides a
  # missing one.
  keys <- c("sample", "dgp", "relative_effect", "replication")
  grid <- dplyr::bind_rows(chunks) |>
    dplyr::mutate(current_stream_key = purrr::map_chr(stream, stream_key)) |>
    dplyr::select(dplyr::all_of(keys), current_stream_key)
  task_keys <- dplyr::select(replication_results, dplyr::all_of(keys))
  duplicated_tasks <- nrow(task_keys) - nrow(dplyr::distinct(task_keys))
  matched <- replication_results |>
    dplyr::select(dplyr::all_of(keys), stream_key) |>
    dplyr::left_join(grid, by = keys)
  unknown <- sum(is.na(matched$current_stream_key))
  missing_tasks <- nrow(dplyr::anti_join(grid, matched, by = keys))
  # A missing stored key is a mismatch too; rows outside the grid are already
  # counted in `unknown`.
  in_grid <- !is.na(matched$current_stream_key)
  wrong_stream <- sum(
    in_grid &
      (is.na(matched$stream_key) |
        matched$stream_key != matched$current_stream_key)
  )
  if (duplicated_tasks + unknown + missing_tasks + wrong_stream > 0L) {
    rlang::abort(glue::glue(
      "The income-power chunks do not cover the task grid of {nrow(grid)} ",
      "tasks exactly once: {duplicated_tasks} duplicated, {unknown} not in ",
      "the grid, {missing_tasks} missing, {wrong_stream} with a missing or ",
      "different stream key."
    ))
  }
  settings <- anchors$settings
  anchor_summaries <- purrr::map(anchors$anchors, function(anchor) {
    anchor[c(
      "label",
      "effect_term",
      "bin_resolution",
      "n_clusters",
      "theta",
      "null_coefficients",
      "full_coefficients",
      "real_modifier",
      "real_modifier_refit",
      "null_removals",
      "full_removals",
      "dropped_effect_rows",
      "seconds_null_fit",
      "seconds_full_fit"
    )]
  })
  power <- summarise_paper_c_power(replication_results)
  settings$workers <- paper_c_workers()
  settings$theta <- purrr::map(anchors$anchors, "theta")
  settings$timings <- list()
  report <- list(
    design_summary = income_design$summary,
    anchor_summaries = anchor_summaries,
    reviewed_inference = anchors$reviewed_inference,
    power_analytic = analytic_paper_c_power(
      anchors$anchors,
      settings$relative_effects
    ),
    power = power,
    replication_results = replication_results,
    settings = settings,
    session = paper_c_session(),
    se_comparison = compare_paper_c_se(anchors$anchors, power)
  )
  check_no_saitama(
    report,
    saitama_departments(panel),
    "The income-power report"
  )
  report
}
