# NOTE: readable analysis code from the private project. It needs the fire-department
# transport records, which are not distributed, and is not sourced by the public
# pipeline (_targets.R). Kept so that the published estimates can be audited.

# Result 2's baseline for paper C (paper_strategy 7.2): the alert main effect
# with no modifier, on the model whose power was computed in 7.12.
#
#   transports_analysis ~ alert_day + alert_trial |
#     department + date + area_bin_t + area_bin_lag + local_bin
#
# 7.12 measured the clustered CR1 t test's size at 8.4 % (full) and 10.4 %
# (stratified) against a nominal 5 %, so the CR1 p-value is reported but is
# not the test. Three readings are computed:
#
#   * Randomization inference, the primary reading: the area-day alert labels
#     are permuted inside area x year x area-maximum WBGT (1 degC) strata and
#     the model refitted. It approximates the assignment mechanism (the alert
#     is issued on the forecast, and the strata condition on realized heat),
#     and being a
#     permutation test of the sharp null it does not rely on the clustered
#     variance estimator being unbiased.
#   * A leave-one-area-out jackknife (CV3J) covariance, which is the remedy for
#     exactly the downward bias in CR1 that 7.12 measured, and which carries
#     the interval randomization inference cannot give for a count outcome.
#   * A score wild cluster bootstrap (Kline & Santos 2012), kept as a
#     heuristic sensitivity rather than as 7.2's wild cluster bootstrap.
#     With Rademacher weights its studentized denominator is constant across
#     replications (w^2 == 1), so it reduces to a normal approximation and
#     cannot reproduce the finite-sample inflation the bootstrap is wanted
#     for. It is off by default.
#
# Every function is pure. The targets of R/paper_c_result2_targets.R call
# them in order (Issue #50).

# Model formula from an explicit term vector, so a specification that has no
# 2020 trial rows can drop `alert_trial` instead of carrying an all-zero
# column that fixest would silently collinearly remove.
# `drop_local = TRUE` removes the department's local WBGT bins and keeps the
# area-level bins: the specification that isolates what local conditioning
# adds to the alert coefficient (Issue #43, item 1). The default reproduces
# the primary model byte for byte.
paper_c_result2_formula <- function(
  terms = c("alert_day", "alert_trial"),
  bin_resolution = c("0.1", "1"),
  drop_local = FALSE
) {
  bin_resolution <- match.arg(bin_resolution)
  suffix <- if (bin_resolution == "1") "_1c" else ""
  # An empty term vector is the fully restricted model: fixed effects only.
  rhs <- if (length(terms) == 0L) "1" else paste(terms, collapse = " + ")
  local <- if (drop_local) "" else paste0(" + local_bin", suffix)
  stats::as.formula(paste0(
    "transports_analysis ~ ",
    rhs,
    " | department + date + area_bin_t",
    suffix,
    " + area_bin_lag",
    suffix,
    local
  ))
}

paper_c_result2_fixef <- function(
  bin_resolution = c("0.1", "1"),
  drop_local = FALSE
) {
  bin_resolution <- match.arg(bin_resolution)
  suffix <- if (bin_resolution == "1") "_1c" else ""
  bins <- c("area_bin_t", "area_bin_lag", if (!drop_local) "local_bin")
  c("department", "date", paste0(bins, suffix))
}

# Cluster score contributions and the Hessian of the concentrated Poisson
# likelihood, at whatever coefficient vector produced `mu`. With the fixed
# effects profiled out, the score of observation i for regressor j is
# (y_i - mu_i) * x~_ij, where x~ is x weighted-demeaned by the fixed effects
# with weights mu, and the Hessian is x~' diag(mu) x~. Computing them here
# rather than reading `fit$scores` is what makes the null-imposed bootstrap
# possible: the restricted fit does not carry the tested regressor.
paper_c_score_pieces <- function(
  design,
  terms,
  mu,
  fixef_vars,
  cluster,
  tol = 1e-8
) {
  if (nrow(design) != length(mu)) {
    rlang::abort("`mu` must have one value per design row.")
  }
  x <- as.matrix(design[, terms, drop = FALSE])
  x_demeaned <- fixest::demean(
    x,
    design[, fixef_vars, drop = FALSE],
    weights = mu,
    tol = tol
  )
  residual <- design$transports_analysis - mu
  scores <- x_demeaned * residual
  colnames(scores) <- terms
  cluster_scores <- rowsum(scores, cluster, reorder = TRUE)
  list(
    terms = terms,
    scores = scores,
    cluster_scores = cluster_scores,
    hessian = crossprod(x_demeaned * sqrt(mu)),
    n_clusters = nrow(cluster_scores)
  )
}

# CR1 covariance from the pieces, times a scalar small-sample factor. The
# factor is not derived here: it is read off the fitted model once (see
# paper_c_ssc_scalar()) so that the reported standard errors are fixest's own
# and the bootstrap statistics use exactly the same estimator.
paper_c_cluster_vcov <- function(pieces, ssc_scalar = 1) {
  bread <- solve(pieces$hessian)
  vcov <- bread %*% crossprod(pieces$cluster_scores) %*% bread * ssc_scalar
  dimnames(vcov) <- list(pieces$terms, pieces$terms)
  vcov
}

# The scalar that turns the unadjusted CR1 covariance into fixest's. It must
# be one number for every coefficient; if it is not, the two estimators differ
# in shape and the bootstrap cannot inherit fixest's standard errors.
paper_c_ssc_scalar <- function(pieces, fit, tolerance = 1e-4) {
  fixest_se <- fixest::se(
    fit,
    cluster = ~area_primary,
    ssc = paper_c_power_ssc()
  )[pieces$terms]
  manual_se <- sqrt(diag(paper_c_cluster_vcov(pieces)))
  ratio <- (fixest_se / manual_se)^2
  # A relative tolerance: the weighted demeaning behind `manual_se` converges
  # to a finite precision, so the ratios agree to about 1e-6 rather than
  # exactly.
  if (diff(range(ratio)) / mean(ratio) > tolerance) {
    rlang::abort(glue::glue(
      "The small-sample factor is not a scalar: ",
      "{paste(signif(ratio, 8), collapse = ', ')}."
    ))
  }
  mean(ratio)
}

# Guard: the manually computed pieces must reproduce what fixest carries for
# the same (unrestricted) fit. Checked on every specification, because the
# restricted pieces have nothing to be checked against.
check_score_pieces <- function(
  pieces,
  fit,
  score_tolerance = 1e-4,
  hessian_tolerance = 1e-6
) {
  score_gap <- max(abs(pieces$scores - fit$scores)) / max(abs(fit$scores))
  hessian_gap <- max(abs(pieces$hessian - fit$hessian)) / max(abs(fit$hessian))
  if (score_gap > score_tolerance || hessian_gap > hessian_tolerance) {
    rlang::abort(glue::glue(
      "Manual scores/Hessian disagree with fixest: relative gaps ",
      "{signif(score_gap, 3)} / {signif(hessian_gap, 3)}."
    ))
  }
  tibble::tibble(score_gap = score_gap, hessian_gap = hessian_gap)
}

# Rademacher weights, one per cluster per replication.
rademacher_weights <- function(n_clusters, replications) {
  matrix(
    sample(c(-1, 1), n_clusters * replications, replace = TRUE),
    nrow = n_clusters,
    ncol = replications
  )
}

# Studentized statistic for a linear combination `contrast` of the terms,
# one value per column of `weights`. With the unrestricted pieces and unit
# weights this is identically zero (the scores sum to zero by the first-order
# conditions), so the observed statistic always comes from the fit itself.
score_bootstrap_statistics <- function(
  pieces,
  contrast,
  ssc_scalar,
  weights
) {
  bread <- solve(pieces$hessian)
  cluster_scores <- pieces$cluster_scores
  direction <- as.numeric(crossprod(bread, contrast))
  # Cluster contribution to the estimate's numerator, one value per cluster.
  numerator_parts <- as.numeric(cluster_scores %*% direction)
  estimates <- as.numeric(crossprod(numerator_parts, weights))
  variances <- as.numeric(crossprod(numerator_parts^2, weights^2)) * ssc_scalar
  estimates / sqrt(variances)
}

# Score wild cluster bootstrap for one hypothesis. `pieces_null` carries the
# scores of the restricted fit (the null imposed); `statistic` is the
# unrestricted studentized statistic the bootstrap distribution is compared
# with. Two-sided symmetric p-value with the usual (b + 1) / (B + 1) rule.
score_wild_cluster_p <- function(
  pieces_null,
  contrast,
  ssc_scalar,
  statistic,
  replications,
  seed
) {
  withr::with_seed(seed, {
    weights <- rademacher_weights(pieces_null$n_clusters, replications)
  })
  bootstrap <- score_bootstrap_statistics(
    pieces_null,
    contrast,
    ssc_scalar,
    weights
  )
  finite <- is.finite(bootstrap)
  tibble::tibble(
    p_value = (sum(abs(bootstrap[finite]) >= abs(statistic)) + 1) /
      (sum(finite) + 1),
    replications = sum(finite),
    replications_dropped = sum(!finite),
    # The bootstrap distribution's spread against the t(G - 1) the clustered
    # test assumes: a value below one says the CR1 critical value is
    # conservative in this sample, above one that it is liberal.
    bootstrap_sd = stats::sd(bootstrap[finite])
  )
}

# Symmetric studentized interval from the unrestricted score bootstrap: the
# 1 - alpha quantile of |t*| replaces the t(G - 1) critical value.
score_wild_cluster_interval <- function(
  pieces,
  contrast,
  ssc_scalar,
  estimate,
  se,
  replications,
  seed,
  level = 0.95
) {
  withr::with_seed(seed, {
    weights <- rademacher_weights(pieces$n_clusters, replications)
  })
  bootstrap <- score_bootstrap_statistics(pieces, contrast, ssc_scalar, weights)
  bootstrap <- bootstrap[is.finite(bootstrap)]
  critical <- stats::quantile(abs(bootstrap), probs = level, names = FALSE)
  tibble::tibble(
    critical_value = critical,
    conf_low = estimate - critical * se,
    conf_high = estimate + critical * se,
    replications = length(bootstrap)
  )
}

# Randomization inference -----------------------------------------------------

# Strata for the permutation. The alert is an area-day treatment issued on the
# forecast, so labels are exchanged only between days of the same area, the
# same year, and the same realized area-maximum WBGT degree: permuting across
# temperatures would compare an alert with days that could never have carried
# one. Strata whose days are all alerted or all unalerted contribute no
# variation and are left as they are.
build_alert_permutation_strata <- function(
  design,
  strata_vars = c("area_primary", "year", "area_bin_t_1c"),
  frozen_years = integer(0)
) {
  design |>
    dplyr::distinct(dplyr::pick(dplyr::all_of(c(
      strata_vars,
      "year",
      "date",
      "alert_day"
    )))) |>
    dplyr::mutate(
      stratum = dplyr::cur_group_id(),
      .by = dplyr::all_of(strata_vars)
    ) |>
    # A frozen year keeps its actual assignment. The global sharp null makes
    # every coefficient zero at once, so a p-value from an unrestricted
    # permutation tests "the alert did nothing anywhere" and not "this one
    # coefficient is zero". Freezing 2020 leaves the trial rows at their
    # observed assignment, so the trial coefficient is still estimated from
    # real data while the nationwide coefficient is tested on its own.
    dplyr::mutate(frozen = year %in% frozen_years)
}

# Count what the permutation can and cannot move, so the note can report the
# fraction of treated area-days that carries randomization variation.
summarise_permutation_strata <- function(area_days) {
  frozen <- area_days |> dplyr::filter(frozen)
  area_days <- area_days |> dplyr::filter(!frozen)
  by_stratum <- area_days |>
    dplyr::summarise(
      area_days = dplyr::n(),
      alerted = sum(alert_day),
      .by = stratum
    )
  mixed <- by_stratum |> dplyr::filter(alerted > 0, alerted < area_days)
  tibble::tibble(
    strata = nrow(by_stratum),
    mixed_strata = nrow(mixed),
    area_days = sum(by_stratum$area_days),
    alerted_area_days = sum(by_stratum$alerted),
    alerted_in_mixed_strata = sum(mixed$alerted),
    median_mixed_stratum_size = stats::median(mixed$area_days),
    frozen_area_days = nrow(frozen),
    frozen_alerted_area_days = sum(frozen$alert_day)
  )
}

# One permuted assignment: shuffle the alert labels inside every stratum,
# then rebuild the design's treatment columns from the area-day labels.
permute_alert_design <- function(
  design,
  area_days,
  seed,
  terms = names(design)
) {
  registry <- paper_c_alert_term_registry()
  allowed_alert_terms <- c("alert_day", names(registry))
  unknown_alert_terms <- setdiff(
    unique(c(
      grep("^alert_", terms, value = TRUE),
      grep("^alert_", names(design), value = TRUE)
    )),
    allowed_alert_terms
  )
  if (length(unknown_alert_terms) > 0L) {
    rlang::abort(paste(
      "Cannot permute unknown alert term(s):",
      paste(unknown_alert_terms, collapse = ", ")
    ))
  }
  rebuild_terms <- unique(c(
    "alert_trial",
    intersect(terms, names(registry))
  ))
  purrr::walk(rebuild_terms, function(term) {
    input <- registry[[term]]$input
    if (!input %in% names(design)) {
      rlang::abort(glue::glue(
        "`{input}` is required to rebuild `{term}`."
      ))
    }
  })
  permuted <- withr::with_seed(seed, {
    area_days |>
      dplyr::mutate(
        alert_day = if (dplyr::first(frozen)) {
          alert_day
        } else {
          alert_day[sample.int(dplyr::n())]
        },
        .by = stratum
      ) |>
      dplyr::select(area_primary, date, alert_permuted = alert_day)
  })
  # Every registry column is dropped here and only those in `terms` are put
  # back. A term the caller left out therefore disappears rather than keeping
  # its pre-permutation values, so a formula that names one fails at the fit
  # instead of reading a column built from the observed assignment.
  result <- design |>
    dplyr::select(
      !dplyr::any_of(c(
        "alert_day",
        names(registry)
      ))
    ) |>
    dplyr::left_join(
      permuted,
      by = dplyr::join_by(area_primary, date),
      relationship = "many-to-one"
    ) |>
    dplyr::mutate(alert_day = alert_permuted, alert_permuted = NULL)
  rebuild_alert_terms(result, rebuild_terms)
}

# Alert runs against the permutation (Issue #43, item 5). Alerts come in runs
# of consecutive days, and the shuffle of permute_alert_design() exchanges
# single days, so it does not preserve them. This counts (i) the observed
# runs, (ii) the share of consecutive days that fall in the same stratum, and
# (iii) the runs the shuffle produces over the RI's own seed grid. The shuffle
# is the expression of permute_alert_design() applied to the area-day table,
# and the first seed is checked against permute_alert_design() itself, so the
# permuted runs are those of the production draws. Nothing is fitted: an odds
# ratio of an alert on the previous day's alert within strata was tried and
# dropped, because with the lagged outcome as regressor neither the
# fixed-effects nor the conditional logit removes the stratum effects (both
# fell to about 0.5 for a true odds ratio of 1 in a simulation; PR #58).
summarise_alert_runs <- function(inputs) {
  settings <- inputs$settings
  area_days <- inputs$area_days
  seeds <- randomization_seeds(settings$seed, settings$replications)
  shuffle <- function(seed) {
    withr::with_seed(seed, {
      area_days |>
        dplyr::mutate(
          alert_day = if (dplyr::first(frozen)) {
            alert_day
          } else {
            alert_day[sample.int(dplyr::n())]
          },
          .by = stratum
        )
    })
  }
  reference <- permute_alert_design(
    inputs$design,
    area_days,
    seeds[[1]],
    terms = settings$terms
  ) |>
    dplyr::distinct(area_primary, date, alert_reference = alert_day)
  first_draw <- shuffle(seeds[[1]]) |>
    dplyr::inner_join(reference, by = dplyr::join_by(area_primary, date))
  if (
    nrow(first_draw) != nrow(area_days) ||
      !identical(
        as.numeric(first_draw$alert_day),
        as.numeric(first_draw$alert_reference)
      )
  ) {
    rlang::abort(glue::glue(
      "The run summary's shuffle does not reproduce permute_alert_design() ",
      "for {settings$label} at seed {seeds[[1]]}."
    ))
  }

  # Order once; a run continues only within an area-year and across
  # consecutive calendar days. Frozen years keep their assignment and are
  # left out of every count.
  ordering <- order(area_days$area_primary, area_days$date)
  sorted <- area_days[ordering, ]
  continues <- c(
    FALSE,
    sorted$area_primary[-1] == sorted$area_primary[-nrow(sorted)] &
      sorted$year[-1] == sorted$year[-nrow(sorted)] &
      as.integer(diff(sorted$date)) == 1L
  )
  active <- !sorted$frozen
  run_lengths <- function(alert) {
    alert <- alert == 1 & active
    starts <- alert & !(continues & c(FALSE, alert[-length(alert)]))
    tabulate(cumsum(starts)[alert])
  }
  observed <- run_lengths(sorted$alert_day)

  pairs <- tibble::tibble(
    area_primary = sorted$area_primary,
    stratum = sorted$stratum,
    alert_day = sorted$alert_day,
    previous_alert = c(NA, sorted$alert_day[-nrow(sorted)]),
    previous_stratum = c(NA, sorted$stratum[-nrow(sorted)])
  )[continues & active & c(FALSE, active[-length(active)]), ]
  alerted_pairs <- dplyr::filter(pairs, alert_day == 1, previous_alert == 1)

  permuted <- purrr::map(seeds, \(seed) {
    lengths <- run_lengths(shuffle(seed)$alert_day[ordering])
    c(
      runs = length(lengths),
      mean_length = mean(lengths),
      max_length = max(lengths)
    )
  }) |>
    do.call(what = rbind)

  tibble::tibble(
    area_days = sum(active),
    alerted_area_days = sum(observed),
    runs = length(observed),
    runs_length_1 = sum(observed == 1L),
    runs_length_2 = sum(observed == 2L),
    runs_length_3 = sum(observed == 3L),
    runs_length_4 = sum(observed == 4L),
    runs_length_5 = sum(observed == 5L),
    runs_length_6_plus = sum(observed >= 6L),
    share_alerted_in_runs_2_plus = sum(observed[observed >= 2L]) /
      sum(observed),
    mean_run_length = mean(observed),
    max_run_length = max(observed),
    consecutive_pairs = nrow(pairs),
    share_pairs_same_stratum = mean(pairs$stratum == pairs$previous_stratum),
    alert_pairs = nrow(alerted_pairs),
    share_alert_pairs_same_stratum = mean(
      alerted_pairs$stratum == alerted_pairs$previous_stratum
    ),
    permutations = length(seeds),
    permuted_runs_mean = mean(permuted[, "runs"]),
    permuted_runs_min = min(permuted[, "runs"]),
    permuted_runs_max = max(permuted[, "runs"]),
    permuted_mean_length_mean = mean(permuted[, "mean_length"]),
    permuted_max_length_mean = mean(permuted[, "max_length"]),
    permuted_runs_at_most_observed = sum(permuted[, "runs"] <= length(observed))
  )
}

# Acceptance check for the modifier-aware permutation. Both interaction
# columns must be exact row by row, and at least one seed in the smoke grid
# must move each interaction away from its observed assignment.
check_permuted_modifier_design <- function(
  design,
  replications,
  seed,
  strata_vars = c("area_primary", "year", "area_bin_t_1c")
) {
  registry <- paper_c_alert_term_registry()
  checked_terms <- c("alert_trial", "alert_high", "alert_income")
  required <- c(
    "alert_day",
    checked_terms,
    # Unnamed: a named vector inside `dplyr::all_of()` renames the column.
    unname(purrr::map_chr(registry[checked_terms], "input"))
  )
  missing <- setdiff(required, names(design))
  if (length(missing) > 0L) {
    rlang::abort(paste(
      "Cannot check modifier rebuilding; the design lacks:",
      paste(missing, collapse = ", ")
    ))
  }
  design <- design |>
    dplyr::select(dplyr::all_of(unique(c(
      "area_primary",
      "date",
      "year",
      required,
      strata_vars
    ))))
  area_days <- build_alert_permutation_strata(design, strata_vars)
  seeds <- randomization_seeds(seed, replications)
  observed_high <- sum(design$alert_high)
  checks <- purrr::map(seeds, function(replication_seed) {
    permuted <- permute_alert_design(
      design,
      area_days,
      replication_seed,
      terms = c("alert_day", "alert_trial", "alert_high", "alert_income")
    )
    high_gap <- max(abs(
      permuted$alert_high -
        as.numeric(permuted$alert_day == 1 & permuted$local_high)
    ))
    income_gap <- max(abs(
      permuted$alert_income - permuted$alert_day * permuted$income_z
    ))
    tibble::tibble(
      seed = replication_seed,
      high_row_gap = high_gap,
      income_row_gap = income_gap,
      alert_high_observed = observed_high,
      alert_high_permuted = sum(permuted$alert_high),
      alert_high_sum_changed = alert_high_permuted != alert_high_observed,
      alert_high_rows_changed = sum(permuted$alert_high != design$alert_high),
      alert_income_rows_changed = sum(
        permuted$alert_income != design$alert_income
      )
    )
  }) |>
    purrr::list_rbind()
  if (
    any(checks$high_row_gap != 0) ||
      any(checks$income_row_gap != 0) ||
      !any(checks$alert_high_rows_changed > 0) ||
      !any(checks$alert_income_rows_changed > 0)
  ) {
    rlang::abort(
      "Modifier permutation failed its row-wise rebuild or movement check."
    )
  }
  # The witness is the seed that moved the most alert_high rows; a row-wise
  # move is what the abort above requires, whereas the aggregate count can
  # stay fixed when the permuted alerts land on as many local-high days.
  witness <- checks |>
    dplyr::slice_max(alert_high_rows_changed, n = 1L, with_ties = FALSE)
  message(
    "Modifier permutation acceptance: row gaps high=",
    max(checks$high_row_gap),
    ", income=",
    max(checks$income_row_gap),
    "; seed ",
    witness$seed,
    " moved ",
    witness$alert_high_rows_changed,
    " alert_high rows and ",
    witness$alert_income_rows_changed,
    " alert_income rows; sum(alert_high) changed from ",
    witness$alert_high_observed,
    " to ",
    witness$alert_high_permuted
  )
  list(witness = witness, checks = checks)
}

# One randomization replication: permute, refit, and return the studentized
# statistics of every reported contrast.
run_alert_permutation <- function(
  design,
  area_days,
  terms,
  contrasts,
  bin_resolution,
  seed,
  outcome = "transports_analysis"
) {
  if (!outcome %in% names(design)) {
    rlang::abort(glue::glue("Outcome column `{outcome}` is absent."))
  }
  design$transports_analysis <- design[[outcome]]
  permuted <- permute_alert_design(design, area_days, seed, terms)
  fit <- fit_paper_c_power(
    permuted,
    paper_c_result2_formula(terms, bin_resolution)
  )
  # The permuted fits need only the unrestricted clustered covariance, which
  # fixest computes from its own scores; the weighted demeaning of
  # paper_c_score_pieces() is what the null-imposed bootstrap needs and would
  # double the cost of a replication here.
  vcov <- stats::vcov(
    fit,
    cluster = ~area_primary,
    ssc = paper_c_power_ssc()
  )[terms, terms, drop = FALSE]
  estimates <- stats::coef(fit)[terms]
  purrr::imap(contrasts, function(contrast, label) {
    estimate <- sum(contrast * estimates)
    se <- sqrt(as.numeric(crossprod(contrast, vcov %*% contrast)))
    tibble::tibble(
      contrast = label,
      seed = seed,
      estimate = estimate,
      se = se,
      statistic = estimate / se
    )
  }) |>
    purrr::list_rbind()
}

# The whole randomization distribution, one forked worker per replication
# chunk. Seeds are drawn once in the parent from an L'Ecuyer stream, so the
# result depends on `seed` alone and not on the worker count.
# Content hash of the columns an estimate depends on. A resume that reuses
# work computed from different data is the failure this exists to prevent, and
# summary statistics do not detect it: the same rows in a different order, or
# two outcomes exchanged, produce identical counts and totals.
paper_c_design_hash <- function(design, columns = NULL) {
  columns <- columns %||% sort(names(design))
  missing <- setdiff(columns, names(design))
  if (length(missing) > 0L) {
    rlang::abort(paste(
      "Cannot fingerprint; the design lacks:",
      paste(missing, collapse = ", ")
    ))
  }
  list(
    rows = nrow(design),
    columns = columns,
    hash = rlang::hash(design[, columns, drop = FALSE])
  )
}

# Seeds for the replications, drawn once from an L'Ecuyer stream so that the
# grid depends on `seed` alone and not on the worker count, the chunking, or a
# resume. Because the draw is sequential, a shorter grid is a prefix of a
# longer one with the same seed: the first 499 seeds of a 4,999-replication
# run are exactly the 499 of a 499-replication run. The partial file still
# refuses a resume across different `replications`, which is stricter than it
# needs to be but costs nothing.
randomization_seeds <- function(seed, replications) {
  withr::with_preserve_seed({
    set.seed(seed, kind = "L'Ecuyer-CMRG")
    sample.int(.Machine$integer.max, replications)
  })
}

# The settings a partial randomization file must share with the run resuming
# from it. A permutation is only comparable with another if every one of
# these is the same.
randomization_partial_settings <- function(
  label,
  terms,
  bin_resolution,
  strata_vars,
  seed,
  replications,
  design = NULL,
  frozen_years = integer(0),
  outcome = "transports_analysis"
) {
  settings <- list(
    label = label,
    terms = terms,
    bin_resolution = bin_resolution,
    strata_vars = strata_vars,
    seed = seed,
    replications = replications,
    frozen_years = frozen_years,
    outcome = outcome
  )
  # A fingerprint of the rows themselves, so that a partial file cannot be
  # resumed against a rebuilt store that changed the design under the same
  # label. Counts and totals are not enough — swapping two rows' outcomes
  # leaves every one of them unchanged — so this hashes the content of every
  # column the model and the strata use.
  if (!is.null(design)) {
    settings$design <- paper_c_design_hash(design)
  }
  settings
}

# Reject invalid draws before saving a chunk or summarising a distribution.
assert_finite_statistics <- function(chunk_results, label) {
  invalid <- chunk_results[!is.finite(chunk_results$statistic), , drop = FALSE]
  if (nrow(invalid) > 0L) {
    affected <- paste0("seed=", invalid$seed, ", contrast=", invalid$contrast)
    rlang::abort(glue::glue(
      "Non-finite randomization statistics for {label}: ",
      "{paste(affected, collapse = '; ')}. ",
      "A resume from the partial file will hit the same replications because ",
      "seeds are deterministic; fix the cause rather than restart the run."
    ))
  }
  invisible(chunk_results)
}

check_randomization_finiteness <- function() {
  finite <- tibble::tibble(
    seed = c(101L, 102L),
    contrast = c("alert_day", "alert_trial"),
    statistic = c(1, 2)
  )
  invalid <- finite
  invalid$statistic <- c(NaN, Inf)
  chunk_error <- tryCatch(
    assert_finite_statistics(invalid, "synthetic chunk"),
    error = identity
  )
  summary_error <- tryCatch(
    summarise_randomization_inference(invalid[1L, ], finite),
    error = identity
  )
  expected <- c(
    "Non-finite randomization statistics for synthetic chunk",
    "seed=101, contrast=alert_day",
    "seed=102, contrast=alert_trial",
    "seeds are deterministic",
    "fix the cause rather than restart"
  )
  checks <- tibble::tibble(
    check = c(
      "finite chunk",
      "NaN and Inf chunk rejected",
      "NaN summary rejected"
    ),
    passed = c(
      identical(assert_finite_statistics(finite, "finite chunk"), finite),
      inherits(chunk_error, "error") &&
        all(vapply(
          expected,
          grepl,
          logical(1),
          x = conditionMessage(chunk_error),
          fixed = TRUE
        )),
      inherits(summary_error, "error") &&
        grepl(
          "Non-finite randomization statistics for summary",
          conditionMessage(summary_error),
          fixed = TRUE
        ) &&
        grepl(
          "seed=101, contrast=alert_day",
          conditionMessage(summary_error),
          fixed = TRUE
        )
    )
  )
  if (!all(checks$passed)) {
    rlang::abort("Randomization finiteness synthetic checks failed.")
  }
  checks
}

# Inputs of one randomization distribution: the design narrowed to the
# columns the model and the strata need, the permutation strata, and the
# settings every stored draw must share with the current inputs, including a
# fingerprint of the rows themselves.
prepare_randomization_inputs <- function(
  design,
  terms,
  bin_resolution,
  replications,
  seed,
  strata_vars = c("area_primary", "year", "area_bin_t_1c"),
  frozen_years = integer(0),
  label = "",
  outcome = "transports_analysis"
) {
  if (!outcome %in% names(design)) {
    rlang::abort(glue::glue("Outcome column `{outcome}` is absent."))
  }
  design$transports_analysis <- design[[outcome]]
  # Only the columns the model and the strata need are carried into the forked
  # workers; a permutation copies the design, and the panel's other columns
  # would multiply that copy by the worker count for nothing.
  registry <- paper_c_alert_term_registry()
  # Unnamed: a named vector inside `dplyr::all_of()` renames the column.
  modifier_inputs <- unname(purrr::map_chr(
    registry[intersect(terms, names(registry))],
    "input"
  ))
  design <- design |>
    dplyr::select(dplyr::all_of(unique(c(
      "transports_analysis",
      outcome,
      "area_primary",
      "date",
      "year",
      terms,
      modifier_inputs,
      paper_c_result2_fixef(bin_resolution),
      strata_vars
    ))))
  area_days <- build_alert_permutation_strata(
    design,
    strata_vars,
    frozen_years
  )
  settings <- randomization_partial_settings(
    label,
    terms,
    bin_resolution,
    strata_vars,
    seed,
    replications,
    design = design,
    frozen_years = frozen_years,
    outcome = outcome
  )
  list(design = design, area_days = area_days, settings = settings)
}

# The replications for one block of seeds, one forked worker per seed.
run_randomization_block <- function(inputs, seeds, contrasts, workers) {
  settings <- inputs$settings
  chunk <- parallel::mclapply(
    seeds,
    function(replication_seed) {
      run_alert_permutation(
        inputs$design,
        inputs$area_days,
        settings$terms,
        contrasts,
        settings$bin_resolution,
        replication_seed,
        outcome = settings$outcome
      )
    },
    mc.cores = workers,
    mc.preschedule = TRUE
  )
  # A forked worker that dies returns NULL rather than a condition, and
  # list_rbind() would drop it without a word, so the shape is checked too.
  failed <- purrr::map_lgl(chunk, function(x) {
    inherits(x, "try-error") ||
      !is.data.frame(x) ||
      nrow(x) != length(contrasts)
  })
  if (any(failed)) {
    rlang::abort(glue::glue(
      "{sum(failed)} randomization replications failed or came back short; ",
      "first: {paste(utils::capture.output(print(chunk[failed][[1]]))[1], collapse = '')}"
    ))
  }
  chunk_results <- purrr::list_rbind(chunk)
  assert_finite_statistics(chunk_results, settings$label)
}

# The whole randomization distribution, in chunks, with the accumulated
# replications written out after every chunk. At four to five thousand
# replications the run takes more than a day, so it has to survive being
# interrupted: a resume reads the partial file, checks that it belongs to
# exactly this grid, and runs only the seeds it does not already hold.
run_randomization_inference <- function(
  design,
  terms,
  contrasts,
  bin_resolution,
  replications,
  seed,
  workers,
  strata_vars = c("area_primary", "year", "area_bin_t_1c"),
  frozen_years = integer(0),
  label = "",
  partial_path = NULL,
  chunk_size = 200L,
  outcome = "transports_analysis"
) {
  inputs <- prepare_randomization_inputs(
    design,
    terms,
    bin_resolution,
    replications,
    seed,
    strata_vars = strata_vars,
    frozen_years = frozen_years,
    label = label,
    outcome = outcome
  )
  area_days <- inputs$area_days
  settings <- inputs$settings
  seeds <- randomization_seeds(seed, replications)

  completed <- NULL
  if (!is.null(partial_path) && file.exists(partial_path)) {
    partial <- readRDS(partial_path)
    stored <- partial[[label]]
    if (!is.null(stored)) {
      if (!identical(stored$settings, settings)) {
        rlang::abort(glue::glue(
          "The partial randomization file for {label} was written under ",
          "different settings; delete it or fix the arguments."
        ))
      }
      if (!all(stored$replications$seed %in% seeds)) {
        rlang::abort(glue::glue(
          "The partial randomization file for {label} holds seeds this grid ",
          "does not assign."
        ))
      }
      # A partial file written before the chunk-level check may still carry
      # a non-finite draw; refuse it here rather than after the remaining
      # chunks have been computed.
      completed <- assert_finite_statistics(stored$replications, label)
      message(
        "Resuming ",
        label,
        " from ",
        dplyr::n_distinct(completed$seed),
        " completed replications"
      )
    }
  }
  remaining <- setdiff(seeds, unique(completed$seed))

  results <- completed
  chunks <- split(remaining, ceiling(seq_along(remaining) / chunk_size))
  for (index in seq_along(chunks)) {
    chunk_results <- run_randomization_block(
      inputs,
      chunks[[index]],
      contrasts,
      workers
    )
    results <- dplyr::bind_rows(results, chunk_results)
    if (!is.null(partial_path)) {
      store <- if (file.exists(partial_path)) readRDS(partial_path) else list()
      store[[label]] <- list(settings = settings, replications = results)
      save_rds_atomic(store, partial_path)
    }
    message(
      label,
      ": ",
      dplyr::n_distinct(results$seed),
      " of ",
      replications,
      " replications (chunk ",
      index,
      " of ",
      length(chunks),
      ")"
    )
  }

  list(
    strata = summarise_permutation_strata(area_days),
    strata_vars = settings$strata_vars,
    frozen_years = settings$frozen_years,
    replications = results
  )
}

# Two-sided randomization p-value per contrast, against the observed
# statistics.
summarise_randomization_inference <- function(replications, observed) {
  assert_finite_statistics(replications, "summary")
  replications |>
    dplyr::inner_join(
      observed |> dplyr::select(contrast, observed = statistic),
      by = dplyr::join_by(contrast)
    ) |>
    dplyr::summarise(
      replications = dplyr::n(),
      p_value = (sum(abs(statistic) >= abs(dplyr::first(observed))) + 1) /
        (dplyr::n() + 1),
      .by = contrast
    )
}

# Specification runner ---------------------------------------------------------

# Contrast vectors over `terms`, named as they are reported. The 2020 total is
# the sum of the nationwide coefficient and the trial deviation, which is what
# the trial areas actually experienced in 2020.
paper_c_result2_contrasts <- function(terms) {
  make_contrast <- function(values) {
    contrast <- stats::setNames(numeric(length(terms)), terms)
    contrast[names(values)] <- values
    contrast
  }
  if (!"alert_day" %in% terms) {
    rlang::abort("Result 2 contrasts require `alert_day`.")
  }
  contrasts <- list(alert_day = make_contrast(c(alert_day = 1)))
  if ("alert_trial" %in% terms) {
    contrasts$alert_trial <- make_contrast(c(alert_trial = 1))
    contrasts$alert_2020_total <- make_contrast(c(
      alert_day = 1,
      alert_trial = 1
    ))
  }
  if ("alert_high" %in% terms) {
    contrasts$alert_high <- make_contrast(c(alert_high = 1))
  }
  # Lead and lag of the alert (Issue #43 item 6): unit contrasts, reported
  # alongside alert_day and never combined with it.
  for (term in intersect(
    c("alert_lead1", "alert_lag1", "emergency", "priority"),
    terms
  )) {
    contrasts[[term]] <- make_contrast(stats::setNames(1, term))
  }
  if ("alert_income" %in% terms) {
    contrasts$alert_income <- make_contrast(c(alert_income = 1))
  }
  contrasts
}

# Restricted design and formula for the hypothesis contrast' beta = 0. The
# regressors are the original ones rotated onto the orthogonal complement of
# the contrast, so the restricted fit spans exactly the coefficient vectors
# the null allows and nothing else.
paper_c_restricted_fit <- function(
  design,
  terms,
  contrast,
  bin_resolution,
  drop_local = FALSE
) {
  basis <- qr.Q(qr(matrix(contrast, ncol = 1L)), complete = TRUE)[,
    -1L,
    drop = FALSE
  ]
  rotated <- as.matrix(design[, terms, drop = FALSE]) %*% basis
  restricted_terms <- paste0("restricted_", seq_len(ncol(rotated)))
  restricted <- design
  for (j in seq_along(restricted_terms)) {
    restricted[[restricted_terms[[j]]]] <- rotated[, j]
  }
  fit_paper_c_power(
    restricted,
    paper_c_result2_formula(restricted_terms, bin_resolution, drop_local)
  )
}

# Estimate one specification: the contrasts, their clustered CR1 statistics,
# and optionally the score wild cluster bootstrap (off by default; see the
# header for why it is only a heuristic sensitivity here). The randomization
# p-value and the jackknife interval are added by the driver, because both
# cost refits and are computed for the primary samples only.
run_paper_c_result2_spec <- function(
  design,
  label,
  terms = c("alert_day", "alert_trial"),
  bin_resolution = "0.1",
  bootstrap_replications = 0L,
  seed = 20260910L,
  outcome = "transports_analysis",
  drop_local = FALSE
) {
  if (!outcome %in% names(design)) {
    rlang::abort(glue::glue("Outcome column `{outcome}` is absent."))
  }
  design$transports_analysis <- design[[outcome]]
  formula <- paper_c_result2_formula(terms, bin_resolution, drop_local)
  fixef_vars <- paper_c_result2_fixef(bin_resolution, drop_local)
  seconds_fit <- system.time({
    fit <- fit_paper_c_power(design, formula)
  })[["elapsed"]]
  kept <- fit_kept_rows(fit, nrow(design))
  design_kept <- design[kept, ]
  mu <- as.numeric(fit$fitted.values)
  pieces <- paper_c_score_pieces(
    design_kept,
    terms,
    mu,
    fixef_vars,
    design_kept$area_primary
  )
  piece_check <- check_score_pieces(pieces, fit)
  ssc_scalar <- paper_c_ssc_scalar(pieces, fit)
  vcov <- paper_c_cluster_vcov(pieces, ssc_scalar)
  estimates <- stats::coef(fit)[terms]
  n_clusters <- pieces$n_clusters
  critical_value <- stats::qt(0.975, df = n_clusters - 1L)
  contrasts <- paper_c_result2_contrasts(terms)

  coefficients <- purrr::imap(contrasts, function(contrast, contrast_label) {
    estimate <- sum(contrast * estimates)
    se <- sqrt(as.numeric(crossprod(contrast, vcov %*% contrast)))
    statistic <- estimate / se
    row <- tibble::tibble(
      specification = label,
      contrast = contrast_label,
      estimate = estimate,
      se = se,
      statistic = statistic,
      clusters = n_clusters,
      cr1_critical_value = critical_value,
      p_value_cr1 = 2 * stats::pt(-abs(statistic), df = n_clusters - 1L),
      conf_low = estimate - critical_value * se,
      conf_high = estimate + critical_value * se
    )
    if (bootstrap_replications > 0L) {
      restricted <- paper_c_restricted_fit(
        design_kept,
        terms,
        contrast,
        bin_resolution,
        drop_local
      )
      if (
        restricted$nobs != nrow(design_kept) ||
          length(stats::coef(restricted)) != length(terms) - 1L ||
          !is.null(restricted$collin.var)
      ) {
        rlang::abort(glue::glue(
          "The restricted fit for {contrast_label} kept {restricted$nobs} of ",
          "{nrow(design_kept)} rows with ",
          "{length(stats::coef(restricted))} coefficients; expected ",
          "{nrow(design_kept)} and {length(terms) - 1L} with no collinearity."
        ))
      }
      pieces_null <- paper_c_score_pieces(
        design_kept,
        terms,
        as.numeric(restricted$fitted.values),
        fixef_vars,
        design_kept$area_primary
      )
      # The restricted fit's own first-order conditions: the total score must
      # vanish in every direction the null leaves free. Nothing else checks
      # the restricted branch, since it has no fixest object to agree with.
      free_score <- crossprod(
        qr.Q(qr(matrix(contrast, ncol = 1L)), complete = TRUE)[,
          -1L,
          drop = FALSE
        ],
        colSums(pieces_null$scores)
      )
      if (max(abs(free_score)) > 1e-4 * max(abs(pieces_null$scores))) {
        rlang::abort(glue::glue(
          "The restricted fit for {contrast_label} does not satisfy its own ",
          "first-order conditions (max free score {signif(max(abs(free_score)), 3)})."
        ))
      }
      bootstrap <- score_wild_cluster_p(
        pieces_null,
        contrast,
        ssc_scalar,
        statistic,
        bootstrap_replications,
        seed
      )
      interval <- score_wild_cluster_interval(
        pieces,
        contrast,
        ssc_scalar,
        estimate,
        se,
        bootstrap_replications,
        seed + 1L
      )
      row <- dplyr::bind_cols(
        row,
        tibble::tibble(
          p_value_wcb = bootstrap$p_value,
          wcb_replications = bootstrap$replications,
          wcb_bootstrap_sd = bootstrap$bootstrap_sd,
          wcb_critical_value = interval$critical_value
        )
      )
    }
    row |>
      dplyr::mutate(
        irr_minus_1 = expm1(estimate),
        irr_conf_low = expm1(conf_low),
        irr_conf_high = expm1(conf_high)
      )
  }) |>
    purrr::list_rbind()

  list(
    label = label,
    terms = terms,
    design_hash = paper_c_design_hash(design),
    estimates = estimates,
    bin_resolution = bin_resolution,
    rows = nrow(design),
    rows_estimated = nrow(design_kept),
    clusters = n_clusters,
    alerted_clusters = dplyr::n_distinct(
      design_kept$area_primary[design_kept$alert_day == 1]
    ),
    treated_rows = sum(design_kept$alert_day == 1),
    treated_high_rows = sum(
      design_kept$alert_day == 1 & design_kept$local_high
    ),
    transports = sum(design_kept$transports_analysis),
    removals = fit_removals(fit),
    ssc_scalar = ssc_scalar,
    piece_check = piece_check,
    coefficients = coefficients,
    seconds_fit = seconds_fit,
    fit = fit,
    design_kept_rows = kept
  )
}

# Leave-one-area-out jackknife -------------------------------------------------

# CV3J: refit the model once per cluster with that cluster dropped, and build
# the covariance from the spread of the leave-one-out estimates around their
# own mean (MacKinnon, Nielsen and Webb's CV3J; centring on the full-sample
# estimate instead changes the standard errors here by under a percent). This is the
# correction for the downward bias in CR1 that 7.12 measured under
# independent draws, and it is what carries the reported interval, since a
# permutation test of a count outcome inverts into no interval.
paper_c_jackknife_vcov <- function(
  design,
  terms,
  bin_resolution,
  workers = 1L,
  outcome = "transports_analysis"
) {
  if (!outcome %in% names(design)) {
    rlang::abort(glue::glue("Outcome column `{outcome}` is absent."))
  }
  design$transports_analysis <- design[[outcome]]
  clusters <- sort(unique(design$area_primary))
  estimates <- parallel::mclapply(
    clusters,
    function(cluster) {
      fit <- fit_paper_c_power(
        design[design$area_primary != cluster, ],
        paper_c_result2_formula(terms, bin_resolution)
      )
      stats::coef(fit)[terms]
    },
    mc.cores = workers,
    mc.preschedule = TRUE
  )
  failed <- purrr::map_lgl(estimates, function(x) {
    inherits(x, "try-error") || length(x) != length(terms) || anyNA(x)
  })
  if (any(failed)) {
    rlang::abort(glue::glue(
      "{sum(failed)} leave-one-area-out fits failed or lost a coefficient."
    ))
  }
  estimate_matrix <- do.call(rbind, estimates)
  colnames(estimate_matrix) <- terms
  rownames(estimate_matrix) <- clusters
  n_clusters <- length(clusters)
  centered <- sweep(estimate_matrix, 2L, colMeans(estimate_matrix))
  vcov <- crossprod(centered) * (n_clusters - 1L) / n_clusters
  list(
    vcov = vcov,
    estimates = estimate_matrix,
    n_clusters = n_clusters
  )
}

# CV3 statistics for the same contrasts, next to the CR1 ones. `estimates` is
# the full-sample coefficient vector; the jackknife supplies only the spread.
summarise_jackknife <- function(jackknife, estimates, contrasts) {
  critical_value <- stats::qt(0.975, df = jackknife$n_clusters - 1L)
  purrr::imap(contrasts, function(contrast, label) {
    estimate <- sum(contrast * estimates)
    se <- sqrt(as.numeric(crossprod(
      contrast,
      jackknife$vcov %*% contrast
    )))
    statistic <- estimate / se
    tibble::tibble(
      contrast = label,
      se_cv3 = se,
      statistic_cv3 = statistic,
      p_value_cv3 = 2 *
        stats::pt(-abs(statistic), df = jackknife$n_clusters - 1L),
      conf_low_cv3 = estimate - critical_value * se,
      conf_high_cv3 = estimate + critical_value * se,
      irr_conf_low_cv3 = expm1(estimate - critical_value * se),
      irr_conf_high_cv3 = expm1(estimate + critical_value * se)
    )
  }) |>
    purrr::list_rbind()
}
