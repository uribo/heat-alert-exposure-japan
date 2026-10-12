# NOTE: readable analysis code from the private project. It needs the fire-department
# transport records, which are not distributed, and is not sourced by the public
# pipeline (_targets.R). Kept so that the published estimates can be audited.

# Multiplicity helpers for paper C's five pre-specified outcomes.

holm_adjust <- function(p) {
  if (is.null(names(p)) || any(names(p) == "") || anyDuplicated(names(p))) {
    rlang::abort("`p` must be a uniquely named numeric vector.")
  }
  tibble::tibble(
    outcome = names(p),
    p_raw = as.numeric(p),
    p_holm = stats::p.adjust(as.numeric(p), method = "holm")
  )
}

westfall_young_stepdown <- function(observed, permuted) {
  required_observed <- c("outcome", "statistic")
  required_permuted <- c("outcome", "seed", "statistic")
  if (!all(required_observed %in% names(observed))) {
    rlang::abort("`observed` must contain outcome and statistic.")
  }
  if (!all(required_permuted %in% names(permuted))) {
    rlang::abort("`permuted` must contain outcome, seed, and statistic.")
  }
  if (
    nrow(observed) == 0L ||
      anyDuplicated(observed$outcome) ||
      anyNA(observed[, required_observed]) ||
      anyNA(permuted[, required_permuted]) ||
      anyDuplicated(permuted[c("outcome", "seed")])
  ) {
    rlang::abort("Multiplicity inputs must be complete and uniquely keyed.")
  }
  if (!setequal(observed$outcome, unique(permuted$outcome))) {
    rlang::abort("Observed and permuted outcomes differ.")
  }
  seed_sets <- split(permuted$seed, permuted$outcome) |>
    lapply(sort)
  if (!all(vapply(seed_sets[-1L], identical, logical(1), seed_sets[[1L]]))) {
    rlang::abort("Every outcome must use the same permutation seeds.")
  }

  ordered <- observed |>
    dplyr::mutate(statistic = abs(statistic)) |>
    dplyr::arrange(dplyr::desc(statistic)) |>
    dplyr::mutate(rank = dplyr::row_number())
  ranked_permuted <- permuted |>
    dplyr::mutate(statistic = abs(statistic)) |>
    dplyr::left_join(
      ordered |> dplyr::select(outcome, rank),
      by = dplyr::join_by(outcome),
      relationship = "many-to-one"
    )
  b <- length(seed_sets[[1L]])
  raw <- ranked_permuted |>
    dplyr::inner_join(
      ordered |> dplyr::select(outcome, observed = statistic),
      by = dplyr::join_by(outcome),
      relationship = "many-to-one"
    ) |>
    dplyr::summarise(
      p_raw = (sum(statistic >= dplyr::first(observed)) + 1) / (b + 1),
      .by = outcome
    )
  adjusted <- purrr::map_dfr(seq_len(nrow(ordered)), function(j) {
    threshold <- ordered$statistic[[j]]
    maxima <- ranked_permuted |>
      dplyr::filter(rank >= j) |>
      dplyr::summarise(maximum = max(statistic), .by = seed)
    tibble::tibble(
      rank = j,
      p_wy_unconstrained = (sum(maxima$maximum >= threshold) + 1) / (b + 1)
    )
  }) |>
    dplyr::mutate(p_wy = cummax(p_wy_unconstrained))

  ordered |>
    dplyr::left_join(raw, by = dplyr::join_by(outcome)) |>
    dplyr::left_join(adjusted, by = dplyr::join_by(rank)) |>
    dplyr::select(outcome, statistic, p_raw, p_wy)
}

check_shared_permutations <- function(hashes) {
  required <- c("design_hash", "strata_hash", "seeds_hash")
  if (
    length(hashes) == 0L ||
      !all(vapply(
        hashes,
        function(x) identical(sort(names(x)), sort(required)),
        logical(1)
      ))
  ) {
    rlang::abort("Each outcome must supply the three permutation hashes.")
  }
  reference <- hashes[[1L]]
  matches <- vapply(hashes, identical, logical(1), reference)
  if (!all(matches)) {
    rlang::abort(
      "Permutation design, strata, or seed hashes differ by outcome."
    )
  }
  tibble::tibble(
    check = required,
    value = vapply(
      reference,
      function(x) paste(x, collapse = ":"),
      character(1)
    ),
    passed = TRUE
  )
}

check_westfall_young <- function() {
  outcomes <- c("a", "b", "c")
  observed <- tibble::tibble(
    outcome = outcomes,
    statistic = c(2.4, 1.7, 0.9)
  )
  permuted <- tidyr::expand_grid(outcome = outcomes, seed = seq_len(50L)) |>
    dplyr::mutate(
      # Local seed: the check must not move the caller's RNG state.
      statistic = withr::with_seed(20260919L, stats::rnorm(dplyr::n()))
    )
  result <- westfall_young_stepdown(observed, permuted)
  single <- westfall_young_stepdown(
    observed |> dplyr::slice(1L),
    permuted |> dplyr::filter(outcome == outcomes[[1L]])
  )
  expected_first <- (sum(
    permuted |>
      dplyr::summarise(maximum = max(abs(statistic)), .by = seed) |>
      dplyr::pull(maximum) >=
      max(abs(observed$statistic))
  ) +
    1) /
    51
  zero <- westfall_young_stepdown(
    observed,
    permuted |> dplyr::mutate(statistic = 0)
  )
  tolerance <- sqrt(.Machine$double.eps)
  if (
    abs(single$p_wy - single$p_raw) > tolerance ||
      any(result$p_wy + tolerance < result$p_raw) ||
      any(diff(result$p_wy) < -tolerance) ||
      abs(result$p_wy[[1L]] - expected_first) > tolerance ||
      any(abs(zero$p_wy - 1 / 51) > tolerance)
  ) {
    rlang::abort("Westfall-Young synthetic checks failed.")
  }
  tibble::tibble(
    check = c(
      "single outcome",
      "adjusted not below raw",
      "step-down monotonicity",
      "first-step maximum",
      "zero permutations"
    ),
    passed = TRUE
  )
}

blind_outcome_table <- function(x) {
  forbidden <- grepl(
    "(^|_)(estimate|statistic)(_|$)|(^|_)conf_(low|high)(_|$)|^ci_|^irr_minus_1$|^p_value_(cr1|cv3)$",
    names(x)
  )
  x[, !forbidden, drop = FALSE]
}
