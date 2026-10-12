# NOTE: readable analysis code from the private project. It needs the fire-department
# transport records, which are not distributed, and is not sourced by the public
# pipeline (_targets.R). Kept so that the published estimates can be audited.

# Registry for alert-treatment columns used in paper-C designs.

paper_c_alert_term_registry <- function() {
  list(
    alert_trial = list(
      input = "year",
      expression = rlang::expr(as.numeric(alert_day == 1 & year == 2020))
    ),
    alert_high = list(
      input = "local_high",
      expression = rlang::expr(as.numeric(alert_day == 1 & local_high))
    ),
    alert_income = list(
      input = "income_z",
      expression = rlang::expr(alert_day * income_z)
    ),
    alert_income_resident = list(
      input = "income_resident_z",
      expression = rlang::expr(alert_day * income_resident_z)
    )
  )
}

rebuild_alert_terms <- function(design, terms) {
  registry <- paper_c_alert_term_registry()
  unknown <- setdiff(terms, names(registry))
  if (length(unknown) > 0L) {
    rlang::abort(paste(
      "Cannot rebuild unknown alert term(s):",
      paste(unknown, collapse = ", ")
    ))
  }
  # Every expression reads `alert_day`, but none declares it as an input.
  # Without this check `rlang::eval_tidy()` would fall back to the calling
  # environment and build a treatment column from a stray object of that name.
  if (!"alert_day" %in% names(design)) {
    rlang::abort("`alert_day` is required to rebuild alert terms.")
  }
  for (term in terms) {
    input <- registry[[term]]$input
    if (!input %in% names(design)) {
      rlang::abort(glue::glue(
        "`{input}` is required to rebuild `{term}`."
      ))
    }
    design[[term]] <- rlang::eval_tidy(
      registry[[term]]$expression,
      data = design
    )
  }
  design
}
