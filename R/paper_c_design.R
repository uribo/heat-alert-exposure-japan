# NOTE: readable analysis code from the private project. It needs the fire-department
# transport records, which are not distributed, and is not sourced by the public
# pipeline (_targets.R). Kept so that the published estimates can be audited.

# The analysis design for paper C's result 2 (paper_strategy 7.2, 7.12).
#
# The design assembly used to live in `src/paper_c_power.R`. It is here so that
# the power calculation and the baseline estimation build the same rows from
# the same guards: one object, two drivers. Every function is pure; the drivers
# read the store and pass the targets in.

# The counts of the pre-specified 23-prefecture sample (7.12). They are the
# default for the scripts in src/ that ran on that sample, which restrict the
# store's panel to paper_c_prespecified_prefectures(); the targets pass the
# counts of the current sample (Issue #59).
paper_c_design_expected_23 <- function() {
  list(
    panel_rows = 200173L,
    departments = 391L,
    full = list(
      rows = 200161L,
      cells = c(184508L, 5542L, 8478L, 1633L),
      areas = 30L,
      alerted_areas = 22L
    ),
    stratified = list(
      rows = 127544L,
      cells = c(118973L, 3155L, 4635L, 781L),
      areas = 21L,
      alerted_areas = 13L
    )
  )
}

# The 23 prefectures of the pre-specified sample (7.8). Checks against records
# made on that sample (the outcome scan) refit on these rows (Issue #59).
paper_c_prespecified_prefectures <- function() {
  c(
    "北海道",
    "岩手県",
    "山形県",
    "福島県",
    "栃木県",
    "群馬県",
    "埼玉県",
    "神奈川県",
    "新潟県",
    "福井県",
    "山梨県",
    "長野県",
    "静岡県",
    "愛知県",
    "三重県",
    "兵庫県",
    "和歌山県",
    "島根県",
    "山口県",
    "徳島県",
    "香川県",
    "愛媛県",
    "佐賀県"
  )
}

# The counts of the 25-prefecture sample that replaced it (Issue #59): 富山 and
# 沖縄 add 15,452 department-days with WBGT, four alert areas that all carry
# alerts, and 26 departments, and every added row has a severity breakdown,
# so the stratified cells grow by the same 14,431 / 110 / 857 / 54.
paper_c_design_expected_25 <- function() {
  list(
    panel_rows = 215625L,
    departments = 417L,
    full = list(
      rows = 215613L,
      cells = c(198939L, 5652L, 9335L, 1687L),
      areas = 34L,
      alerted_areas = 26L
    ),
    stratified = list(
      rows = 142996L,
      cells = c(133404L, 3265L, 5492L, 835L),
      areas = 25L,
      alerted_areas = 17L
    )
  )
}

# Assemble the estimation design from the store's targets and guard it against
# the counts in `expected`. Returns the two samples, the audit trail, and the
# summary the reports print.
build_paper_c_analysis_design <- function(
  panel,
  wbgt_daily,
  wbgt_stations,
  heat_alert,
  severity_prefectures,
  taxable_income = NULL,
  income_fiscal_year = 2020L,
  expected = paper_c_design_expected_23()
) {
  if (
    nrow(panel) != expected$panel_rows ||
      dplyr::n_distinct(panel$department) != expected$departments
  ) {
    rlang::abort(glue::glue(
      "Expected {expected$panel_rows} panel rows / {expected$departments} ",
      "departments, got {nrow(panel)} / {dplyr::n_distinct(panel$department)}."
    ))
  }

  station_areas <- assign_wbgt_station_alert_area(wbgt_stations, heat_alert)
  # WBGT series stations that have no point geometry cannot be placed; they
  # are excluded from the area maximum and reported. None of them is a
  # department's assigned station.
  stations_without_geometry <- setdiff(
    unique(as.character(wbgt_daily$station_no)),
    station_areas$assignment$station_no
  )
  if (any(panel$near_wbgt_station %in% stations_without_geometry)) {
    rlang::abort("A department's station has no geometry.")
  }
  nonfinite_wbgt_rows <- sum(!is.finite(wbgt_daily$wbgt_max))
  area_wbgt_max <- build_area_wbgt_max(wbgt_daily, station_areas$assignment)
  design_full <- build_paper_c_power_design(
    panel,
    area_wbgt_max,
    taxable_income = taxable_income,
    income_fiscal_year = income_fiscal_year,
    expected_departments = expected$departments
  )
  design_stratified <- design_full |>
    dplyr::filter(prefecture %in% severity_prefectures)

  check_paper_c_power_sample(
    design_full,
    expected_rows = expected$full$rows,
    expected_cells = expected$full$cells,
    expected_areas = expected$full$areas,
    expected_alerted_areas = expected$full$alerted_areas,
    label = "Full"
  )
  check_paper_c_power_sample(
    design_stratified,
    expected_rows = expected$stratified$rows,
    expected_cells = expected$stratified$cells,
    expected_areas = expected$stratified$areas,
    expected_alerted_areas = expected$stratified$alerted_areas,
    label = "Stratified"
  )

  designs <- list(full = design_full, stratified = design_stratified)
  bin_levels <- purrr::map(designs, function(design) {
    design |>
      dplyr::summarise(dplyr::across(
        c(
          area_bin_t,
          area_bin_lag,
          local_bin,
          area_bin_t_1c,
          area_bin_lag_1c,
          local_bin_1c
        ),
        dplyr::n_distinct
      ))
  })
  stations_per_area <- station_areas$assignment |>
    dplyr::count(alert_area, name = "stations") |>
    dplyr::mutate(in_sample = alert_area %in% design_full$area_primary)
  sample_prefectures <- unique(panel$prefecture)
  coarse_in_sample <- c("鹿児島県", "沖縄県") %in% sample_prefectures

  income_standardization <- NULL
  if (!is.null(taxable_income)) {
    income_standardization <- design_full |>
      dplyr::distinct(
        prefecture,
        department,
        income_per_taxpayer,
        income_per_resident,
        income_z,
        income_resident_z
      ) |>
      dplyr::summarise(
        departments = dplyr::n(),
        taxpayer_mean = mean(income_per_taxpayer),
        taxpayer_sd = stats::sd(income_per_taxpayer),
        resident_mean = mean(income_per_resident),
        resident_sd = stats::sd(income_per_resident),
        taxpayer_z_mean = mean(income_z),
        taxpayer_z_sd = stats::sd(income_z),
        resident_z_mean = mean(income_resident_z),
        resident_z_sd = stats::sd(income_resident_z)
      )
    if (
      income_standardization$departments != expected$departments ||
        max(abs(c(
          income_standardization$taxpayer_z_mean,
          income_standardization$resident_z_mean
        ))) >
          1e-12 ||
        max(abs(
          c(
            income_standardization$taxpayer_z_sd,
            income_standardization$resident_z_sd
          ) -
            1
        )) >
          1e-12
    ) {
      rlang::abort(
        "Income standardization did not use the sample's departments once each."
      )
    }
  }

  summary <- list(
    rows = purrr::map_int(designs, nrow),
    cells = purrr::map(designs, function(design) {
      design |>
        dplyr::count(alert_day, local_high) |>
        dplyr::arrange(alert_day, local_high)
    }),
    clusters = purrr::map_int(designs, function(design) {
      dplyr::n_distinct(design$area_primary)
    }),
    alerted_clusters = purrr::map_int(designs, function(design) {
      design |>
        dplyr::filter(alert_day == 1) |>
        dplyr::pull(area_primary) |>
        dplyr::n_distinct()
    }),
    bin_levels = bin_levels,
    stations_per_area = stations_per_area,
    station_audit = station_areas$audit,
    stations_without_geometry = stations_without_geometry,
    nonfinite_wbgt_rows = nonfinite_wbgt_rows,
    wbgt_date_range = range(wbgt_daily$date),
    severity_prefectures = severity_prefectures,
    income_fiscal_year = if (is.null(taxable_income)) {
      NULL
    } else {
      income_fiscal_year
    },
    income_standardization = income_standardization,
    coarse_label_prefectures_in_sample = stats::setNames(
      coarse_in_sample,
      c("鹿児島県", "沖縄県")
    )
  )

  list(designs = designs, summary = summary)
}
