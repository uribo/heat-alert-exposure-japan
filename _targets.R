library(targets)

tar_option_set(
  packages = c(
    "curl",
    "digest",
    "dplyr",
    "glue",
    "jmastats",
    "jsonlite",
    "lubridate",
    "pointblank",
    "purrr",
    "readr",
    "rlang",
    "sf",
    "stringr",
    "tibble",
    "tidyr",
    "zipangu"
  ),
  format = "rds"
)

tar_source(c(
  "R/public_pipeline.R",
  "R/alert_area.R",
  "R/moe_wbgt_api.R",
  "R/moe_forecast_threshold.R"
))

list(
  tar_target(
    source_hashes_file,
    "data-raw/source_sha256.tsv",
    format = "file"
  ),
  tar_target(source_hashes, read_source_hashes(source_hashes_file)),
  tar_target(
    n03_archive,
    download_verified(
      "https://nlftp.mlit.go.jp/ksj/gml/data/N03/N03-2024/N03-20240101_GML.zip",
      "data-raw/downloads/N03-20240101_GML.zip",
      source_hash(source_hashes, "N03-20240101_GML.zip")
    ),
    format = "file"
  ),
  tar_target(
    jma_alert_archive,
    download_verified(
      "https://www.data.jma.go.jp/developer/gis/20190125_AreaForecastLocalM_prefecture_GIS.zip",
      "data-raw/downloads/20190125_AreaForecastLocalM_prefecture_GIS.zip",
      source_hash(
        source_hashes,
        "20190125_AreaForecastLocalM_prefecture_GIS.zip"
      )
    ),
    format = "file"
  ),
  tar_target(
    n03_files,
    extract_zip(n03_archive, "data-raw/downloads/n03-2024"),
    format = "file"
  ),
  tar_target(
    jma_alert_files,
    extract_jma_zip(jma_alert_archive, "data-raw/downloads/jma-alert-20190125"),
    format = "file"
  ),
  tar_target(n03_geometry, read_n03_geometry(n03_files)),
  tar_target(
    fd_municipality_file,
    "data-raw/release/fd_municipality.csv",
    format = "file"
  ),
  tar_target(
    fd_municipality,
    readr::read_csv(fd_municipality_file, col_types = "ccc")
  ),
  tar_target(
    sample_departments_file,
    "data-raw/release/sample_departments.csv",
    format = "file"
  ),
  tar_target(
    sample_departments,
    readr::read_csv(sample_departments_file, col_types = "cc")
  ),
  tar_target(fd_geometry, build_fd_geometry(n03_geometry, fd_municipality)),
  tar_target(
    alert_records_file,
    "data-raw/release/moe_alert_historical.csv",
    format = "file"
  ),
  tar_target(
    alert_records,
    read_alert_records(alert_records_file)
  ),
  tar_target(alert_geometry, read_jma_alert_geometry(jma_alert_files)),
  tar_target(
    jma_stations,
    {
      if (utils::packageVersion("jmastats") != "0.3.0") {
        rlang::abort("The pipeline requires jmastats 0.3.0.")
      }
      jmastats::stations |>
        dplyr::select(
          area,
          station_no,
          station_type,
          station_name,
          pref_code,
          block_no
        ) |>
        dplyr::distinct(area, station_no, .keep_all = TRUE) |>
        pointblank::row_count_match(1286L)
    }
  ),
  tar_target(
    wbgt_station_roster_file,
    "data-raw/wbgt_stations2023.csv",
    format = "file"
  ),
  tar_target(
    wbgt_station_roster,
    readr::read_csv(wbgt_station_roster_file, col_types = "cdcc")
  ),
  tar_target(
    wbgt_station_geometry,
    jma_stations |>
      dplyr::mutate(
        station_no = dplyr::case_when(
          station_name == "銚子" ~ 45147,
          .default = station_no
        )
      ) |>
      dplyr::inner_join(
        wbgt_station_roster,
        by = c("area", "station_no", "station_name", "pref_code")
      ) |>
      dplyr::mutate(station_no = as.character(station_no)) |>
      # Keep the WGS 84 coordinates of jmastats::stations: the private
      # pipeline measures department-to-station distances in that CRS.
      pointblank::row_count_match(840L)
  ),
  tar_target(
    wbgt_daily_file,
    "data-raw/release/moe_wbgt_daily_agg.csv.gz",
    format = "file"
  ),
  tar_target(
    wbgt_daily,
    readr::read_csv(wbgt_daily_file, show_col_types = FALSE)
  ),
  tar_target(station_coverage, build_station_coverage(wbgt_daily)),
  tar_target(
    fd_alert_area,
    assign_fd_alert_area(
      sample_departments = sample_departments,
      fd_geometry = fd_geometry,
      alert_geometry = alert_geometry,
      wbgt_station_geometry = wbgt_station_geometry,
      station_coverage = station_coverage,
      straddle_share = 0.01,
      expected_departments = 417L,
      expected_straddles = 1L,
      geometry_prefectures = c("北海道", "沖縄県")
    )
  ),
  tar_target(
    fd_alert_area_expected_file,
    "data-raw/release/fd_alert_area.csv",
    format = "file"
  ),
  tar_target(
    fd_alert_area_expected,
    readr::read_csv(fd_alert_area_expected_file, show_col_types = FALSE)
  ),
  tar_target(
    fd_alert_area_validation,
    validate_fd_alert_area(fd_alert_area, fd_alert_area_expected)
  ),
  tar_target(
    fd_population_exclusions_file,
    "data-raw/release/fd_population_exclusions.csv",
    format = "file"
  ),
  tar_target(
    fd_population_exclusions,
    readr::read_csv(fd_population_exclusions_file, col_types = "ccc")
  ),
  tar_target(
    fd_shared_municipality_weights_file,
    "data-raw/release/fd_shared_municipality_weights.csv",
    format = "file"
  ),
  tar_target(
    fd_shared_municipality_weights,
    readr::read_csv(fd_shared_municipality_weights_file, show_col_types = FALSE)
  ),
  tar_target(
    population_2020_file,
    "data-raw/release/population_2020.csv",
    format = "file"
  ),
  tar_target(
    population_2020,
    readr::read_csv(population_2020_file, show_col_types = FALSE)
  ),
  tar_target(
    fd_population,
    build_fd_population(
      fd_municipality,
      fd_population_exclusions,
      population_2020,
      fd_shared_municipality_weights
    )
  ),
  tar_target(
    taxable_income_checksum_file,
    "data-raw/frozen/SHA256SUMS",
    format = "file"
  ),
  tar_target(
    taxable_income_integrity,
    verify_checksum_file(taxable_income_checksum_file, "data-raw/frozen")
  ),
  tar_target(
    forecast_0517_checksum_file,
    "data-raw/frozen/moe_forecast_SHA256SUMS",
    format = "file"
  ),
  tar_target(
    forecast_0517_integrity,
    verify_checksum_file(
      forecast_0517_checksum_file,
      "data-raw/frozen/moe_forecast"
    )
  ),
  tar_target(
    forecast_0416_checksum_file,
    "data-raw/frozen/moe_forecast_0416_SHA256SUMS",
    format = "file"
  ),
  tar_target(
    forecast_0416_integrity,
    verify_checksum_file(
      forecast_0416_checksum_file,
      "data-raw/frozen/moe_forecast"
    )
  ),
  tar_target(
    forecast_0517_manifest_file,
    "data-raw/frozen/moe_forecast_manifest.tsv",
    format = "file"
  ),
  tar_target(
    forecast_daily_0517,
    read_moe_forecast_set(
      forecast_0517_manifest_file,
      expected_files = 609L
    )
  ),
  tar_target(
    forecast_0416_manifest_file,
    "data-raw/frozen/moe_forecast_0416_manifest.tsv",
    format = "file"
  ),
  tar_target(
    forecast_daily_0416,
    read_moe_forecast_set(
      forecast_0416_manifest_file,
      expected_files = 608L
    )
  ),
  tar_target(
    alert_area_roster_file,
    "data-raw/release/alert_area_roster.csv",
    format = "file"
  ),
  tar_target(
    alert_area_roster,
    readr::read_csv(alert_area_roster_file, show_col_types = FALSE)
  ),
  tar_target(
    forecast_area_day,
    build_moe_forecast_area_day(
      daily_sets = list(
        "05/17" = forecast_daily_0517,
        "04/16" = forecast_daily_0416
      ),
      wbgt_stations = wbgt_station_geometry,
      heat_alert = alert_geometry,
      alert_records = alert_records,
      area_roster = alert_area_roster
    ) |>
      pointblank::row_count_match(16592L) |>
      pointblank::rows_distinct(
        columns = c(set, alert_area, target_date)
      )
  ),
  tar_target(
    forecast_area_day_expected_file,
    "data-raw/release/moe_forecast_area_day.csv",
    format = "file"
  ),
  tar_target(
    forecast_area_day_expected,
    readr::read_csv(forecast_area_day_expected_file, show_col_types = FALSE)
  ),
  tar_target(
    forecast_area_day_validation,
    validate_forecast_area_day(forecast_area_day, forecast_area_day_expected)
  )
)
