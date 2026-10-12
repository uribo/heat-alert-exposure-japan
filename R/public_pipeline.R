# Public exposure-side pipeline helpers.

sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

verify_sha256 <- function(path, expected) {
  if (
    length(expected) != 1L ||
      is.na(expected) ||
      !grepl("^[0-9a-f]{64}$", expected)
  ) {
    rlang::abort(glue::glue(
      "No final SHA256 is recorded for {basename(path)}."
    ))
  }
  actual <- sha256_file(path)
  if (!identical(actual, expected)) {
    rlang::abort(glue::glue(
      "SHA256 mismatch for {basename(path)}: expected {expected}, got {actual}."
    ))
  }
  path
}

download_verified <- function(url, path, expected) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(path)) {
    temp <- tempfile(tmpdir = dirname(path), fileext = ".download")
    on.exit(unlink(temp), add = TRUE)
    curl::curl_download(url, temp, quiet = FALSE)
    verify_sha256(temp, expected)
    if (!file.rename(temp, path)) {
      rlang::abort(glue::glue(
        "Could not move the verified download to {path}."
      ))
    }
  }
  verify_sha256(path, expected)
}

extract_zip <- function(path, directory) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  files <- utils::unzip(path, exdir = directory)
  normalizePath(files, mustWork = TRUE)
}

extract_jma_zip <- function(path, directory) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  status <- system2("/usr/bin/ditto", c("-x", "-k", path, directory))
  if (!identical(status, 0L)) {
    rlang::abort("Could not extract the JMA forecast-area archive.")
  }
  files <- list.files(directory, full.names = TRUE)
  if (length(files) != 3L) {
    rlang::abort("The JMA forecast-area archive should contain three files.")
  }
  normalizePath(files, mustWork = TRUE)
}

read_source_hashes <- function(path = "data-raw/source_sha256.tsv") {
  readr::read_tsv(
    path,
    col_types = readr::cols(
      source = readr::col_character(),
      sha256 = readr::col_character()
    ),
    na = "TO_BE_REPLACED",
    show_col_types = FALSE
  )
}

verify_checksum_file <- function(
  checksum_path,
  directory = dirname(checksum_path)
) {
  lines <- readr::read_lines(checksum_path)
  fields <- stringr::str_match(lines, "^([0-9a-f]{64})  (.+)$")
  if (anyNA(fields[, 2:3, drop = FALSE])) {
    rlang::abort(glue::glue("Malformed checksum line in {checksum_path}."))
  }
  purrr::walk2(fields[, 3], fields[, 2], function(file, expected) {
    verify_sha256(file.path(directory, file), expected)
  })
  tibble::tibble(
    checksum_file = checksum_path,
    files = nrow(fields),
    valid = TRUE
  )
}

source_hash <- function(hashes, source) {
  value <- hashes$sha256[match(source, hashes$source)]
  if (length(value) != 1L || is.na(value)) {
    rlang::abort(glue::glue(
      "The final SHA256 for {source} is pending because this environment cannot access the official download."
    ))
  }
  value
}

read_n03_geometry <- function(files) {
  # The private pipeline reads the GeoJSON member of the same archive (WGS 84
  # text coordinates), not the shapefile (JGD2011 doubles); reading the
  # shapefile moved every centroid by about 0.5 mm and changed the released
  # distances beyond 1e-6 m. st_make_valid() as the private pipeline does.
  geojson <- files[grepl("/N03-20240101\\.geojson$", files)]
  if (length(geojson) != 1L) {
    rlang::abort("The N03 archive does not contain one N03-20240101.geojson.")
  }
  sf::st_read(geojson, as_tibble = TRUE, quiet = TRUE) |>
    pointblank::row_count_match(124133L) |>
    sf::st_make_valid() |>
    dplyr::filter(!is.na(N03_007)) |>
    dplyr::mutate(N03_007 = as.character(N03_007))
}

build_fd_geometry <- function(n03, crosswalk) {
  n03 |>
    dplyr::select(area_code = N03_007, geometry) |>
    dplyr::inner_join(
      dplyr::select(crosswalk, fd_id, prefecture, area_code),
      by = dplyr::join_by(area_code),
      relationship = "many-to-many"
    ) |>
    sf::st_as_sf() |>
    dplyr::summarise(
      geometry = sf::st_union(geometry),
      .by = c(fd_id, prefecture)
    )
}

read_jma_alert_geometry <- function(files) {
  shape <- files[grepl("\\.shp$", files, ignore.case = TRUE)]
  if (length(shape) != 1L) {
    rlang::abort("The JMA archive does not contain exactly one shapefile.")
  }
  sf::st_read(shape, quiet = TRUE) |>
    pointblank::row_count_match(64L) |>
    sf::st_set_crs(6668) |>
    dplyr::mutate(
      name = dplyr::replace_values(
        name,
        "上川地方" ~ "上川・留萌地方",
        "留萌地方" ~ "上川・留萌地方",
        "石狩地方" ~ "石狩・空知・後志地方",
        "空知地方" ~ "石狩・空知・後志地方",
        "後志地方" ~ "石狩・空知・後志地方",
        "釧路地方" ~ "釧路・根室地方",
        "根室地方" ~ "釧路・根室地方",
        "胆振地方" ~ "胆振・日高地方",
        "日高地方" ~ "胆振・日高地方",
        "渡島地方" ~ "渡島・檜山地方",
        "檜山地方" ~ "渡島・檜山地方"
      )
    ) |>
    dplyr::group_by(name) |>
    dplyr::mutate(code = min(as.numeric(code))) |>
    dplyr::group_by(code, name) |>
    dplyr::summarise(.groups = "drop") |>
    pointblank::row_count_match(58L) |>
    dplyr::arrange(code) |>
    dplyr::transmute(
      code,
      area = stringr::str_remove(name, "地方$")
    ) |>
    dplyr::mutate(
      area = zipangu::harmonize_prefecture_name(area, to = "short") |>
        stringr::str_replace("鹿児島県", "鹿児島")
    )
}

build_station_coverage <- function(wbgt_daily) {
  wbgt_daily |>
    dplyr::mutate(
      year = lubridate::year(date),
      month = lubridate::month(date)
    ) |>
    dplyr::filter(
      year %in% 2018:2022,
      month %in% 6:9,
      is.finite(wbgt_max)
    ) |>
    dplyr::summarise(
      days = dplyr::n(),
      .by = c(station_no, year)
    ) |>
    tidyr::complete(station_no, year = 2018:2022, fill = list(days = 0L)) |>
    dplyr::summarise(
      years_full = sum(days >= 100L),
      min_days = min(days),
      .by = station_no
    ) |>
    pointblank::row_count_match(846L) |>
    pointblank::rows_distinct(columns = station_no)
}

assign_wbgt_station_alert_area <- function(
  wbgt_station_geometry,
  alert_geometry
) {
  stations <- wbgt_station_geometry |>
    sf::st_transform(6933)
  areas <- alert_geometry |>
    sf::st_make_valid() |>
    sf::st_transform(6933)
  within <- sf::st_within(stations, areas)
  area_index <- vapply(
    within,
    function(index) {
      if (length(index) == 1L) index[[1L]] else NA_integer_
    },
    integer(1)
  )
  fallback <- is.na(area_index)
  area_index[fallback] <- sf::st_nearest_feature(stations[fallback, ], areas)
  list(
    assignment = stations |>
      sf::st_drop_geometry() |>
      dplyr::transmute(
        station_no,
        station_name,
        station_label = area,
        alert_area = areas$area[area_index],
        nearest_fallback = fallback
      ),
    audit = tibble::tibble(
      stations = nrow(stations),
      within_polygon = sum(!fallback),
      nearest_fallback = sum(fallback)
    )
  )
}

build_fd_population <- function(crosswalk, exclusions, population, weights) {
  population_key <- population |>
    tidyr::separate_longer_delim(n03_area_codes, delim = ";") |>
    dplyr::transmute(
      prefecture,
      population_area_code,
      area_code = as.character(n03_area_codes)
    )
  population_values <- population |>
    dplyr::select(!n03_area_codes)
  population_weights <- weights |>
    dplyr::transmute(
      fd_id,
      prefecture,
      population_area_code = area_code,
      weight
    )
  crosswalk |>
    dplyr::anti_join(exclusions, by = dplyr::join_by(fd_id, area_code)) |>
    dplyr::inner_join(
      population_key,
      by = dplyr::join_by(prefecture, area_code),
      relationship = "many-to-one"
    ) |>
    dplyr::distinct(fd_id, prefecture, population_area_code) |>
    dplyr::inner_join(
      population_values,
      by = dplyr::join_by(prefecture, population_area_code),
      relationship = "many-to-one"
    ) |>
    dplyr::left_join(
      population_weights,
      by = dplyr::join_by(fd_id, prefecture, population_area_code),
      relationship = "many-to-one"
    ) |>
    dplyr::mutate(weight = dplyr::coalesce(weight, 1)) |>
    dplyr::summarise(
      population = sum(population * weight),
      age_class1 = sum(age_class1 * weight),
      age_class2 = sum(age_class2 * weight),
      age_class3 = sum(age_class3 * weight),
      .by = c(fd_id, prefecture)
    ) |>
    dplyr::mutate(
      prop_age_class1 = age_class1 / population * 100,
      prop_age_class2 = age_class2 / population * 100,
      prop_age_class3 = age_class3 / population * 100
    )
}

read_alert_records <- function(path) {
  readr::read_csv(
    path,
    col_types = readr::cols(
      area = readr::col_character(),
      datetime = readr::col_datetime(format = "%Y-%m-%dT%H:%M:%S%z"),
      alert = readr::col_logical(),
      special_alert = readr::col_logical()
    ),
    show_col_types = FALSE
  )
}

validate_forecast_area_day <- function(actual, expected) {
  actual <- actual |>
    dplyr::arrange(set, alert_area, target_date)
  expected <- expected |>
    dplyr::arrange(set, alert_area, target_date)
  equal_columns <- identical(names(actual), names(expected)) &&
    all(vapply(
      names(expected),
      function(column) {
        isTRUE(all.equal(
          actual[[column]],
          expected[[column]],
          tolerance = 0,
          check.attributes = FALSE
        ))
      },
      logical(1)
    ))
  if (nrow(actual) != 16592L || !equal_columns) {
    rlang::abort(
      "The rebuilt forecast area-day table differs from the released table."
    )
  }
  tibble::tibble(rows = nrow(actual), equal = TRUE)
}

validate_fd_alert_area <- function(actual, expected) {
  # Rebuilt from the same N03 GeoJSON, jmastats 0.3.0 coordinates and the
  # same assignment code, the table reproduces the released one field by
  # field: every identifier, area, station and count exactly, the three area
  # shares within 1e-9 (measured 9e-12) and the distances within 1e-6 m
  # (measured 1.8e-7 m, the centroid rounding of the union).
  actual <- dplyr::arrange(actual, prefecture, fd_id)
  expected <- dplyr::arrange(expected, prefecture, fd_id)
  if (!identical(names(actual), names(expected))) {
    rlang::abort("The rebuilt assignment table has different columns.")
  }
  share_fields <- c("share_primary", "share_secondary", "share_covered")
  distance_fields <- c("dist_m", "dist_m_any")
  exact_fields <- setdiff(names(expected), c(share_fields, distance_fields))
  equal_exact <- vapply(
    exact_fields,
    function(column) {
      x <- actual[[column]]
      y <- expected[[column]]
      if (column %in% c("near_wbgt_station", "near_wbgt_station_any")) {
        x <- as.character(x)
        y <- as.character(y)
      }
      identical(is.na(x), is.na(y)) && all(x == y, na.rm = TRUE)
    },
    logical(1)
  )
  if (!all(equal_exact)) {
    rlang::abort(glue::glue(
      "Exact assignment fields differ from the released table: ",
      "{paste(exact_fields[!equal_exact], collapse = ', ')}."
    ))
  }
  max_abs_error <- function(fields) {
    x <- unlist(actual[fields])
    y <- unlist(expected[fields])
    if (!identical(is.na(x), is.na(y))) {
      rlang::abort(glue::glue(
        "Missing values differ in {paste(fields, collapse = ', ')}."
      ))
    }
    max(abs(x - y), na.rm = TRUE)
  }
  max_share_error <- max_abs_error(share_fields)
  max_distance_error <- max_abs_error(distance_fields)
  if (max_share_error > 1e-9) {
    rlang::abort(
      "A recomputed area share exceeds the absolute tolerance of 1e-9."
    )
  }
  if (max_distance_error > 1e-6) {
    rlang::abort(
      "A recomputed distance exceeds the absolute tolerance of 1e-6 m."
    )
  }
  tibble::tibble(
    rows = nrow(actual),
    exact_fields = length(exact_fields),
    exact_fields_equal = TRUE,
    max_absolute_share_error = max_share_error,
    max_distance_error_m = max_distance_error
  )
}
