# Assign anonymized service areas in the fixed paper-C sample to heat-alert areas.
#
# Prefecture names determine the assignment outside `geometry_prefectures`.
# Geometry provides the split of the prefectures that hold several alert areas
# (Hokkaido's eight, and Okinawa's four since Issue #59) and audits the other
# prefectures. Area shares use
# EPSG:6933, a global equal-area projection, because the source boundaries use
# different geographic CRSs and raw intersections create boundary slivers.
#
# The nearest station is chosen among stations that carry daily WBGT in every
# sample year (`station_coverage`, one row per station with `years_full`), not
# among the whole 2023 roster: a station opened in 2021 or 2023 is otherwise
# the geometric nearest and silently loses the earlier seasons (issue #22).
# The distance-only choice is kept as `near_wbgt_station_any` for the audit,
# and the coverage of both choices is carried as `station_years_covered` /
# `station_years_covered_any` (NA for the departments without a station).
assign_fd_alert_area <- function(
  sample_departments,
  fd_geometry,
  alert_geometry,
  wbgt_station_geometry,
  station_coverage,
  required_years = 5L,
  straddle_share = 0.01,
  expected_departments,
  expected_straddles,
  geometry_prefectures = c("北海道", "沖縄県")
) {
  sample_departments <- sample_departments |>
    dplyr::distinct(prefecture, fd_id)

  if (nrow(sample_departments) != expected_departments) {
    rlang::abort(glue::glue(
      "Expected {expected_departments} anonymized units in the paper-C sample, ",
      "got {nrow(sample_departments)}."
    ))
  }
  if (
    nrow(alert_geometry) != 58L || dplyr::n_distinct(alert_geometry$area) != 58L
  ) {
    rlang::abort("Expected 58 unique heat-alert areas.")
  }
  if (
    is.na(sf::st_crs(fd_geometry)) ||
      is.na(sf::st_crs(alert_geometry)) ||
      is.na(sf::st_crs(wbgt_station_geometry))
  ) {
    rlang::abort(
      "Fire-department, alert-area, and WBGT-station geometries need a CRS."
    )
  }

  fd_keys <- fd_geometry |>
    sf::st_drop_geometry() |>
    dplyr::distinct(prefecture, fd_id)
  missing_exact <- sample_departments |>
    dplyr::anti_join(
      fd_keys,
      by = dplyr::join_by(prefecture, fd_id)
    )

  if (nrow(missing_exact) != 0L) {
    rlang::abort(glue::glue(
      "Every sample unit should have geometry; ",
      "{nrow(missing_exact)} lack it: ",
      "{paste(missing_exact$fd_id, collapse = \", \")}."
    ))
  }

  exact_geometry <- sample_departments |>
    dplyr::inner_join(
      fd_geometry,
      by = dplyr::join_by(prefecture, fd_id)
    ) |>
    sf::st_as_sf() |>
    dplyr::transmute(
      prefecture,
      fd_id,
      geometry_source = "department",
      geometry
    )

  fd_with_geometry <- exact_geometry
  fd_centroids <- suppressWarnings(sf::st_centroid(fd_with_geometry))

  # Distance-only choice over the whole roster, kept for the audit against
  # sf_fd_merge_all, which uses the same rule.
  nearest_any <- sf::st_nearest_feature(fd_centroids, wbgt_station_geometry)
  # Coverage-aware choice: the nearest station with daily WBGT in every
  # sample year.
  covered_stations <- station_coverage |>
    dplyr::filter(years_full >= required_years) |>
    dplyr::pull(station_no) |>
    as.character()
  eligible <- wbgt_station_geometry |>
    dplyr::filter(as.character(station_no) %in% covered_stations)
  if (nrow(eligible) < 0.9 * nrow(wbgt_station_geometry)) {
    rlang::abort(glue::glue(
      "Only {nrow(eligible)} of {nrow(wbgt_station_geometry)} stations cover every sample year; check station_coverage."
    ))
  }
  nearest_covered <- sf::st_nearest_feature(fd_centroids, eligible)
  station_years <- function(station) {
    years <- station_coverage$years_full[
      match(station, as.character(station_coverage$station_no))
    ]
    dplyr::coalesce(as.integer(years), 0L)
  }
  station_assignment <- fd_with_geometry |>
    sf::st_drop_geometry() |>
    dplyr::transmute(
      prefecture,
      fd_id,
      near_wbgt_station = as.character(eligible$station_no[nearest_covered]),
      dist_m = as.numeric(sf::st_distance(
        fd_centroids,
        eligible[nearest_covered, ],
        by_element = TRUE
      )),
      near_wbgt_station_any = as.character(
        wbgt_station_geometry$station_no[nearest_any]
      ),
      dist_m_any = as.numeric(sf::st_distance(
        fd_centroids,
        wbgt_station_geometry[nearest_any, ],
        by_element = TRUE
      )),
      station_reassigned = near_wbgt_station != near_wbgt_station_any,
      # Declared coverage of the chosen station (issue #22). station_coverage
      # only lists stations with at least one 2018-2022 row, so a station that
      # opened later (秩父別 15156) matches nothing and is coalesced to 0.
      station_years_covered = station_years(near_wbgt_station),
      station_years_covered_any = station_years(near_wbgt_station_any)
    )
  # Two anonymized units move to stations with full temporal coverage. Any
  # other change means the station roster or coverage calculation moved.
  reassigned <- station_assignment |>
    dplyr::filter(station_reassigned) |>
    dplyr::pull(fd_id)
  if (!setequal(reassigned, c("FD026", "FD119"))) {
    rlang::abort(paste(
      "Coverage-aware reassignment changed:",
      paste(reassigned, collapse = ", ")
    ))
  }

  projected_fd <- fd_with_geometry |>
    sf::st_make_valid() |>
    sf::st_transform(6933) |>
    dplyr::mutate(fd_area_m2 = as.numeric(sf::st_area(geometry)))
  projected_alert <- alert_geometry |>
    dplyr::select(area) |>
    sf::st_make_valid() |>
    sf::st_transform(6933)

  shares <- suppressWarnings(
    sf::st_intersection(projected_fd, projected_alert)
  ) |>
    dplyr::mutate(
      intersection_m2 = as.numeric(sf::st_area(geometry)),
      share = intersection_m2 / fd_area_m2
    ) |>
    sf::st_drop_geometry() |>
    dplyr::filter(share > 0) |>
    dplyr::arrange(prefecture, fd_id, dplyr::desc(share)) |>
    dplyr::mutate(
      share_rank = dplyr::row_number(),
      .by = c(prefecture, fd_id)
    )

  share_summary <- shares |>
    dplyr::summarise(
      area_weighted = area[which.max(share)],
      area_secondary = dplyr::nth(
        area[order(share, decreasing = TRUE)],
        2L,
        default = NA_character_
      ),
      share_primary = max(share),
      share_secondary = dplyr::nth(
        sort(share, decreasing = TRUE),
        2L,
        default = 0
      ),
      share_covered = sum(share),
      n_areas_intersected = dplyr::n(),
      n_areas_touched = sum(share >= straddle_share),
      .by = c(prefecture, fd_id, geometry_source)
    )

  result <- sample_departments |>
    dplyr::left_join(
      share_summary,
      by = dplyr::join_by(prefecture, fd_id)
    ) |>
    dplyr::left_join(
      station_assignment |>
        dplyr::select(!dplyr::any_of("stored_near_wbgt_station")),
      by = dplyr::join_by(prefecture, fd_id)
    ) |>
    dplyr::mutate(
      # Every sample department has department geometry since the crosswalk
      # repair, so this fill is a tripwire, not a fallback: the guard below
      # aborts if it ever fires.
      geometry_source = tidyr::replace_na(geometry_source, "none"),
      area_by_name = dplyr::if_else(
        prefecture %in% geometry_prefectures,
        NA_character_,
        zipangu::harmonize_prefecture_name(prefecture, to = "short")
      ),
      area_primary = dplyr::if_else(
        prefecture %in% geometry_prefectures,
        area_weighted,
        area_by_name
      ),
      agree = dplyr::if_else(
        is.na(area_by_name) | is.na(area_weighted),
        NA,
        area_by_name == area_weighted
      ),
      straddles = dplyr::if_else(
        is.na(n_areas_touched),
        NA,
        n_areas_touched > 1L
      )
    ) |>
    dplyr::select(
      prefecture,
      fd_id,
      area_primary,
      near_wbgt_station,
      dist_m,
      near_wbgt_station_any,
      dist_m_any,
      station_reassigned,
      station_years_covered,
      station_years_covered_any,
      share_primary,
      share_secondary,
      share_covered,
      n_areas_touched,
      n_areas_intersected,
      area_by_name,
      area_weighted,
      area_secondary,
      agree,
      straddles,
      geometry_source
    )

  if (
    anyNA(result$area_primary) ||
      !all(result$area_primary %in% alert_geometry$area)
  ) {
    rlang::abort(
      "Every sample unit must map to one of the 58 alert areas."
    )
  }
  n_straddles <- sum(result$straddles, na.rm = TRUE)
  if (n_straddles != expected_straddles) {
    rlang::abort(glue::glue(
      "Expected {expected_straddles} unit(s) to cross the ",
      "{straddle_share * 100}% area-share threshold, got {n_straddles}: ",
      "{paste(result$fd_id[result$straddles %in% TRUE], collapse = ', ')}."
    ))
  }
  if (any(result$agree == FALSE, na.rm = TRUE)) {
    rlang::abort(
      "A name-assigned unit's geometry disagrees with its prefecture name."
    )
  }
  if (!all(result$geometry_source == "department")) {
    rlang::abort(glue::glue(
      "Every sample unit should carry geometry; ",
      "{sum(result$geometry_source != 'department')} do not."
    ))
  }

  result
}


# Convert MOE issuance timestamps to the dates targeted by each bulletin.
# A 17:00 issuance targets the following day; a 05:00 issuance targets its date.
#
# `datetime` carries JST wall-clock values under a UTC tzone label, because the
# vendored reader builds it with as_datetime() from a naive string. Reading it
# with tz = "Asia/Tokyo" therefore treats 17:00 as 02:00 JST the next day and
# advances the date a second time, putting every 17:00 bulletin two days out.
# Take the date in the object's own clock and never convert.
make_moe_alert_days <- function(alert_records, valid_areas) {
  record_areas <- sort(unique(alert_records$area))
  valid_areas <- sort(unique(valid_areas))
  if (
    length(record_areas) != 58L ||
      length(valid_areas) != 58L ||
      !identical(record_areas, valid_areas)
  ) {
    rlang::abort(
      "The MOE record and spatial layer must contain the same 58 areas."
    )
  }

  issue_hour <- lubridate::hour(alert_records$datetime)
  if (!all(issue_hour %in% c(5L, 17L))) {
    rlang::abort("MOE issuance times must be 05:00 or 17:00.")
  }

  targeted <- alert_records |>
    dplyr::mutate(
      issue_hour = lubridate::hour(datetime),
      issue_date = lubridate::as_date(datetime),
      date = issue_date + dplyr::if_else(issue_hour == 17L, 1, 0)
    )

  # Pin the offset against the wall-clock date so that reintroducing a timezone
  # conversion fails here instead of silently shifting the 17:00 bulletins.
  expected_offset <- dplyr::if_else(targeted$issue_hour == 17L, 1, 0)
  if (
    !all(
      as.numeric(targeted$date - lubridate::as_date(targeted$datetime)) ==
        expected_offset
    )
  ) {
    rlang::abort(
      "A target date is not its issuance wall-clock date plus 0 (05:00) or 1 (17:00)."
    )
  }

  targeted |>
    dplyr::summarise(
      alert_day = any(alert),
      alert_5am = any(alert & issue_hour == 5L),
      alert_17pm = any(alert & issue_hour == 17L),
      .by = c(area, date)
    ) |>
    pointblank::rows_distinct(columns = c(area, date))
}
