# NOTE: readable analysis code from the private project. It needs the fire-department
# transport records, which are not distributed, and is not sourced by the public
# pipeline (_targets.R). Kept so that the published estimates can be audited.

# Geometry for paper C's map (Issue #43 item 15): the 58 Heatstroke Alert
# forecast areas, simplified, with a flag for the areas holding sample
# departments. Written as a GeoJSON that the manuscript reads without the
# store, so it is a file target: a change in the areas, the assignment or the
# panel rebuilds it. No department geometry or name enters the file.
write_paper_c_alert_areas <- function(areas, assignment, panel, path) {
  departments <- unique(panel$department)
  sample_areas <- assignment |>
    dplyr::filter(department %in% departments) |>
    dplyr::summarise(
      prefectures = paste(sort(unique(prefecture)), collapse = "; "),
      .by = area_primary
    )
  counts <- panel |>
    dplyr::distinct(area_primary, department) |>
    dplyr::count(area_primary, name = "departments")
  sample_areas <- sample_areas |>
    dplyr::inner_join(counts, by = dplyr::join_by(area_primary))
  if (nrow(sample_areas) != 34L || sum(sample_areas$departments) != 417L) {
    rlang::abort(glue::glue(
      "Expected 34 sample areas and 417 departments, got ",
      "{nrow(sample_areas)} and {sum(sample_areas$departments)}."
    ))
  }
  out <- areas |>
    sf::st_transform(6668) |>
    sf::st_simplify(dTolerance = 500) |>
    sf::st_transform(4326) |>
    dplyr::left_join(sample_areas, by = dplyr::join_by(area == area_primary)) |>
    dplyr::mutate(in_sample = !is.na(departments)) |>
    dplyr::select(code, area, in_sample, departments, prefectures)
  if (
    nrow(out) != 58L || sum(out$in_sample) != 34L || !all(sf::st_is_valid(out))
  ) {
    rlang::abort(
      "The simplified alert areas are not 58 valid polygons with 30 in the sample."
    )
  }
  if (file.exists(path)) {
    unlink(path)
  }
  sf::st_write(out, path, quiet = TRUE)
  path
}
