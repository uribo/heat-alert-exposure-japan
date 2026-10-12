# Heat-alert exposure in Japan

This repository contains the reproducible exposure-side data pipeline for a study of heat alerts in Japan. Public service-area identifiers are anonymized as `FD001`–`FD728`. The pipeline reconstructs service-area geometry, assigns heat-alert areas and WBGT stations, computes population denominators, verifies frozen Ministry of the Environment inputs, and reconstructs the 2021–2022 forecast area-by-day table.

Restricted health outcome records are not included because their providers' data-use conditions do not permit redistribution. The `R/paper_c_*.R` files document the downstream analysis but are intentionally not sourced by the public pipeline and cannot run without those records. No restricted outcome totals are included at any geographic or temporal resolution.

## Rebuild

The project uses R 4.6.1, renv, and targets. Copy `Renviron.example` to `.Renviron`, then run:

```r
renv::restore()
targets::tar_make()
```

The build downloads only the versioned 2024 National Land Numerical Information administrative-area archive and the versioned 2019-01-25 Japan Meteorological Agency forecast-area archive. On the first fetch, the official archive bytes were pinned in `data-raw/source_sha256.tsv`; every build verifies those hashes, including when a local copy already exists. WBGT station coordinates come from the `jmastats` 0.3.0 package data and are joined to the shipped 2023 station roster. This version is pinned because a later package version drops coordinates for a station used by the study. All other inputs are shipped and checked locally. See `data/README.md` for file-level provenance and licensing.

The main acceptance objects are the `fd_alert_area_validation` and `forecast_area_day_validation` targets. Both abort on a mismatch. The former checks every identifier, alert area, station assignment and integer field of `fd_alert_area.csv` exactly, the three area shares within an absolute `1e-9` (measured `9e-12`), and the station distances within `1e-6` m (measured `1.8e-7` m, the centroid rounding of the polygon union). The latter checks all 16,592 rows exactly. The department polygons are rebuilt from the GeoJSON member of the N03 archive (the file the original analysis read); the shapefile member carries the same geometry at a different coordinate precision and moves every centroid by about 0.5 mm, which the distance check rejects.

## Sources and licenses

Code is released under the MIT License. Shipped data are released under CC BY 4.0 with source attribution, except crosswalk-derived tables, which are released under CC BY-SA 4.0. Upstream sources include MLIT, the Japan Meteorological Agency, the Ministry of the Environment, e-Stat, and the Japanese Wikipedia service-area list. File-level details are in `data/README.md`.

## Citation

Please cite the associated article and archived repository release:

> Uryu, S. (2026). *Heat-alert exposure in Japan*. Zenodo. DOI: **[TO BE ADDED]**.
