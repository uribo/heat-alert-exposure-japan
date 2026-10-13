# Data sources, provenance, and licenses

All identifiers beginning with `FD` are anonymized. The public tables contain no service-organization names and no restricted health outcome records.

## Release tables

| File | Contents and source | License | SHA256 |
|---|---|---|---|
| `data-raw/release/alert_area_roster.csv` | Heat-alert areas in the study sample, derived from Ministry of the Environment alert records | CC BY 4.0 | `3737beea490f39b538253462452fb828bd7d1b3a2a2b2bbdf0b1976ac0098fe4` |
| `data-raw/release/fd_alert_area.csv` | Anonymized service-area to alert-area and WBGT-station assignment; joins the crosswalk below to official spatial sources | CC BY-SA 4.0 | `efd2306477a5d1a9c8f5c53f4f41ce630ca13c3bad02e64d12748b0c344d96e7` |
| `data-raw/release/fd_municipality.csv` | Anonymized service-area to municipality-code crosswalk derived from the Japanese Wikipedia service-area list | CC BY-SA 4.0 | `9ac3e7b1457c6ee0185be1dfad3eaec3f161407ad2901c15a116d77abeaf95a9` |
| `data-raw/release/fd_population_exclusions.csv` | Explicit 2020 crosswalk exclusions derived from the crosswalk audit | CC BY-SA 4.0 | `805ce4a4654ab472028628038bedbfd180efdfb20ffe70672e15fe18d1eeb305` |
| `data-raw/release/fd_shared_municipality_weights.csv` | Weights for shared municipalities, derived from the crosswalk and official historical populations | CC BY-SA 4.0 | `1b16919934161fcacd5f1829db2897acb1f5e345a45936bbbadda16980040645` |
| `data-raw/release/moe_alert_historical.csv` | Historical heat-alert bulletins scraped from Ministry of the Environment pages | CC BY 4.0 | `ccaddb4298ff4f6ab171dc31893f2818f78df1f6375c33f01fa8106a572410fa` |
| `data-raw/release/moe_forecast_area_day.csv` | Area-by-day reconstruction derived from the frozen Ministry of the Environment forecast bodies | CC BY 4.0 | `1148bb89a09a2e90c81ab4ac488a00a9126b49566159b0667e34f34f19ffe791` |
| `data-raw/release/moe_wbgt_daily_agg.csv.gz` | Station-by-day WBGT summary derived from Ministry of the Environment final station-month files | CC BY 4.0 | `d870916389aa98d5c727152d2ad1a941ad04c58424da22b3d5166ee6ff616db6` |
| `data-raw/release/population_2020.csv` | 2020 municipality population and age classes from e-Stat census table 0003445139, mapped to 2024 N03 codes | CC BY 4.0 | `b5e2fc7c1aa440d4b5bca7cbf07a96a31abc9e0d54a25a930a65a74122ffa30c` |
| `data-raw/release/sample_departments.csv` | Anonymized study sample derived from the service-area crosswalk | CC BY-SA 4.0 | `5f4fc90d3523f3a3e7db08aa9a4fb1f4e87d1cec352c77ee1d43bbd6617bb0ce` |

**Export provenance.** These ten files were built on the canonical targets store (suMS, 2026-10-13, private repository at `be81bc8`) from an up-to-date upstream (only the 14 `release_*` targets were outdated; no `shortcut`). They are byte-identical to the provisional set built earlier on the author's laptop, so the SHA256 values above are the final-release values. The department-name scan (`src/department_name_scan.R` in the private repository) exits 0 on the set; the only stem match is the place name 佐賀 (Saga), which is also an alert-area name.

## Other shipped inputs

| File or directory | Source and role | License and integrity |
|---|---|---|
| `data-raw/frozen/estat_0000020203_taxable_income.csv` | e-Stat taxable-income table 0000020203 | CC BY 4.0; checked by `data-raw/frozen/SHA256SUMS` |
| `data-raw/frozen/moe_forecast/*.json.gz` | Ministry of the Environment forecast API responses used for the 05:00/17:00 and 04:00/16:00 reconstructions | CC BY 4.0; 1,217 files checked by the two shipped SHA256SUMS files and documented by the two manifests, whose `count` column is the number of forecast records in each file (the manifests are kept byte-identical to the private provenance record) |
| `data-raw/wbgt_stations2023.csv` | Ministry of the Environment 2023 WBGT station roster | CC BY 4.0 |

WBGT station coordinates are package data from `jmastats` 0.3.0; the pipeline de-duplicates them by area and station number, applies the Choshi station-number correction used to build the released table, and joins them to the shipped 2023 roster. Version 0.3.0 is pinned because a later package version drops coordinates for a station used by the study. No AMeDAS metadata file is downloaded.

For `fd_alert_area.csv`, the rebuild reproduces every identifier, area, station and integer field exactly; `share_primary`, `share_secondary` and `share_covered` differ from the released values by at most `1e-9` in absolute terms (measured `9e-12`), and `dist_m` / `dist_m_any` by at most `1e-6` m (measured `1.8e-7` m). The department polygons are rebuilt from `N03-20240101.geojson` inside the N03 archive, the member the original analysis read; station coordinates keep the WGS 84 coordinates of `jmastats::stations`, in which the released distances were measured.

The MOE station-month source can be reconstructed with `R/read_moe_wbgt.R`, `R/wbgt.R`, and `R/moe_wbgt_api.R`. Historical alert pages can be reconstructed with `R/moe_alert.R`. These live-source rebuilds are not run by default because the repository ships the frozen derived tables used by the study.

## Downloaded spatial sources

The targets pipeline downloads and verifies these sources:

- MLIT National Land Numerical Information N03 administrative areas, 2024 (`N03-20240101_GML.zip`). Municipality joins use `N03_007` directly; ward codes are not collapsed.
- JMA forecast-area GIS, published 2019-01-25 (`20190125_AreaForecastLocalM_prefecture_GIS.zip`).

Both archives are versioned official sources. Their bytes were pinned at first fetch, and `data-raw/source_sha256.tsv` records the verified SHA256 values. The `n03_archive` and `jma_alert_archive` targets skip downloading when those files already exist but still verify their bytes. The `fd_alert_area_validation` and `forecast_area_day_validation` targets abort if reconstructed outputs differ from the shipped comparison tables.
