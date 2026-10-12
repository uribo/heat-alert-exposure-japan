# Client for the MOE forecast WBGT archive (getForecastData), gate 10.
#
# Spec: https://www.wbgt.env.go.jp/man15NH/wbgt_data_api_service_manual.pdf
# (API 仕様書 1.1, 2026-06-24). Facts established by probing on 2026-09-29
# (notes/paper_strategy.md 6.7):
#   - A bulletin is issued every hour. Each carries 3-hourly values from the
#     next 3-hour slot to 24:00 two days ahead.
#   - forecast_val is a string holding ten times the WBGT, and every value
#     seen is a multiple of 10: forecasts are rounded to 1 degree C.
#   - fixed_time = "HHMMSS" restricts a date_search_type = 1 query to the
#     bulletins issued at that clock time; the result equals the matching
#     subset of the hourly pull.
#   - One call returns at most 25,000 rows. Past that the API returns
#     status "error" with no data, so a query is split in two only on that
#     message and every other error aborts.
#   - The manual defines `flag` for getSurveyData only. Every probed value
#     was 0; any other value aborts the read until it is understood.
#
# Pulls are by prefecture code (location_type = 2) rather than by station,
# so that the station roster is the one the API held at the time, which is
# the roster the alert rule was applied to.

moe_forecast_endpoint <- "https://www.wbgt.env.go.jp/api/v1/getForecastData"

# Substring of the errMsg returned when a query exceeds 25,000 rows.
moe_forecast_cap_message <- "検索結果が閾値以上"

# A query by station list (`stations = TRUE`) is named by the prefecture it
# stands in for, so that its files cannot be mistaken for a prefecture pull.
moe_forecast_file_name <- function(
  pref_cd,
  from_date,
  to_date,
  fixed_time,
  stations = FALSE
) {
  glue::glue(
    "forecast_{if (stations) 'stations' else 'pref'}{pref_cd}_",
    "{format(from_date, '%Y%m%d')}_",
    "{format(to_date, '%Y%m%d')}_{fixed_time}.json.gz"
  )
}

# Bulletins issued at `fixed_time` on every date from `from_date` to
# `to_date` (inclusive, JST), for all stations of one prefecture code, or for
# the stations in `wbgt_nos` (location_type = 1, one wbgt_nos parameter per
# station; a comma-separated list is rejected). Okinawa needs the latter: no
# prefecture code the API accepts covers its stations (Issue #59, 2026-10-05).
moe_forecast_url <- function(
  pref_cd,
  from_date,
  to_date,
  fixed_time,
  wbgt_nos = NULL
) {
  location <- if (is.null(wbgt_nos)) {
    c(location_type = "2", date_search_type = "1", pref_cds = pref_cd)
  } else {
    c(
      location_type = "1",
      date_search_type = "1",
      stats::setNames(wbgt_nos, rep("wbgt_nos", length(wbgt_nos)))
    )
  }
  query <- c(
    location,
    range_date_from = paste0(format(from_date, "%Y%m%d"), "000000"),
    range_date_to = paste0(format(to_date, "%Y%m%d"), "235959"),
    fixed_time = fixed_time
  )
  paste0(
    moe_forecast_endpoint,
    "?",
    paste(names(query), query, sep = "=", collapse = "&")
  )
}

# Fetch one query and freeze the response body, as received, gzipped under
# `dir`. The sha256 is taken over the decompressed body. Returns one manifest
# row per file written; a query the API cannot serve whole (row cap, repeated
# 5xx) is split into halves until each part succeeds or one day fails. The
# split is a loop over pending ranges, not recursion: targets reads R/ into
# its dependency graph, and recursion between functions there is a cycle.
# `on_row`, when given, is called with each manifest row as soon as its file
# is written, so that a caller can record the parts of a split query that
# succeeded before a later part fails.
fetch_moe_forecast <- function(
  pref_cd,
  from_date,
  to_date,
  fixed_time,
  dir,
  pause = 1,
  on_row = NULL,
  wbgt_nos = NULL
) {
  pending <- list(c(from_date, to_date))
  rows <- list()
  while (length(pending) > 0L) {
    range <- pending[[1L]]
    pending <- pending[-1L]
    row <- fetch_moe_forecast_range(
      pref_cd,
      range[1L],
      range[2L],
      fixed_time,
      dir,
      pause,
      wbgt_nos = wbgt_nos
    )
    if (is.null(row)) {
      mid <- range[1L] + floor(as.numeric(range[2L] - range[1L]) / 2)
      pending <- c(list(c(range[1L], mid), c(mid + 1, range[2L])), pending)
    } else {
      rows <- c(rows, list(row))
      if (!is.null(on_row)) {
        on_row(row)
      }
    }
  }
  dplyr::bind_rows(rows)
}

# One request. Returns a manifest row, or NULL when the range must be split.
fetch_moe_forecast_range <- function(
  pref_cd,
  from_date,
  to_date,
  fixed_time,
  dir,
  pause,
  wbgt_nos = NULL
) {
  url <- moe_forecast_url(pref_cd, from_date, to_date, fixed_time, wbgt_nos)
  handle <- curl::new_handle(
    useragent = "uribo/23_Heatstroke forecast-archive fetcher (R curl)"
  )
  # A 5xx is retried twice with a growing wait: on 2026-09-29 one request
  # returned 502 and the same request succeeded in 2.3 s on retry. Queries
  # near the row cap (about 23,500 rows) returned 502 three times running,
  # so a third 5xx splits the range like the cap does; a single day that
  # still fails aborts. Any other non-200 status aborts.
  for (attempt in 1:3) {
    Sys.sleep(pause * attempt^2)
    res <- curl::curl_fetch_memory(url, handle = handle)
    if (res$status_code < 500L) {
      break
    }
    message("HTTP ", res$status_code, " (attempt ", attempt, ") for ", url)
  }
  if (res$status_code >= 500L && from_date < to_date) {
    message("Splitting after repeated HTTP ", res$status_code, ": ", url)
    return(NULL)
  }
  if (res$status_code != 200L) {
    rlang::abort(glue::glue("HTTP {res$status_code} for {url}"))
  }
  body <- res$content
  parsed <- jsonlite::fromJSON(rawToChar(body), simplifyVector = FALSE)

  if (!identical(parsed$status, "success")) {
    msg <- parsed$errMsg %||% "(no errMsg)"
    if (!grepl(moe_forecast_cap_message, msg, fixed = TRUE)) {
      rlang::abort(glue::glue("API error for {url}: {msg}"))
    }
    if (from_date == to_date) {
      rlang::abort(glue::glue("One day exceeds the row cap: {url}"))
    }
    return(NULL)
  }

  n_data <- length(parsed$data)
  if (!identical(as.integer(parsed$count), n_data) || n_data == 0L) {
    rlang::abort(glue::glue(
      "count {parsed$count} and {n_data} data rows for {url}"
    ))
  }

  file <- moe_forecast_file_name(
    pref_cd,
    from_date,
    to_date,
    fixed_time,
    stations = !is.null(wbgt_nos)
  )
  tmp <- tempfile(tmpdir = dir, fileext = ".tmp")
  con <- gzfile(tmp, "wb")
  writeBin(body, con)
  close(con)
  path <- file.path(dir, file)
  if (!file.rename(tmp, path)) {
    rlang::abort(glue::glue("Could not move the response to {path}."))
  }
  con <- gzfile(path, "rb")
  written <- readBin(con, "raw", n = length(body) + 1L)
  close(con)
  if (!identical(written, body)) {
    rlang::abort(glue::glue("{path} does not read back as the response body."))
  }

  tibble::tibble(
    file = file,
    pref_cd = pref_cd,
    from_date = from_date,
    to_date = to_date,
    fixed_time = fixed_time,
    url = url,
    fetched_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S JST", tz = "Asia/Tokyo"),
    records = n_data,
    bytes = length(body),
    sha256 = digest::digest(body, algo = "sha256", serialize = FALSE)
  )
}

# Read frozen files listed in `manifest`, checking each body against its
# sha256 and row count. forecast_val stays in tenths of a degree (integer).
read_moe_forecast <- function(manifest, dir) {
  purrr::pmap(
    dplyr::select(manifest, file, bytes, records, sha256, fixed_time),
    \(file, bytes, records, sha256, fixed_time) {
      con <- gzfile(file.path(dir, file), "rb")
      body <- readBin(con, "raw", n = bytes + 1L)
      close(con)
      got <- digest::digest(body, algo = "sha256", serialize = FALSE)
      if (length(body) != bytes || got != sha256) {
        rlang::abort(glue::glue(
          "{file}: sha256 or size differs from the manifest."
        ))
      }
      parsed <- jsonlite::fromJSON(rawToChar(body))
      df <- tibble::as_tibble(parsed$data)
      if (nrow(df) != records) {
        rlang::abort(glue::glue("{file}: {nrow(df)} rows, manifest {records}."))
      }
      # The API has to have honoured fixed_time: every bulletin in the file
      # is issued at that clock time.
      issued <- unique(substr(df$reference_time, 12L, 19L))
      expected <- paste(
        substring(fixed_time, c(1L, 3L, 5L), c(2L, 4L, 6L)),
        collapse = ":"
      )
      if (!identical(issued, expected)) {
        rlang::abort(glue::glue(
          "{file}: bulletins issued at {paste(issued, collapse = ', ')}, ",
          "expected {expected} only."
        ))
      }
      df
    }
  ) |>
    purrr::list_rbind() |>
    dplyr::transmute(
      wbgt_no = as.character(wbgt_no),
      reference_time = lubridate::ymd_hms(reference_time, tz = "Asia/Tokyo"),
      forecast_time = lubridate::ymd_hms(forecast_time, tz = "Asia/Tokyo"),
      forecast_val = as.integer(forecast_val),
      flag = as.integer(flag)
    ) |>
    moe_forecast_check()
}

moe_forecast_check <- function(df) {
  if (anyNA(df)) {
    rlang::abort("Forecast rows with missing fields.")
  }
  if (any(df$flag != 0L)) {
    flags <- table(df$flag)
    rlang::abort(glue::glue(
      "Non-zero flag values ({paste(names(flags), flags, sep = ': ', collapse = ', ')}); ",
      "getForecastData does not define them."
    ))
  }
  pointblank::rows_distinct(
    df,
    columns = c(wbgt_no, reference_time, forecast_time)
  )
}

# Daily maximum forecast WBGT per station, bulletin and target date. An
# afternoon bulletin (17:00, or 16:00 in the sensitivity analysis of 6.9.1)
# targets the following day and a morning one (05:00, or 04:00) its own day,
# as in make_moe_alert_days(). A day covers the 3-hourly slots 03:00 to
# 24:00, so a morning bulletin has 7 slots (06:00 to 24:00) and an afternoon
# one 8.
moe_forecast_daily_max <- function(forecast) {
  issue_hour <- lubridate::hour(forecast$reference_time)
  if (!all(issue_hour %in% c(4L, 5L, 16L, 17L))) {
    rlang::abort(
      "Only the 04:00, 05:00, 16:00 and 17:00 bulletins are expected."
    )
  }
  daily <- forecast |>
    dplyr::mutate(
      issue_hour = lubridate::hour(reference_time),
      target_date = lubridate::as_date(reference_time) +
        dplyr::if_else(issue_hour >= 12L, 1L, 0L),
      slot_date = lubridate::as_date(forecast_time - 1)
    ) |>
    dplyr::filter(slot_date == target_date) |>
    dplyr::summarise(
      n_slots = dplyr::n(),
      forecast_max = max(forecast_val) / 10,
      .by = c(wbgt_no, reference_time, issue_hour, target_date)
    )
  short <- dplyr::filter(
    daily,
    n_slots != dplyr::if_else(issue_hour >= 12L, 8L, 7L)
  )
  if (nrow(short) > 0L) {
    rlang::abort(glue::glue(
      "{nrow(short)} station-bulletins lack 3-hourly slots for the target day, ",
      "e.g. station {short$wbgt_no[1]} at {short$reference_time[1]}."
    ))
  }
  daily
}

# SHA256SUMS over the frozen .json.gz files themselves, in `shasum -a 256`
# format, so that any copy of the directory can be checked with
# `shasum -a 256 -c` without R. The manifest's sha256 is over the decompressed
# body and is what read_moe_forecast() verifies.
write_moe_forecast_sums <- function(manifest, dir, path) {
  files <- sort(manifest$file)
  sums <- vapply(
    file.path(dir, files),
    \(x) digest::digest(x, algo = "sha256", file = TRUE),
    character(1)
  )
  writeLines(paste0(sums, "  ", files), path)
  invisible(path)
}
