# Plate reader IO and growth curve plotting helpers.
# Sourced by the notebooks in scripts/.

# Plotting ----------------------------------------------------------------

# for plotting growth curves of a 96-well plate
plotplate <- function(
  df,
  dfxy,
  unsmoothed = TRUE,
  predicted = FALSE,
  plate,
  rows,
  cols,
  page,
  scales = "free_y"
) {
  dffilt <- dplyr::filter(df, plate_name == {{ plate }})
  xyfilt <- if (!is.null(dfxy)) {
    dplyr::inner_join(
      dfxy,
      dplyr::distinct(dffilt, well, plate_name),
      by = dplyr::join_by(well, plate_name)
    ) %>%
      tidyr::drop_na()
  }

  ggplot2::ggplot(dffilt, ggplot2::aes(x = hours)) +
    list(
      ggplot2::geom_line(ggplot2::aes(y = OD600_rollmean), color = "blue"),
      if (unsmoothed) {
        ggplot2::geom_line(ggplot2::aes(y = OD600), color = "orange", lty = 2)
      },
      if (predicted) {
        ggplot2::geom_line(ggplot2::aes(y = predicted), color = "orange")
      },
      if (!is.null(dfxy)) {
        ggplot2::geom_point(
          data = xyfilt,
          ggplot2::aes(x = x, y = y),
          color = "red",
          size = 2
        )
      },
      ggplot2::labs(x = "Hours", y = "OD600"),
      ggplot2::scale_x_continuous(
        breaks = seq(0, 48, 12),
        labels = seq(0, 48, 12)
      ),
      ggforce::facet_wrap_paginate(
        ~well,
        nrow = rows,
        ncol = cols,
        page = page,
        scales = scales
      ),
      ggplot2::theme(axis.text = ggplot2::element_text(size = 5))
    )
}


# Plate reader IO ---------------------------------------------------------

# Each plate reader keeps its own clock, so the raw elapsed time for what is
# meant to be the same nominal timepoint can drift by a second or more between
# plates/runs. Left as-is, that drift makes exact-value grouping (e.g.
# `group_by(hours)` when pooling replicates across plates) silently fracture
# into multiple near-duplicate groups instead of one, which shows up
# downstream as an artificially jagged mean. Snap elapsed hours onto the
# nominal read-interval grid (default 20 min / 1/3 hour, matching the design
# of these growth curve runs) so timepoints line up exactly across plates.
snap_hours <- function(hours, interval_hours = 1 / 3) {
  round(hours / interval_hours) * interval_hours
}

# Shared tidying for the wide "Time in column 1, wells across the top" layout
# that both plate readers export. Sets the first interval as time zero, pivots
# to long form, and zero-pads the well ids (A1 -> A01) so they join to the
# samplesheets.
tidy_plate_wide <- function(df) {
  df %>%
    # set interval start to be first cell and make all intervals relative to that
    # use time_length to just create an hours variable of type numeric
    dplyr::mutate(
      seconds = lubridate::time_length(
        lubridate::interval(Time[1], Time),
        unit = "second"
      )
    ) %>%
    tidyr::pivot_longer(
      c(-seconds, -Time),
      names_to = "well",
      values_to = "OD600"
    ) %>%
    # snap to the nominal read-interval grid so replicate plates align exactly
    dplyr::mutate(
      hours = snap_hours(lubridate::time_length(seconds, unit = "hours"))
    ) %>%
    # converting the well format so it matches the samplesheet
    dplyr::mutate(
      well = paste0(
        stringr::str_extract(well, "^[A-H]"),
        stringr::str_pad(
          stringr::str_extract(well, "\\d+"),
          width = 2,
          pad = "0",
          side = "left"
        )
      )
    ) %>%
    dplyr::select(seconds, hours, well, OD600) %>%
    dplyr::mutate(OD600 = as.numeric(OD600))
}

# Logphase 600 plate reader reading function. `path` is the full path to the
# xlsx export; a single export often holds several plates on different sheets.
read_logphase_xlsx <- function(path, sheet, skip) {
  readxl::read_xlsx(path, sheet = sheet, skip = skip) %>%
    tidy_plate_wide()
}

# Synergy H1 multimode plate reader reading function. Same wide layout as the
# Logphase 600 export, but the raw data sits on its own sheet.
read_synergy_xlsx <- function(path, sheet = 3, skip = 1) {
  readxl::read_xlsx(path, sheet = sheet, skip = skip) %>%
    tidy_plate_wide()
}

# Logphase 600 plate reader reading function but for the text file export.
# Unlike the xlsx export these files carry a header block that must be skipped,
# and the plate identity lives in the filename rather than in a sheet name.
read_logphase_txt <- function(file, skip) {
  readr::read_tsv(file, skip = skip) %>%
    dplyr::mutate(
      seconds = lubridate::time_length(
        lubridate::interval(Time[1], Time, tzone = "UTC"),
        unit = "second"
      ),
    ) %>%
    tidyr::pivot_longer(
      c(-seconds, -Time),
      names_to = "well",
      values_to = "OD600"
    ) %>%
    # snap to the nominal read-interval grid so replicate plates align exactly
    dplyr::mutate(
      hours = snap_hours(lubridate::time_length(seconds, unit = "hours"))
    ) %>%
    # converting the well format so it matches the samplesheet
    dplyr::mutate(
      well = paste0(
        stringr::str_extract(well, "^[A-H]"),
        stringr::str_pad(
          stringr::str_extract(well, "\\d+"),
          width = 2,
          pad = "0",
          side = "left"
        )
      )
    ) %>%
    dplyr::select(seconds, hours, well, OD600) %>%
    dplyr::mutate(OD600 = as.numeric(OD600)) %>%
    dplyr::mutate(plate_file = fs::path_file(file))
}


# Growth curve summaries --------------------------------------------------

# makes a summary table for inspecting 96-well growth curve plates
plate_summary <- function(df, plate_id) {
  df %>%
    dplyr::filter(plate_name == {{ plate_id }}) %>%
    dplyr::select(
      plate_name,
      column,
      evo_hist,
      strainID,
      carbon_source,
      streptomycin_ug_ml,
      plate_file
    ) %>%
    dplyr::distinct()
}

# function for calculating AUC using trapezoid rule
trap_auc <- function(x, y) {
  sum(diff(x) * (utils::head(y, -1) + utils::tail(y, -1)) / 2)
}

# Function for calculating bootstrapped mean and 95% CI
boot_mean_ci <- function(x, n_boot = N_BOOT, ci = 0.95) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(tibble::tibble(mn = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_))
  }
  alpha <- (1 - ci) / 2
  boot_means <- replicate(n_boot, mean(sample(x, length(x), replace = TRUE)))
  tibble::tibble(
    mn = mean(x),
    ci_lo = stats::quantile(boot_means, alpha),
    ci_hi = stats::quantile(boot_means, 1 - alpha)
  )
}

# To estimate a 95% confidence interval, the formula generally takes
# the form of: Estimate ~(1.96 * std.error)
dirty_mean_ci <- function(x, ci = 0.95) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(tibble::tibble(mn = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_))
  }
  y <- mean(x)
  se <- sd(x) / sqrt(length(x))
  tibble::tibble(
    mn = y,
    ci_lo = y - 1.96 * se,
    ci_hi = y + 1.96 * se
  )
}

tidy_bayes_boot <- function(
  x,
  estimator,
  use.weights = TRUE,
  n_boot = 4000,
  ci = 0.95
) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(tibble::tibble(mn = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_))
  }
  alpha <- (1 - ci) / 2

  boot_estimator <- pull(bayesboot::bayesboot(
    x,
    {{ estimator }},
    R = n_boot,
    R1 = n_boot,
    use.weights = use.weights
  ))

  tibble::tibble(
    md = stats::quantile(boot_estimator, 0.5),
    ci_lo = stats::quantile(boot_estimator, alpha),
    ci_hi = stats::quantile(boot_estimator, 1 - alpha)
  )
}
