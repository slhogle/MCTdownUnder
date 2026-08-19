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
    dplyr::mutate(hours = lubridate::time_length(seconds, unit = "hours")) %>%
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
