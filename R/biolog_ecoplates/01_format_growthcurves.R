library(here)
library(tidyverse)
library(readxl)
library(stringr)
library(lubridate)
library(fs)
library(ggforce)
library(slider)
source(here::here("R", "utils_gcurves.R"))

data_raw <- here::here("_data_raw", "biolog_ecoplates")
data <- here::here("data", "biolog_ecoplates")
figs <- here::here("figs", "biolog_ecoplates")

# make processed data and figs directories if they don't exist
fs::dir_create(data)
fs::dir_create(figs)

# Read and tidy Synergy output -------------------------------------------

# Tidying growth curves from Synergy H1 multimode microplate reader

plate01 <- readxl::read_xlsx(
  here::here(data_raw, "biolog_ecoplate_1287A", "1287A_raw_data.xlsx"),
  sheet = 3,
  skip = 1
) %>%
  # set interval start to be first cell and make all intervals relative to that
  # use time_length to just create an hours variable of type numeric
  mutate(
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
  mutate(hours = lubridate::time_length(seconds, unit = "hours")) %>%
  # converting the well format so it matches the samplesheet
  mutate(
    well = paste0(
      str_extract(well, "^[A-H]"),
      str_pad(str_extract(well, "\\d+"), width = 2, pad = "0", side = "left")
    )
  ) %>%
  dplyr::select(seconds, hours, well, OD600) %>%
  # create a plate variable for later combining
  mutate(plate_name = "plate01")

plate02 <- readxl::read_xlsx(
  here::here(
    data_raw,
    "biolog_ecoplate_1287E_1977A",
    "1287E_1977A_rawdata.xlsx"
  ),
  sheet = 3,
  skip = 1
) %>%
  # set interval start to be first cell and make all intervals relative to that
  # use time_length to just create an hours variable of type numeric
  mutate(
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
  mutate(hours = lubridate::time_length(seconds, unit = "hours")) %>%
  # converting the well format so it matches the samplesheet
  mutate(
    well = paste0(
      str_extract(well, "^[A-H]"),
      str_pad(str_extract(well, "\\d+"), width = 2, pad = "0", side = "left")
    )
  ) %>%
  dplyr::select(seconds, hours, well, OD600) %>%
  # create a plate variable for later combining
  mutate(plate_name = "plate02")

plate03 <- readxl::read_xlsx(
  here::here(
    data_raw,
    "biolog_ecoplate_1977A_1977E",
    "1977A_1977E_rawdata.xlsx"
  ),
  sheet = 3,
  skip = 1
) %>%
  # set interval start to be first cell and make all intervals relative to that
  # use time_length to just create an hours variable of type numeric
  mutate(
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
  mutate(hours = lubridate::time_length(seconds, unit = "hours")) %>%
  # converting the well format so it matches the samplesheet
  mutate(
    well = paste0(
      str_extract(well, "^[A-H]"),
      str_pad(str_extract(well, "\\d+"), width = 2, pad = "0", side = "left")
    )
  ) %>%
  dplyr::select(seconds, hours, well, OD600) %>%
  # create a plate variable for later combining
  mutate(plate_name = "plate03")


# Combine, tidy, format --------------------------------------------------

## Read sample metadata

samplesheet01 <- readxl::read_xlsx(here::here(
  data_raw,
  "biolog_ecoplate_1287A",
  "samplesheet_1287A.xlsx"
)) %>%
  mutate(strainID = paste0("HAMBI_", strain)) %>%
  mutate(plate_name = "plate01")

samplesheet02 <- readxl::read_xlsx(here::here(
  data_raw,
  "biolog_ecoplate_1287E_1977A",
  "samplesheet_1287E_1977A.xlsx"
)) %>%
  mutate(strainID = paste0("HAMBI_", strain)) %>%
  mutate(plate_name = "plate02")

samplesheet03 <- readxl::read_xlsx(here::here(
  data_raw,
  "biolog_ecoplate_1977A_1977E",
  "samplesheet_1977A_1977E.xlsx"
)) %>%
  mutate(strainID = paste0("HAMBI_", strain)) %>%
  mutate(plate_name = "plate03")

## Join with metadata to remove ununsed samples

ecoplate_gcurves_sm <- bind_rows(plate01, plate02, plate03) %>%
  left_join(
    bind_rows(samplesheet01, samplesheet02, samplesheet03),
    by = join_by(well, plate_name)
  ) %>%
  dplyr::group_by(plate_name, well) %>%
  dplyr::mutate(
    OD600_rollmean = slider::slide_dbl(OD600, mean, .before = 2, .after = 2)
  ) %>%
  ungroup() %>%
  relocate(OD600_rollmean, .after = "OD600")

# save output for later use
readr::write_tsv(
  ecoplate_gcurves_sm,
  here::here(data, "gcurves_smooth.tsv")
)

# Plot growth curves for later inspection --------------------------------

# quick function for saving files

plotter <- function(id) {
  ggsave(
    here::here(figs, paste0(id, ".svg")),
    plotplate(
      ecoplate_gcurves_sm,
      dfxy = NULL,
      unsmoothed = TRUE,
      predicted = FALSE,
      plate = id,
      rows = 8,
      cols = 12,
      page = 1,
      scales = "fixed"
    ),
    device = "svg",
    width = 12,
    height = 9
  )
}

paste0("plate", str_pad(seq(1:3), 2, "left", pad = "0")) %>%
  map(plotter)
