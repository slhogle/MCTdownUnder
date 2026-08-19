library(here)
library(tidyverse)
library(readxl)
library(stringr)
library(lubridate)
library(fs)
library(slider)
source(here::here("R", "utils_gcurves.R"))

# Global vars
data_raw <- here::here(
  "data_raw",
  "serial_transfer_coculture",
  "20260410_serial_transfer"
)
data <- here::here(
  "data",
  "serial_transfer_coculture",
  "20260410_serial_transfer"
)

# make processed data and figs directories if they don't exist
fs::dir_create(data)

# Read data --------------------------------------------------------------

# Read sample metadata
samplesheet_coculture <- readxl::read_xlsx(here::here(
  data_raw,
  "coculture_samplesheet.xlsx"
))

# all dilutions done at 10,000 fold

# Cycle 1 ----------------------------------------------------------------

# 1977
c1_1977 <- read_logphase_xlsx(
  "",
  "Serial_transfer_1_start_10-huhti-2026 12-14-59.xlsx",
  2,
  1
) %>%
  right_join(
    filter(samplesheet_coculture, strain == "H1977"),
    by = join_by(well)
  ) %>%
  mutate(dilution = 10000, plate_name = "c1", start_hrs_from_0 = 0)

# 1287
c1_1287 <- read_logphase_xlsx(
  "",
  "Serial_transfer_1_start_10-huhti-2026 12-14-59.xlsx",
  5,
  1
) %>%
  right_join(
    filter(samplesheet_coculture, strain == "H1287"),
    by = join_by(well)
  ) %>%
  mutate(dilution = 10000, plate_name = "c1", start_hrs_from_0 = 0)

# Cycle 2 ----------------------------------------------------------------

# 1977
c2_1977 <- read_logphase_xlsx(
  "",
  "Serial_transfer_2_10.4_12-huhti-2026 13-12-32.xlsx",
  2,
  1
) %>%
  right_join(
    filter(samplesheet_coculture, strain == "H1977"),
    by = join_by(well)
  ) %>%
  mutate(dilution = 10000, plate_name = "c2", start_hrs_from_0 = 48)

# 1287
c2_1287 <- read_logphase_xlsx(
  "",
  "Serial_transfer_2_10.4_12-huhti-2026 13-12-32.xlsx",
  5,
  1
) %>%
  right_join(
    filter(samplesheet_coculture, strain == "H1287"),
    by = join_by(well)
  ) %>%
  mutate(dilution = 10000, plate_name = "c2", start_hrs_from_0 = 48)

# Cycle 3 ----------------------------------------------------------------

# 1977
c3_1977 <- read_logphase_xlsx(
  "",
  "Serial_transfer_3_12.4_14-huhti-2026 13-05-36.xlsx",
  2,
  1
) %>%
  right_join(
    filter(samplesheet_coculture, strain == "H1977"),
    by = join_by(well)
  ) %>%
  mutate(dilution = 10000, plate_name = "c3", start_hrs_from_0 = 96)

# 1287
c3_1287 <- read_logphase_xlsx(
  "",
  "Serial_transfer_3_12.4_14-huhti-2026 13-05-36.xlsx",
  5,
  1
) %>%
  right_join(
    filter(samplesheet_coculture, strain == "H1287"),
    by = join_by(well)
  ) %>%
  mutate(dilution = 10000, plate_name = "c3", start_hrs_from_0 = 96)

# Cycle 4 ----------------------------------------------------------------

# 1977
c4_1977 <- read_logphase_xlsx(
  "",
  "Serial_transfer_4_16-huhti-2026 13-33-18.xlsx",
  2,
  1
) %>%
  right_join(
    filter(samplesheet_coculture, strain == "H1977"),
    by = join_by(well)
  ) %>%
  mutate(dilution = 10000, plate_name = "c4", start_hrs_from_0 = 144)

# 1287
c4_1287 <- read_logphase_xlsx(
  "",
  "Serial_transfer_4_16-huhti-2026 13-33-18.xlsx",
  5,
  1
) %>%
  right_join(
    filter(samplesheet_coculture, strain == "H1287"),
    by = join_by(well)
  ) %>%
  mutate(dilution = 10000, plate_name = "c4", start_hrs_from_0 = 144)

# Cycle 5 ----------------------------------------------------------------

# 1977
c5_1977 <- read_logphase_xlsx(
  "",
  "Serial_transfer_5_20-huhti-2026 07-42-52.xlsx",
  2,
  1
) %>%
  right_join(
    filter(samplesheet_coculture, strain == "H1977"),
    by = join_by(well)
  ) %>%
  mutate(dilution = 10000, plate_name = "c5", start_hrs_from_0 = 192)

# 1287
c5_1287 <- read_logphase_xlsx(
  "",
  "Serial_transfer_5_20-huhti-2026 07-42-52.xlsx",
  5,
  1
) %>%
  right_join(
    filter(samplesheet_coculture, strain == "H1287"),
    by = join_by(well)
  ) %>%
  mutate(dilution = 10000, plate_name = "c5", start_hrs_from_0 = 192)

# Combine, tidy, format --------------------------------------------------

# combine all samples, group by plate + well, calculate rolling mean. This
# reduces the jaggedness of the curves and also makes the resulting saved data
# file smaller

# 1977
coculture_gcurves_sm_1977 <- bind_rows(
  c1_1977,
  c2_1977,
  c3_1977,
  c4_1977,
  c5_1977
) %>%
  dplyr::group_by(plate_name, well) %>%
  arrange(plate_name, well) %>%
  dplyr::mutate(
    OD600_rollmean = slider::slide_dbl(OD600, mean, .before = 2, .after = 2)
  ) %>%
  ungroup() %>%
  relocate(OD600_rollmean, .after = "OD600")

# 1287
coculture_gcurves_sm_1287 <- bind_rows(
  c1_1287,
  c2_1287,
  c3_1287,
  c4_1287,
  c5_1287
) %>%
  dplyr::group_by(plate_name, well) %>%
  arrange(plate_name, well) %>%
  dplyr::mutate(
    OD600_rollmean = slider::slide_dbl(OD600, mean, .before = 2, .after = 2)
  ) %>%
  ungroup() %>%
  relocate(OD600_rollmean, .after = "OD600")

# save result for later
readr::write_tsv(
  coculture_gcurves_sm_1977,
  here::here(data, "1977_coculture_gcurves_smooth.tsv")
)

# save result for later
readr::write_tsv(
  coculture_gcurves_sm_1287,
  here::here(data, "1287_coculture_gcurves_smooth.tsv")
)
