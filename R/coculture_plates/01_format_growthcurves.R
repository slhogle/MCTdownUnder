library(here)
library(tidyverse)
library(readxl)
library(stringr)
library(lubridate)
library(fs)
library(ggforce)
library(slider)
source(here::here("R", "utils_gcurves.R"))

# Global vars
data_raw <- here::here("_data_raw", "coculture_plates")
data <- here::here("data", "coculture_plates")
figs <- here::here("figs", "coculture_plates")

# make processed data and figs directories if they don't exist
fs::dir_create(data)
fs::dir_create(figs)

# Batch01 - 2025.09.20 ---------------------------------------------------

batch <- "20250920_batch01"

# Read sample metadata
samplesheet_diff_test <- readxl::read_xlsx(here::here(
  data_raw,
  batch,
  "diffusion_test_2_samplesheet.xlsx"
))
samplesheet_coculture <- readxl::read_xlsx(here::here(
  data_raw,
  batch,
  "coculture_samplesheet.xlsx"
))

# Read crowth curves
# Diffusion test plate01
plate01 <- read_logphase_xlsx(
  batch,
  "diffusion_test_2_Co-culture_1977_1287AE_12.9_15-syys-2025 09-00-26.xlsx",
  2,
  1
) %>%
  # remove rows A and H
  filter(str_detect(well, "^A|^H", negate = TRUE)) %>%
  left_join(samplesheet_diff_test, by = join_by(well)) %>%
  mutate(plate_name = "plate01")

# Experiment plate02 - 0 ug/ml streptomycin
plate02 <- read_logphase_xlsx(
  batch,
  "Co-culture_0_4_8ug_plates_22-syys-2025 08-17-10.xlsx",
  2,
  1
) %>%
  # rows A and H not collected here so they are NA
  drop_na() %>%
  left_join(samplesheet_coculture, by = join_by(well)) %>%
  mutate(streptomycin = 0) %>%
  mutate(plate_name = "plate02")

# Experiment plate03 -  4 ug/ml streptomycin
plate03 <- read_logphase_xlsx(
  batch,
  "Co-culture_0_4_8ug_plates_22-syys-2025 08-17-10.xlsx",
  5,
  1
) %>%
  # rows A and H not collected here so they are NA
  drop_na() %>%
  left_join(samplesheet_coculture, by = join_by(well)) %>%
  mutate(streptomycin = 4) %>%
  mutate(plate_name = "plate03")

# Experiment plate04 - 8 ug/ml streptomycin
plate04 <- read_logphase_xlsx(
  batch,
  "Co-culture_0_4_8ug_plates_22-syys-2025 08-17-10.xlsx",
  8,
  1
) %>%
  # rows A and H not collected here so they are NA
  drop_na() %>%
  left_join(samplesheet_coculture, by = join_by(well)) %>%
  mutate(streptomycin = 8) %>%
  mutate(plate_name = "plate04")


# Batch02 - 2025.09.25 ---------------------------------------------------

batch <- "20250925_batch02"

# Read sample metadata
samplesheet_coculture <- readxl::read_xlsx(here::here(
  data_raw,
  batch,
  "coculture_samplesheet.xlsx"
))

# Read crowth curves
# Experiment plate05 - 12 ug/ml streptomycin
plate05 <- read_logphase_xlsx(
  batch,
  "Co-culture_12_16_24_ug_plates_25-syys-2025 13-42-42.xlsx",
  2,
  1
) %>%
  # rows A and H not collected here so they are NA
  drop_na() %>%
  left_join(samplesheet_coculture, by = join_by(well)) %>%
  mutate(streptomycin = 12) %>%
  mutate(plate_name = "plate05")

# Experiment plate06 =  16 ug/ml streptomycin
plate06 <- read_logphase_xlsx(
  batch,
  "Co-culture_12_16_24_ug_plates_25-syys-2025 13-42-42.xlsx",
  5,
  1
) %>%
  # rows A and H not collected here so they are NA
  drop_na() %>%
  left_join(samplesheet_coculture, by = join_by(well)) %>%
  mutate(streptomycin = 16) %>%
  mutate(plate_name = "plate06")

# Experiment plate07 - 24 ug/ml streptomycin
plate07 <- read_logphase_xlsx(
  batch,
  "Co-culture_12_16_24_ug_plates_25-syys-2025 13-42-42.xlsx",
  8,
  1
) %>%
  # rows A and H not collected here so they are NA
  drop_na() %>%
  left_join(samplesheet_coculture, by = join_by(well)) %>%
  mutate(streptomycin = 24) %>%
  mutate(plate_name = "plate07")


# Batch03 - 2025.10.06 ---------------------------------------------------

batch <- "20251006_batch03"

# Read growth curves
# Experiment plate08
plate08 <- read_logphase_xlsx(
  batch,
  "coculture_plates_8_9_10_11_06-10-2025.xlsx",
  2,
  1
) %>%
  # rows A and H not collected here so they are NA
  drop_na() %>%
  left_join(
    readxl::read_xlsx(here::here(
      data_raw,
      batch,
      "coculture_samplesheet_plate08.xlsx"
    )),
    by = join_by(well)
  ) %>%
  mutate(plate_name = "plate08")

# Experiment plate09
plate09 <- read_logphase_xlsx(
  batch,
  "coculture_plates_8_9_10_11_06-10-2025.xlsx",
  5,
  1
) %>%
  # rows A and H not collected here so they are NA
  drop_na() %>%
  left_join(
    readxl::read_xlsx(here::here(
      data_raw,
      batch,
      "coculture_samplesheet_plate09.xlsx"
    )),
    by = join_by(well)
  ) %>%
  mutate(plate_name = "plate09")

# Experiment plate10
plate10 <- read_logphase_xlsx(
  batch,
  "coculture_plates_8_9_10_11_06-10-2025.xlsx",
  8,
  1
) %>%
  # rows A and H not collected here so they are NA
  drop_na() %>%
  left_join(
    readxl::read_xlsx(here::here(
      data_raw,
      batch,
      "coculture_samplesheet_plate10.xlsx"
    )),
    by = join_by(well)
  ) %>%
  mutate(plate_name = "plate10")

# Experiment plate11
plate11 <- read_logphase_xlsx(
  batch,
  "coculture_plates_8_9_10_11_06-10-2025.xlsx",
  11,
  1
) %>%
  # rows A and H not collected here so they are NA
  drop_na() %>%
  left_join(
    readxl::read_xlsx(here::here(
      data_raw,
      batch,
      "coculture_samplesheet_plate11.xlsx"
    )),
    by = join_by(well)
  ) %>%
  mutate(plate_name = "plate11")


# Batch04 - 2025.10.30 ---------------------------------------------------

batch <- "20251030_batch04"

plate12 <- read_logphase_xlsx(
  batch,
  "coculture_plate12.xlsx",
  2,
  1
) %>%
  # rows A and H not collected here so they are NA
  drop_na() %>%
  left_join(
    readxl::read_xlsx(here::here(
      data_raw,
      batch,
      "coculture_samplesheet_plate12.xlsx"
    )),
    by = join_by(well)
  ) %>%
  mutate(plate_name = "plate12")

# Combine, tidy, format --------------------------------------------------

# combine all samples, group by plate + well, calculate rolling mean. This
# reduces the jaggedness of the curves and also makes the resulting saved data
# file smaller
coculture_gcurves_sm <- bind_rows(
  plate01,
  plate02,
  plate03,
  plate04,
  plate05,
  plate06,
  plate07,
  plate08,
  plate09,
  plate10,
  plate11,
  plate12
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
  coculture_gcurves_sm,
  here::here(data, "coculture_gcurves_smooth.tsv")
)

# Plot growth curves for later inspection --------------------------------

# quick function for saving files
plotter <- function(id) {
  ggsave(
    here::here(figs, paste0(id, ".svg")),
    plotplate(
      coculture_gcurves_sm,
      dfxy = NULL,
      unsmoothed = TRUE,
      predicted = FALSE,
      plate = id,
      rows = 6,
      cols = 12,
      page = 1,
      scales = "fixed"
    ),
    device = "svg",
    width = 12,
    height = 8
  )
}

paste0("plate", str_pad(seq(1:12), 2, "left", pad = "0")) %>%
  map(plotter)
