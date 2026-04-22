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
data_raw <- here::here(
  "data_raw",
  "serial_transfer_coculture",
  "20260202_serial_transfer_pilot"
)
data <- here::here(
  "data",
  "serial_transfer_coculture",
  "20260202_serial_transfer_pilot"
)
figs <- here::here(
  "figs",
  "serial_transfer_coculture",
  "20260202_serial_transfer_pilot"
)

# make processed data and figs directories if they don't exist
fs::dir_create(data)
fs::dir_create(figs)


# Read data --------------------------------------------------------------

# Read sample metadata
samplesheet_coculture <- readxl::read_xlsx(here::here(
  data_raw,
  "coculture_samplesheet.xlsx"
))

# d1 = 100,000 fold dilution
# d2 = 1000 fold dilution
# d3 = 50 fold dilution

# 48 hours, 100000 fold dilution
h48_d1 <- read_logphase_xlsx(
  "",
  "Coculture_45h_04-helmi-2026 11-54-19.xlsx",
  2,
  1
) %>%
  left_join(samplesheet_coculture, by = join_by(well)) %>%
  mutate(dilution = 100000, plate_name = "h48_d1", start_hrs_from_0 = 0)

# 48 hours, 50 fold dilution
h48_d3 <- read_logphase_xlsx(
  "",
  "Coculture_45h_04-helmi-2026 11-54-19.xlsx",
  5,
  1
) %>%
  left_join(samplesheet_coculture, by = join_by(well)) %>%
  mutate(dilution = 50, plate_name = "h48_d3", start_hrs_from_0 = 0)

# 96 hours, 100000 fold dilution
h96_d1 <- read_logphase_xlsx(
  "",
  "Co-culture_96h_06-helmi-2026 11-54-27.xlsx",
  2,
  1
) %>%
  left_join(samplesheet_coculture, by = join_by(well)) %>%
  mutate(dilution = 100000, plate_name = "h96_d1", start_hrs_from_0 = 45)

# 96 hours, 1000 fold dilution
h96_d2 <- read_logphase_xlsx(
  "",
  "Co-culture_96h_06-helmi-2026 11-54-27.xlsx",
  5,
  1
) %>%
  left_join(samplesheet_coculture, by = join_by(well)) %>%
  mutate(dilution = 1000, plate_name = "h96_d2", start_hrs_from_0 = 45)

# 96 hours, 50 fold dilution
h96_d3 <- read_logphase_xlsx(
  "",
  "Co-culture_96h_06-helmi-2026 11-54-27.xlsx",
  8,
  1
) %>%
  left_join(samplesheet_coculture, by = join_by(well)) %>%
  mutate(dilution = 50, plate_name = "h96_d3", start_hrs_from_0 = 45)

# 144 hours, 100000 fold dilution
h144_d1 <- read_logphase_xlsx(
  "",
  "Co-culture_6.2_09-helmi-2026 15-48-52.xlsx",
  2,
  1
) %>%
  left_join(samplesheet_coculture, by = join_by(well)) %>%
  mutate(dilution = 100000, plate_name = "h144_d1", start_hrs_from_0 = 96)

# 144 hours, 1000 fold dilution
h144_d2 <- read_logphase_xlsx(
  "",
  "Co-culture_6.2_09-helmi-2026 15-48-52.xlsx",
  5,
  1
) %>%
  left_join(samplesheet_coculture, by = join_by(well)) %>%
  mutate(dilution = 1000, plate_name = "h144_d2", start_hrs_from_0 = 96)

# 144 hours, 50 fold dilution
h144_d3 <- read_logphase_xlsx(
  "",
  "Co-culture_6.2_09-helmi-2026 15-48-52.xlsx",
  8,
  1
) %>%
  left_join(samplesheet_coculture, by = join_by(well)) %>%
  mutate(
    dilution = 50,
    plate_name = "h144_d3",
    start_hrs_from_0 = 96
  )

# Combine, tidy, format --------------------------------------------------

# combine all samples, group by plate + well, calculate rolling mean. This
# reduces the jaggedness of the curves and also makes the resulting saved data
# file smaller
coculture_gcurves_sm <- bind_rows(
  h48_d1,
  h48_d3,
  h96_d1,
  h96_d2,
  h96_d3,
  h144_d1,
  h144_d2,
  h144_d3,
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


# Plot growth curves with transfers combined -----------------------------

p <- coculture_gcurves_sm %>%
  mutate(hours_agg = hours + start_hrs_from_0) %>%
  mutate(
    strainID = paste0(
      str_to_upper(str_sub(evolution, 1, 1)),
      str_sub(strain, 2, 5)
    ),
    growthenv = paste(media, strep_conc_ugml, sep = "_")
  ) %>%
  mutate(
    strainID = factor(strainID, levels = c("A1287", "E1287", "A1977", "E1977"))
  ) %>%
  summarize(
    y = mean(OD600_rollmean),
    sd = sd(OD600_rollmean),
    .by = c(
      hours_agg,
      strainID,
      growthenv,
      dilution
    )
  ) %>%
  ggplot(aes(
    x = hours_agg,
  )) +
  ggplot2::geom_ribbon(
    aes(ymin = y - sd, ymax = y + sd, fill = factor(dilution)),
    alpha = 0.5,
    show.legend = FALSE
  ) +
  ggplot2::geom_line(aes(y = y, color = factor(dilution)), linewidth = 0.75) +
  facet_grid(growthenv ~ strainID) +
  labs(
    y = "OD",
    x = "Total hours from start",
    color = "Serial dilution\nfactor"
  )

ggsave(
  here::here(figs, "trial_timecourse.png"),
  p,
  device = "png",
  dpi = 300,
  width = 14,
  height = 8
)

# Plot individual growth curves for records ------------------------------

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
      rows = 8,
      cols = 12,
      page = 1,
      scales = "fixed"
    ),
    device = "svg",
    width = 12,
    height = 8
  )
}

map(
  c(
    "h48_d1",
    "h48_d3",
    "h96_d1",
    "h96_d2",
    "h96_d3",
    "h144_d1",
    "h144_d2",
    "h144_d3"
  ),
  plotter
)
