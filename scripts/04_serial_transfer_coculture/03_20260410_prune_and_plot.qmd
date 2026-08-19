# Setup ------------------------------------------------------------------

library(here)
library(tidyverse)
library(fs)
library(ggforce)

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
figs <- here::here(
  "figs",
  "serial_transfer_coculture",
  "20260410_serial_transfer"
)

# make processed data and figs directories if they don't exist
fs::dir_create(data)
fs::dir_create(figs)

# Read data --------------------------------------------------------------

# 1977
coculture_gcurves_sm_1977 <- readr::read_tsv(
  here::here(data, "1977_coculture_gcurves_smooth.tsv")
)

# 1287
coculture_gcurves_sm_1287 <- readr::read_tsv(
  here::here(data, "1287_coculture_gcurves_smooth.tsv")
)

# Plot individual growth curves for records ------------------------------

# quick function for saving files
plotter <- function(id, sp, df) {
  ggsave(
    here::here(figs, paste0(sp, "_", id, ".svg")),
    plotplate(
      df,
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

# 1977
map(
  c(
    "c1",
    "c2",
    "c3",
    "c4",
    "c5"
  ),
  \(x) plotter(x, "1977", coculture_gcurves_sm_1977)
)

# 1287
map(
  c(
    "c1",
    "c2",
    "c3",
    "c4",
    "c5"
  ),
  \(x) plotter(x, "1287", coculture_gcurves_sm_1287)
)


# Blacklists -------------------------------------------------------------

# these are well pairs where there was growth where there should be none...
# identified by inspecting the plots from above

blacklist_1287 <- tribble(
  ~strainID , ~plate_name , ~joined_well_pair ,
  "H1287"   , "c2"        , "C05C06"          ,
  "H1287"   , "c2"        , "E05E06"          ,
  "H1287"   , "c2"        , "F05F06"          ,
  "H1287"   , "c2"        , "E09E10"          ,
  "H1287"   , "c3"        , "C05C06"          ,
  "H1287"   , "c3"        , "E09E10"          ,
  "H1287"   , "c4"        , "C05C06"          ,
  "H1287"   , "c4"        , "E09E10"          ,
  "H1287"   , "c5"        , "C05C06"          ,
  "H1287"   , "c5"        , "E09E10"
)

blacklist_1977 <- tribble(
  ~strainID , ~plate_name , ~joined_well_pair ,
  "H1977"   , "c2"        , "D05D06"          ,
  "H1977"   , "c2"        , "E05E06"          ,
  "H1977"   , "c3"        , "D05D06"          ,
  "H1977"   , "c3"        , "E05E06"          ,
  "H1977"   , "c3"        , "G05G06"          ,
  "H1977"   , "c3"        , "B09B10"          ,
  "H1977"   , "c3"        , "C09C10"          ,
  "H1977"   , "c4"        , "D05D06"          ,
  "H1977"   , "c4"        , "E05E06"          ,
  "H1977"   , "c4"        , "G05G06"          ,
  "H1977"   , "c4"        , "B09B10"          ,
  "H1977"   , "c4"        , "C09C10"          ,
  "H1977"   , "c5"        , "D05D06"          ,
  "H1977"   , "c5"        , "E05E06"          ,
  "H1977"   , "c5"        , "G05G06"          ,
  "H1977"   , "c5"        , "B09B10"          ,
  "H1977"   , "c5"        , "C09C10"
)

# !! now that I think about this there is no need to include monocultures
# on these plates Because we have that information from the previous non-serial
# transfer runs. This might save us some plates in the future.

# write cleaned/filtered data

# 1287
coculture_gcurves_sm_1287 <- coculture_gcurves_sm_1287 |>
  anti_join(blacklist_1287, by = join_by(joined_well_pair, plate_name)) |>
  filter(
    (culture_type == "monoculture" &
      str_detect(well, "[BCD]01|[BCD]05|[BCD]09") |
      str_detect(well, "[EFG]02|[EFG]06|[EFG]10")) |
      culture_type == "coculture"
  )

readr::write_tsv(
  coculture_gcurves_sm_1287,
  here::here(data, "1287_coculture_gcurves_smooth.tsv")
)

# 1977
coculture_gcurves_sm_1977 <- coculture_gcurves_sm_1977 |>
  anti_join(blacklist_1977, by = join_by(joined_well_pair, plate_name)) |>
  filter(
    (culture_type == "monoculture" &
      str_detect(well, "[BCD]01|[BCD]05|[BCD]09") |
      str_detect(well, "[EFG]02|[EFG]06|[EFG]10")) |
      culture_type == "coculture"
  )

readr::write_tsv(
  coculture_gcurves_sm_1977,
  here::here(data, "1977_coculture_gcurves_smooth.tsv")
)


# Plot growth curves with transfers combined -----------------------------

cyclelines <- tibble(
  plate_name = c("c1", "c2", "c3", "c4", "c5"),
  cycle_hour = c(48, 96, 144, 192, 240)
)

p <- coculture_gcurves_sm_1977 |>
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
      evolution,
      growthenv,
      dilution,
      culture_type
    )
  ) %>%
  ggplot(aes(
    x = hours_agg,
  )) +
  geom_vline(data = cyclelines, aes(xintercept = cycle_hour), lty = 3) +
  ggplot2::geom_ribbon(
    aes(ymin = y - sd, ymax = y + sd, fill = factor(evolution)),
    alpha = 0.5,
    show.legend = FALSE
  ) +
  ggplot2::geom_line(aes(y = y, color = factor(evolution)), linewidth = 0.75) +
  facet_grid(growthenv ~ culture_type) +
  labs(
    y = "OD",
    x = "Total hours from start",
    color = "Evolution",
    title = "Pseudomonas 1977"
  )

ggsave(
  here::here(figs, "1977_competition_timecourse.svg"),
  p,
  device = "svg",
  dpi = 300,
  width = 10,
  height = 6
)
