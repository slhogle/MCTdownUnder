# Libraries
library(tidyverse)
library(here)
library(fs)
library(scales)

# Global vars
data <- here::here("data", "bioscreen_strains")
figs <- here::here("figs", "bioscreen_strains")

# make processed data and figs directories if they don't exist
fs::dir_create(data)
fs::dir_create(figs)

# Read data
many_auc_res <- readr::read_tsv(here::here(data, "gcurve_auc_results.tsv"))
many_spline_res <- readr::read_tsv(here::here(
  data,
  "gcurve_spline_results.tsv"
))
gcurves <- readr::read_tsv(here::here(data, "gcurves_formatted_thinned.tsv"))


# Growth curves ----------------------------------------------------------

pgcurves <- gcurves %>%
  mutate(
    hist = str_split_i(sp_hist, "_", 1),
    sp = str_split_i(sp_hist, "_", 2),
    strep_conc = factor(strep_conc, ordered = TRUE)
  ) %>%
  mutate(
    sp = if_else(sp == "1287", "C. koseri 1287", "P. chlororaphis 1977")
  ) %>%
  ggplot(aes(x = hours, y = OD600_rollmean, color = hist)) +
  geom_line(aes(group = interaction(sp_hist, bioscreen_well))) +
  labs(
    y = "Optical density (OD600)",
    x = "Time (hours)",
    color = "Evolutionary history"
  ) +
  facet_grid(strep_conc ~ sp) +
  theme(legend.position = "bottom")

ggsave(
  here::here(figs, "growthcurves_monoculture.svg"),
  plot = pgcurves,
  device = "svg",
  width = 7,
  height = 12
)

# Growth rates vs streptomycin level -------------------------------------

gr_strep <- many_spline_res %>%
  mutate(
    hist = str_split_i(sp_hist, "_", 1),
    sp = str_split_i(sp_hist, "_", 2)
  ) %>%
  summarize(ggplot2::mean_cl_boot(mumax), .by = c(sp, hist, strep_conc)) %>%
  mutate(strep_conc = factor(strep_conc, ordered = TRUE)) %>%
  ggplot(aes(x = strep_conc, y = y)) +
  geom_linerange(aes(ymin = ymin, ymax = ymax, color = hist)) +
  geom_line(aes(color = hist, group = hist), lty = 1) +
  geom_point(aes(color = hist)) +
  labs(
    y = "Maximum per capita growth rate μ (hr-1)",
    x = "Streptomycin conc. (μg/ml)",
    color = "Evolutionary\nhistory"
  ) +
  facet_grid(~sp) +
  scale_x_discrete(guide = guide_axis(angle = 90))

ggsave(
  here::here(figs, "u_streptomycin_monoculture.svg"),
  plot = gr_strep,
  device = "svg",
  width = 7,
  height = 4
)

# Maximum carrying capcaity (K, maximum observed OD) ---------------------

k_strep <- many_auc_res %>%
  mutate(
    hist = str_split_i(sp_hist, "_", 1),
    sp = str_split_i(sp_hist, "_", 2)
  ) %>%
  summarize(ggplot2::mean_cl_boot(max_od), .by = c(sp, hist, strep_conc)) %>%
  mutate(strep_conc = factor(strep_conc, ordered = TRUE)) %>%
  ggplot(aes(x = strep_conc, y = y)) +
  geom_linerange(aes(ymin = ymin, ymax = ymax, color = hist)) +
  geom_line(aes(color = hist, group = hist), lty = 1) +
  geom_point(aes(color = hist)) +
  labs(
    y = "Maximum carrying capacity (K, OD units)",
    x = "Streptomycin conc. (μg/ml)",
    color = "Evolutionary\nhistory"
  ) +
  facet_grid(~sp) +
  scale_x_discrete(guide = guide_axis(angle = 90))

ggsave(
  here::here(figs, "k_streptomycin_monoculture.svg"),
  plot = k_strep,
  device = "svg",
  width = 7,
  height = 4
)

# Area under the curve (AUC) vs streptomycin level ---------------------

# not run
# auc_strep <- many_auc_res %>%
#   mutate(
#     hist = str_split_i(sp_hist, "_", 1),
#     sp = str_split_i(sp_hist, "_", 2)
#   ) %>%
#   summarize(ggplot2::mean_cl_boot(auc), .by = c(sp, hist, strep_conc)) %>%
#   mutate(strep_conc = factor(strep_conc, ordered = TRUE)) %>%
#   ggplot(aes(x = strep_conc, y = y)) +
#   geom_linerange(aes(ymin = ymin, ymax = ymax, color = hist)) +
#   geom_line(aes(color = hist, group = hist), lty = 1) +
#   geom_point(aes(color = hist)) +
#   labs(
#     y = "Total area under the growth curve",
#     x = "Streptomycin conc. (μg/ml)",
#     color = "Evolutionary\nhistory"
#   ) +
#   facet_grid(~sp) +
#   scale_x_discrete(guide = guide_axis(angle = 90))
