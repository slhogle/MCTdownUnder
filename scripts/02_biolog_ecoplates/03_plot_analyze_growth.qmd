library(here)
library(tidyverse)
library(ggforce)
library(errors)
source(here::here("R", "utils_gcurves.R"))

# Global vars
data_raw <- here::here("_data_raw", "biolog_ecoplates")
data <- here::here("data", "biolog_ecoplates")
figs <- here::here("figs", "biolog_ecoplates")

# make processed data and figs directories if they don't exist
fs::dir_create(data)
fs::dir_create(figs)

# Read data --------------------------------------------------------------

# Growth summary data (e.g., AUC and growth rates from )
many_auc_res <- readr::read_tsv(here::here(
  data,
  "gcurves_auc_results.tsv"
))

many_spline_res <- readr::read_tsv(here::here(
  data,
  "gcurves_spline_results.tsv"
))

# Full growth curves
gcurves <- readr::read_tsv(here::here(data, "gcurves_smooth.tsv")) %>%
  arrange(strainID, well, plate_name) %>%
  group_by(strainID, evolution, plate_name, `carbon source`, well) %>%
  mutate(id = dplyr::cur_group_id()) %>%
  group_by(strainID, evolution, `carbon source`) %>%
  mutate(
    replicate = case_when(
      id == min(id) ~ 1,
      id == max(id) ~ n_distinct(id),
      TRUE ~ 2
    )
  ) %>%
  mutate(replicate = LETTERS[replicate]) %>%
  dplyr::select(-id) %>%
  ungroup() %>%
  rename(carbon_source = `carbon source`) %>%
  mutate(
    hist = if_else(str_detect(evolution, "anc"), "Str Sens.", "Str Res.")
  )

# read which wells correspond to which compounds
cc <- readr::read_tsv(here::here(data_raw, "carbon_compound_map.tsv")) %>%
  dplyr::select(compound, class)

# Plot ecoplate growthcurves ---------------------------------------------

# HAMBI_1287
p1287ecoplates_gcurves <- gcurves %>%
  filter(strainID == "HAMBI_1287") %>%
  ggplot() +
  ggplot2::geom_line(aes(
    x = hours,
    y = OD600_rollmean,
    color = hist,
    linetype = replicate,
    group = interaction(well, plate_name)
  )) +
  ggplot2::labs(x = "Hours", y = "OD600") +
  ggplot2::scale_x_continuous(
    breaks = seq(0, 48, 12),
    labels = seq(0, 48, 12)
  ) +
  labs(
    y = "Metabolic performance (absorbance units, λ = 590 nm)",
    x = "Hours",
    color = "",
    lintype = "Replicate",
    title = "Citrobacter koseri HAMBI_1287"
  ) +
  facet_wrap(
    ~carbon_source,
    labeller = labeller(carbon_source = label_wrap_gen(20))
  ) +
  theme(legend.position = "bottom")

ggsave(
  here::here(figs, "H1287_ecoplates_gcurves.svg"),
  plot = p1287ecoplates_gcurves,
  device = "svg",
  width = 7,
  height = 8
)

# HAMBI_1977
p1977ecoplates_gcurves <- gcurves %>%
  filter(strainID == "HAMBI_1977") %>%
  ggplot() +
  ggplot2::geom_line(aes(
    x = hours,
    y = OD600_rollmean,
    color = hist,
    linetype = replicate,
    group = interaction(well, plate_name)
  )) +
  ggplot2::labs(x = "Hours", y = "OD600") +
  ggplot2::scale_x_continuous(
    breaks = seq(0, 48, 12),
    labels = seq(0, 48, 12)
  ) +
  facet_wrap(~carbon_source) +
  labs(
    y = "Metabolic performance (absorbance units, λ = 590 nm)",
    x = "Hours",
    color = "",
    lintype = "Replicate",
    title = "Pseudomonas chlororaphis HAMBI_1977"
  ) +
  facet_wrap(
    ~carbon_source,
    labeller = labeller(carbon_source = label_wrap_gen(20))
  ) +
  theme(legend.position = "bottom")

ggsave(
  here::here(figs, "H1977_ecoplates_gcurves.svg"),
  plot = p1977ecoplates_gcurves,
  device = "svg",
  width = 7,
  height = 8
)

# Plot ecoplate growth summaries -----------------------------------------

# 1. for each replicate take the mean over all measurements within 0.05 units of
#    the max
# 2. subtract water from all values
# 3. calculate mean and bootstrapped 95% CI across biological replicates

# order carbon levels by compound class
carbon_levels <- unique(paste(cc$class, cc$compound, sep = ":")) %>%
  sort() %>%
  str_extract(".*\\:(.*$)", group = 1)

gcurves_metabolic_pref <- gcurves %>%
  mutate(
    strainID2 = paste0(
      toupper(stringr::str_sub(evolution, start = 1, end = 1)),
      strain
    )
  ) %>%
  filter(hours > 24) %>%
  group_by(strainID2, carbon_source, replicate) %>%
  filter(OD600 >= max(OD600) - 0.05) %>%
  summarize(OD600 = mean(OD600)) %>%
  ungroup() %>%
  pivot_wider(names_from = carbon_source, values_from = "OD600") %>%
  pivot_longer(
    cols = c(-strainID2, -replicate, -water),
    names_to = "carbon_source",
    values_to = "OD600"
  ) %>%
  mutate(OD600_norm = if_else((OD600 - water) < 0, 0, OD600 - water)) %>%
  summarize(
    ggplot2::mean_cl_boot(OD600),
    .by = c(strainID2, carbon_source)
  )

gcurves_metabolic_pref_fmt <- gcurves_metabolic_pref %>%
  mutate(
    evolution = if_else(str_starts(strainID2, "A"), "Sens", "Res"),
    strainID2 = str_replace(strainID2, "A|E", "H")
  ) %>%
  pivot_wider(
    id_cols = c(strainID2, carbon_source),
    names_from = evolution,
    values_from = c(y, ymin, ymax)
  ) %>%
  mutate(
    y_Sens = errors::set_errors(y_Sens, ymax_Sens - y_Sens),
    y_Res = errors::set_errors(y_Res, ymax_Res - y_Res)
  ) %>%
  mutate(delta_sens_res = y_Sens - y_Res) %>%
  mutate(delta_sens_res_cl95 = errors::errors(delta_sens_res)) %>%
  mutate(delta_sens_res = errors::drop_errors(delta_sens_res)) %>%
  left_join(cc, by = join_by(carbon_source == compound)) %>%
  mutate(
    tot_met_output = errors::drop_errors(y_Sens + y_Res)
  )

carbon_levels <- gcurves_metabolic_pref_fmt %>%
  filter(strainID2 == "H1287") %>%
  arrange(delta_sens_res) %>%
  arrange(class) %>%
  distinct(carbon_source) %>%
  pull()

# plot
pcarbonindv <- gcurves_metabolic_pref_fmt %>%
  mutate(carbon_source = factor(carbon_source, levels = carbon_levels)) %>%
  ggplot(
    aes(
      x = carbon_source,
      y = delta_sens_res,
      color = class
    ),
  ) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point(aes(size = tot_met_output)) +
  geom_linerange(aes(
    ymin = delta_sens_res - delta_sens_res_cl95,
    ymax = delta_sens_res + delta_sens_res_cl95
  )) +
  labs(
    x = "",
    y = "Δ metabolic output (Sens - Res)",
    color = "Carbon class",
    size = "Total metabolic\noutput (Sens + Res)"
  ) +
  geom_line(aes(group = 1), alpha = 0.5) +
  coord_flip() +
  facet_grid(~strainID2)

ggsave(
  here::here(figs, "summary_metabolic_output.svg"),
  plot = pcarbonindv,
  device = "svg",
  width = 7,
  height = 7
)

# Old code - not needed --------------------------------------------------

# Code for plotting heatmaps of growth parameters
# Plot maximum observed absorbance
# Some clustering to plot by more similar responses

# df2clust <- many_auc_res %>%
#   mutate(strainID2 = paste0(toupper(evolution), "_", strain)) %>%
#   summarize(y = min(max_od), .by = c(strainID2, `carbon source`)) %>%
#   pivot_wider(names_from = "carbon source", values_from = "y") %>%
#   pivot_longer(c(-strainID2, -water)) %>%
#   mutate(value_norm = value - water) %>%
#   pivot_wider(
#     id_cols = c(-water, -value),
#     names_from = "strainID2",
#     values_from = "value_norm"
#   ) %>%
#   as.data.frame() %>%
#   column_to_rownames(var = "name")

# # scale the dataframe for clustering
# df2clust_scaled <- scale(df2clust)

# # perform the hierarchcical clustering using euclidean distance and Ward's D
# hc <- hclust(dist(df2clust_scaled, method = "euclidean"), method = "ward.D2")

# # plot
# df2clust %>%
#   rownames_to_column(var = "carbon") %>%
#   pivot_longer(cols = -carbon) %>%
#   mutate(
#     carbon = factor(carbon, levels = hc$labels[hc$order]),
#     name = factor(
#       name,
#       levels = c("ANC_1287", "EVO_1287", "ANC_1977", "EVO_1977")
#     )
#   ) %>%
#   ggplot(aes(x = name, y = carbon)) +
#   geom_tile(aes(fill = value)) +
#   scale_x_discrete(guide = guide_axis(angle = 90)) +
#   scale_fill_viridis_c() +
#   scale_color_manual(values = c("white", "black"), guide = "none") +
#   labs(y = "Carbon substrate", fill = "Maximum\nabsorbance", x = "") +
#   coord_fixed() +
#   ggplot2::theme(
#     panel.grid = element_blank(),
#     panel.background = element_blank(),
#     strip.background = element_blank(),
#     panel.border = element_blank()
#   )

# many_auc_res %>%
#   mutate(strainID2 = paste0(toupper(evolution), "_", strain)) %>%
#   summarize(y = min(max_od), .by = c(strainID2, `carbon source`)) %>%
#   #summarize(ggplot2::mean_cl_boot(max_od), .by=c(strainID2, `carbon source`)) %>%
#   mutate(
#     strainID2 = factor(
#       strainID2,
#       levels = c("ANC_1287", "EVO_1287", "ANC_1977", "EVO_1977")
#     )
#   ) %>%
#   ggplot(aes(x = strainID2, y = `carbon source`)) +
#   geom_tile(aes(fill = y)) +
#   scale_x_discrete(guide = guide_axis(angle = 90)) +
#   scale_fill_viridis_c() +
#   scale_color_manual(values = c("white", "black"), guide = "none") +
#   labs(y = "Carbon substrate", fill = "Maximum\nabsorbance", x = "") +
#   coord_fixed() +
#   ggplot2::theme(
#     panel.grid = element_blank(),
#     panel.background = element_blank(),
#     strip.background = element_blank(),
#     panel.border = element_blank()
#   )

# # Plot AUC

# many_auc_res %>%
#   mutate(strainID2 = paste0(toupper(evolution), "_", strain)) %>%
#   summarize(
#     ggplot2::mean_cl_boot(auc),
#     .by = c(
#       strainID2,
#       `carbon source`
#     )
#   ) %>%
#   mutate(
#     strainID2 = factor(
#       strainID2,
#       levels = c("ANC_1287", "EVO_1287", "ANC_1977", "EVO_1977")
#     )
#   ) %>%
#   ggplot(aes(
#     x = strainID2,
#     y = `carbon source`
#   )) +
#   geom_tile(aes(fill = y)) +
#   scale_x_discrete(guide = guide_axis(angle = 90)) +
#   scale_fill_viridis_c() +
#   scale_color_manual(values = c("white", "black"), guide = "none") +
#   labs(y = "Carbon substrate", fill = "Area under the\ngrowth curve", x = "") +
#   coord_fixed() +
#   ggplot2::theme(
#     panel.grid = element_blank(),
#     panel.background = element_blank(),
#     strip.background = element_blank(),
#     panel.border = element_blank()
#   )

# # Plot growth rate
# many_spline_res %>%
#   mutate(strainID2 = paste0(toupper(evolution), "_", strain)) %>%
#   summarize(
#     ggplot2::mean_cl_boot(mumax),
#     .by = c(strainID2, `carbon source`)
#   ) %>%
#   mutate(
#     strainID2 = factor(
#       strainID2,
#       levels = c("ANC_1287", "EVO_1287", "ANC_1977", "EVO_1977")
#     )
#   ) %>%
#   ggplot(aes(x = strainID2, y = `carbon source`)) +
#   geom_tile(aes(fill = y)) +
#   scale_x_discrete(guide = guide_axis(angle = 90)) +
#   scale_fill_viridis_c() +
#   scale_color_manual(values = c("white", "black"), guide = "none") +
#   labs(y = "Carbon substrate", fill = "Growth rate\n(hr-1)", x = "") +
#   coord_fixed() +
#   ggplot2::theme(
#     panel.grid = element_blank(),
#     panel.background = element_blank(),
#     strip.background = element_blank(),
#     panel.border = element_blank()
#   )
