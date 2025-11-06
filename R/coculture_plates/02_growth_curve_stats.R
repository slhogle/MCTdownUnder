library(here)
library(tidyverse)
library(growthrates)
library(DescTools)
library(ggforce)
source(here::here("R", "utils_gcurves.R"))

# Using the tool growthrates -
# https://cran.r-project.org/web/packages/growthrates/index.html to estimate
# mu_max. I have found this works a lot better the gcplyr and is more convenient
# than using another tool outside of R. Nonparametric estimate growth rates by
# spline is very fast. Fitting to a model takes more time resources. Generally
# it is best to try multiple approaches and to visualize/check the data to make
# sure it makes sense -
# https://www.frontiersin.org/journals/ecology-and-evolution/articles/10.3389/fevo.2023.1313500/full)

# Global vars
data_raw <- here::here("_data_raw", "coculture_plates")
data <- here::here("data", "coculture_plates")
figs <- here::here("figs", "coculture_plates")

# make processed data and figs directories if they don't exist
fs::dir_create(data)
fs::dir_create(figs)

# save result for later
coculture_gcurves_sm <- readr::read_tsv(
  here::here(data, "coculture_gcurves_smooth.tsv")
) %>%
  # make uniq id
  mutate(id = paste0(plate_name, "|", well))


# Spline based estimates -------------------------------------------------

# Smoothing splines are a quick method to estimate maximum growth. The method is
# called nonparametric, because the growth rate is directly estimated from the
# smoothed data without being restricted to a specific model formula.

# From growthrates documentation:
# https://cran.r-project.org/web/packages/growthrates/growthrates.pdf

# The method was inspired by an algorithm of Kahm et al.
# (2010) (https://www.jstatsoft.org/article/view/v033i07), with different
# settings and assumptions. In the moment, spline fitting is always done with
# log-transformed data, assuming exponential growth at the time point of the
# maximum of the first derivative of the spline fit. All the hard work is done
# by function smooth.spline from package stats, that is highly user
# configurable. Normally, smoothness is automatically determined via
# cross-validation. This works well in many cases, whereas manual adjustment is
# required otherwise, e.g. by setting spar to a fixed value \[0, 1\] that also
# disables cross-validation.

# Fit splines to all growth curves
set.seed(45278)
many_spline <- growthrates::all_splines(
  OD600_rollmean ~ hours | id,
  data = coculture_gcurves_sm,
  spar = 0.5
)

readr::write_rds(many_spline, here::here(data, "coculture_spline_fits"))

# Extract results

many_spline_res <- growthrates::results(many_spline)

# Extract time of maximum specific growth rate for plotting

many_spline_xy <- purrr::map(many_spline@fits, \(x) {
  data.frame(x = x@xy[1], y = x@xy[2])
}) %>%
  purrr::list_rbind(names_to = "id") %>%
  separate_wider_delim(
    id,
    delim = "|",
    names = c("plate_name", "well")
  )

# Extract slope of maximum specific growth rate for plotting

many_spline_fitted <- purrr::map(many_spline@fits, \(x) {
  data.frame(x@FUN(x@obs$time, x@par))
}) %>%
  purrr::list_rbind(names_to = "id") %>%
  dplyr::rename(hours = time, predicted = y) %>%
  dplyr::left_join(coculture_gcurves_sm, by = dplyr::join_by(id, hours)) %>%
  dplyr::group_by(id) %>%
  # this step makes sure we don't plot fits that go outside the range of the data
  dplyr::mutate(
    predicted = dplyr::if_else(
      dplyr::between(predicted, min(OD600_rollmean), max(OD600_rollmean)),
      predicted,
      NA_real_
    )
  ) %>%
  dplyr::ungroup()


# Area under growth curve (AUC) ------------------------------------------
# Calculates AUC using `DescTools` package

many_auc_res <- coculture_gcurves_sm %>%
  dplyr::summarize(
    auc = DescTools::AUC(hours, OD600_rollmean),
    max_od = max(OD600_rollmean),
    min_od = min(OD600_rollmean),
    .by = id
  ) %>%
  dplyr::left_join(
    dplyr::distinct(dplyr::select(coculture_gcurves_sm, strain:id)),
    by = join_by(id)
  ) %>%
  dplyr::select(-id)

# Write auc output
readr::write_tsv(
  many_auc_res,
  here::here(data, "coculture_gcurve_auc_results.tsv")
)

# write spline results with strain info
many_spline_res %>%
  dplyr::left_join(
    dplyr::distinct(dplyr::select(coculture_gcurves_sm, strain:id)),
    by = join_by(id)
  ) %>%
  dplyr::select(-id) %>%
  readr::write_tsv(here::here(data, "coculture_gcurve_spline_results.tsv"))

# Make plots -------------------------------------------------------------

plotter <- function(plate) {
  ggsave(
    here::here(figs, paste0(plate, "_fitted.svg")),
    plotplate(
      many_spline_fitted,
      dfxy = many_spline_xy,
      unsmoothed = FALSE,
      predicted = TRUE,
      plate = plate,
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
