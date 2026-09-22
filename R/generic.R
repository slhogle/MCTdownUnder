# Shared lookups, colors, and plot styling for the MCTdownUnder project.
# Sourced by the notebooks in scripts/.

# Data --------------------------------------------------------------------

# The two species used throughout this project
tax <- tibble::tribble(
  ~strainID    , ~genus        , ~species       , ~label                 ,
  "HAMBI_1287" , "Citrobacter" , "koseri"       , "C. koseri 1287"       ,
  "HAMBI_1977" , "Pseudomonas" , "chlororaphis" , "P. chlororaphis 1977"
)

# Maps the short species number used in samplesheets to a display label
sp_labels <- c(
  "1287" = "C. koseri 1287",
  "1977" = "P. chlororaphis 1977"
)

# Evolutionary history: ancestral (streptomycin sensitive) vs evolved
# (streptomycin resistant)
hist_labels <- c(
  "anc" = "Ancestral (Str sensitive)",
  "evo" = "Evolved (Str resistant)"
)


# Colors and themes -------------------------------------------------------

# Strain identity is the first letter of the evolutionary history (A = ancestral,
# E = evolved) followed by the species number
straincols <- c(
  "A1287" = "orange",
  "E1287" = "purple",
  "A1977" = "dodgerblue",
  "E1977" = "limegreen"
)

# Evolutionary history alone, used where the species is already faceted out.
# The samplesheets use lowercase (anc/evo) but the imported bioscreen tables use
# uppercase (ANC/EVO), so both are keyed here.
histcols <- c(
  "anc" = "orange",
  "evo" = "purple",
  "ANC" = "orange",
  "EVO" = "purple"
)

# Default project plot styling. Returned as a list (rather than chaining
# `+`) because a `theme` object and a `guides` object can't be combined with
# `+` on their own -- only when each is added directly onto a ggplot. Adding
# a list onto a ggplot applies its elements in turn, so `p + theme_mct()`
# still works as expected.
theme_mct <- function() {
  list(
    ggplot2::theme_classic() +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_line(color = "grey90"),
        panel.grid.minor = ggplot2::element_blank(),
        strip.placement = "outside",
        strip.background = ggplot2::element_blank(),
        legend.position = "bottom",
        legend.key.size = ggplot2::unit(7, "mm")
      ),
    capit()
  )
}


# Helpers -----------------------------------------------------------------

# Builds the A1287/E1977 style strain identifier from the samplesheet columns
# `evolution` (anc/evo) and `strain` (which contains the species number)
make_strainID <- function(evolution, strain) {
  paste0(
    stringr::str_to_upper(stringr::str_sub(evolution, start = 1L, end = 1L)),
    stringr::str_extract(strain, "\\d+")
  )
}


capit <- function() {
  ggplot2::guides(
    x = guide_axis(cap = TRUE),
    y = guide_axis(cap = TRUE)
  )
}
