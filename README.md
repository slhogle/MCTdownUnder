# MCTdownUnder

MCT = 'modern coexistence theory.' downUnder refers to Australian collaborators. Also because I am stupid

[Click here to view rendered notebooks of the analysis.](https://slhogle.github.io/MCTdownUnder/)

## Structure:
The `_data_raw` directory should never be touched or modified! It includes raw data files obtained from instruments.

The `data` directory is where processed data projects should go. Usually, in an analysis workflow you will start with raw data, 
clean/organize it, perhaps transform it in some way, then save that product in `data` for later branches of the workflow. 

The `R/Py` directories store analysis code/scripts for the project. I prefer to keep a separate directory for each analysis language I am using in the project, but you 
may, of course, combine all code, regardless of language, into a single directory structure if you prefer.

## Manuscript:

### Published record

TBD

### Preprint

TBD

## Availability

Data and code in this GitHub repository (<https://github.com/slhogle/MCTdownUnder>) are provided under [GNU AGPL3](https://www.gnu.org/licenses/agpl-3.0.html).
The rendered project site is available at <https://slhogle.github.io/MCTdownUnder/>, which has been produced using [Quarto notebooks](https://quarto.org/). 
The content on the rendered site is released under the [CC BY 4.0.](https://creativecommons.org/licenses/by/4.0/)
This repository hosts all code and data for this project, including the code necessary to fully recreate the rendered webpage.

An archived release of the code is available from Zenodo: <https://zenodo.org/records/EVENTUAL_ZENODO_RECORD>

## Reproducibility

The project uses [`renv`](https://rstudio.github.io/renv/index.html) to create a reproducible environment to execute the code in this project. [See here](https://rstudio.github.io/renv/articles/renv.html#collaboration) for a brief overview on collaboration and reproduction of the entire project. 
To get up and running from an established repository, you could do:

``` r
install.packages("renv")
renv::restore()
```

To initiate `renv` for a new project:

``` r
options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/jammy/latest"))
install.packages("renv")
# initialize
renv::init()
# install some new packages
renv::install("tidyverse")
# record those packages in the lockfile
renv::snapshot()
```
