# hogleLabProjTemplate

A simple template for a project using R/Python and [Quarto](https://quarto.org/) notebooks. 

The project is rendered with Quarto and can then be served using GitHub Pages [following these instructions](https://quarto.org/docs/publishing/github-pages.html). 

You can view the HTML documents on the web and share them with collaborators. Once rendered and configured with GitHub pages, 
you would be able to view the website at this URL <https://GHUSERNAME.github.io/REPOSITORYNAME/>

[Click here to view rendered notebooks of the analysis.](https://slhogle.github.io/hogleLabProjTemplate/)

## Structure:
The `_data_raw` directory should never be touched or modified! It includes raw data files obtained from instruments.

The `data` directory is where processed data projects should go. Usually, in an analysis workflow you will start with raw data, 
clean/organize it, perhaps transform it in some way, then save that product in `data` for later branches of the workflow. 

The `R/Py` directories store analysis code/scripts for the project. I prefer to keep a separate directory for each analysis language I am using in the project, but you 
may, of course, combine all code, regardless of language, into a single directory structure if you prefer.

Note: You can create whatever kind of sub-directory structure you prefer within `_data_raw`, `data`, `R`, and `Py`. You can also create other new directories. 
Some steps of the process will create new directories for you (e.g., running `renv::init`, or `quarto render`). This is expected. Renv will create its own `.gitignore` to prevent 
committing huge R environments.

## Manuscript:

◇ Corresponding author

### Published record

**Title XYZ**\
FIRST AUTHOR<sup>◇</sup>, ..., LAST AUTHOR<sup>◇</sup>. *XYZ* (2025) [doi:]()

### Preprint

**Title XYZ**\
FIRST AUTHOR<sup>◇</sup>, ..., LAST AUTHOR<sup>◇</sup>. *BioRxiv* (2025) [doi:]()

## Availability

Data and code in this GitHub repository (<https://github.com/GHUSERNAME/REPOSITORYNAME>) are provided under [GNU AGPL3](https://www.gnu.org/licenses/agpl-3.0.html).
The rendered project site is available at <https://GHUSERNAME.github.io/REPOSITORYNAME/>, which has been produced using [Quarto notebooks](https://quarto.org/). 
The content on the rendered site is released under the [CC BY 4.0.](https://creativecommons.org/licenses/by/4.0/)
This repository hosts all code and data for this project, including the code necessary to fully recreate the rendered webpage.

An archived release of the code is available from Zenodo: <https://zenodo.org/records/EVENTUAL_ZENODO_RECORD>

Raw sequencing data used in the project is available from NCBI Bioproject [PRJNA00000000](https://www.ncbi.nlm.nih

## Reproducibility

The project uses [`renv`](https://rstudio.github.io/renv/index.html) to create a reproducible environment to execute the code in this project. [See here](https://rstudio.github.io/renv/articles/renv.html#collaboration) for a brief overview on collaboration and reproduction of the entire project. 
To get up and running from an established repository, you could do:

``` r
install.packages("renv")
renv::restore()
```

To initiate `renv` for a new project:

``` r
install.packages("renv")
# initialize
renv::init()
# install some new packages
renv::install("tidyverse")
# record those packages in the lockfile
renv::snapshot()
```
