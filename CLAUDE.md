# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

A research compendium (not an R package) for the MCTdownUnder project: modern coexistence theory applied to
streptomycin-resistance evolution in two bacteria. There is no `DESCRIPTION`, no `tests/`, and no `renv` —
`devtools::check()`/`devtools::test()`/`covr` do not apply here. Packages are expected to be available in the
user library: tidyverse, here, fs, readxl, lubridate, slider, ggforce, growthrates, DescTools, scales.

Naming used everywhere in code, data, and prose:

- `1287` = *Citrobacter koseri* HAMBI_1287, `1977` = *Pseudomonas chlororaphis* HAMBI_1977
- `A`/`anc` = ancestral, streptomycin sensitive; `E`/`evo` = evolved, streptomycin resistant
- `strainID` = one-letter evolution prefix + species number, e.g. `A1287`, `E1977`
- `competition_pair` = the two strainIDs joined, e.g. `A1287_E1977`
- Wells are zero-padded (`A01`, not `A1`) so plate-reader output joins to samplesheets

## Core architecture: scripts compute, notebooks only narrate

This is the single most important thing to understand before editing anything.

- `R/**/*.R` does **all** computation: reads `data_raw/`, writes tidy TSVs to `data/`, writes SVGs to `figs/`.
- `notebooks/**/*.qmd` contain **zero executable code chunks**. They are prose plus `![](/figs/...)` image
  references (root-absolute paths) and Quarto cross-references. Rendering the site does not recompute anything.

So a change to an analysis is a three-step loop: edit the R script → run it (regenerates TSV/SVG) → `quarto render`
to rebuild the site. Editing a `.qmd` alone can never change a figure.

## Commands

```bash
# Run a pipeline step (always from the repo root; here::here() anchors all paths)
Rscript R/coculture_plates/01_format_growthcurves.R

# Rebuild the whole website into docs/ (committed; served by GitHub Pages)
quarto render

# Live preview while editing prose
quarto preview

# Render one notebook only
quarto render notebooks/coculture_plates/01_inspect_growthcurves.qmd
```

## Pipeline convention

Each experiment track is a sibling directory under `R/`, `data/`, and `figs/`, with numbered scripts run in order.
Every script opens the same way: `library()` calls, `source(here::here("R", "utils_gcurves.R"))`, then path globals
(`data_raw`, `data`, `figs`) followed by `fs::dir_create()` on the outputs.

| Track | Scripts | Notes |
|---|---|---|
| `biolog_ecoplates` | `01_format_growthcurves` → `02_growth_curve_stats` → `03_plot` | Synergy H1 reader; xlsx parsed inline (not via the shared helper). `03_plot` joins `data_raw/biolog_ecoplates/carbon_compound_map.tsv` |
| `coculture_plates` | `01_format_growthcurves` → `02_growth_curve_stats` → `03_plot` | Logphase 600 membrane co-culture plates; `03_plot` emits the per-pair `competition_*_final.tsv` |
| `serial_transfer_coculture/<YYYYMMDD_name>/` | `01_format_growthcurves` → `02_prune_and_plot` (pilot: `02_growth_curve_stats`) | One dated subdirectory per experiment; no `03_plot` step |
| `bioscreen_strains` | `01_plot_analyze` only | Reads TSVs already present in `data/bioscreen_strains/`; the formatting step lives outside this repo |

Stage 01 reads raw plate-reader files + samplesheet, joins by `well`, tags `plate_name`, adds a 5-point centered
rolling mean (`slider::slide_dbl`, `.before = 2, .after = 2`) as `OD600_rollmean`, and writes `*_gcurves_smooth.tsv`.
Stage 02 fits `growthrates::all_splines(OD600_rollmean ~ hours | id, spar = 0.5)` keyed on `id = "plate_name|well"`,
writes the fit object plus `*_spline_results.tsv`, and computes `DescTools::AUC` into `*_auc_results.tsv`.
Both stages save one paginated 6×12 SVG per plate via the shared `plotplate()` helper (`_fitted.svg` suffix for stage 02).

`R/utils_gcurves.R` holds the only shared code: `read_logphase_xlsx(subdir, filename, sheet, skip)` (a single xlsx
often holds several plates on different sheets) and `plotplate()`. Note `read_logphase_xlsx` resolves paths against
the caller's `data_raw` global — it depends on that variable existing in the calling environment.

## Gotchas

- **`_data_raw` is stale.** `R/biolog_ecoplates/*` and `R/coculture_plates/*` (plus the serial-transfer pilot's
  `02_`) still set `data_raw <- here::here("_data_raw", ...)`, but the directory is `data_raw`. Re-running those
  01/03 scripts fails until the path is fixed; the newer `serial_transfer_coculture/20260410_*` scripts are correct.
- `data_raw/` is committed and **never modified in place** — instrument output is the record.
- `data/`, `figs/`, and `docs/` are committed build products. `docs/` is regenerated wholesale by `quarto render`,
  so review its diff loosely but do commit it (GitHub Pages serves from it).
- `.gitignore` excludes `*.rds`/`*.parquet`, so `growthrates` fit objects are written extensionless
  (e.g. `data/coculture_plates/coculture_spline_fits`) to stay tracked.
- `_notrack/` is a deliberate scratch area for untracked data and code; `wip/` likewise.
- New notebooks must be added to the sidebar `contents:` in `_quarto.yml` or they will not appear on the site.
