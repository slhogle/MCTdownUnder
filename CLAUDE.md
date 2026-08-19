# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

A research compendium (not an R package) for the MCTdownUnder project: modern coexistence theory applied to
streptomycin-resistance evolution in two bacteria. There is no `DESCRIPTION`, no `tests/`, and no `renv` —
`devtools::check()`/`devtools::test()`/`covr` do not apply here. Packages are expected to be available in the
user library: tidyverse, here, fs, readxl, lubridate, slider, ggforce, growthrates, DescTools, errors, scales.

The layout deliberately mirrors its sister project
[hambiEvoEnvCoexist](https://github.com/slhogle/hambiEvoEnvCoexist); when a convention here is ambiguous, that
repo is the reference.

Naming used everywhere in code, data, and prose:

- `1287` = *Citrobacter koseri* HAMBI_1287, `1977` = *Pseudomonas chlororaphis* HAMBI_1977
- `A`/`anc` = ancestral, streptomycin sensitive; `E`/`evo` = evolved, streptomycin resistant
- `strainID` = one-letter evolution prefix + species number, e.g. `A1287`, `E1977` (build it with
  `make_strainID()` from [R/generic.R](R/generic.R))
- `competition_pair` = the two strainIDs joined, e.g. `A1287_E1977`
- Wells are zero-padded (`A01`, not `A1`) so plate-reader output joins to samplesheets

## Core architecture: notebooks compute, R/ holds shared functions

- `scripts/**/*.qmd` are **fully executable**. Each reads from `data/raw/`, writes analysis-ready tables to
  `data/processed/<notebook dir name>/`, caches expensive fits in `data/interim/<notebook dir name>/`, and
  saves publication figures to `output/figures/`. Rendering the site re-runs the analysis.
- `R/*.R` holds only shared helpers — no analysis, no side effects, no paths. Functions there are
  package-qualified (`dplyr::filter()`, not `filter()`).
- Per-plate QC grids are **inline chunk output** (`#| output: true` with a `fig-cap`), not files. Only figures
  meant for the manuscript get `ggsave`d to `output/figures/`.

So changing an analysis means editing the `.qmd` and re-rendering it. There is no separate script to run.

## Commands

```bash
# Render one notebook (re-runs its code, rewrites its outputs)
quarto render scripts/03_coculture_plates/01_format_growthcurves.qmd

# Rebuild the whole site into _site/ (committed; uploaded by .github/workflows/static.yml)
quarto render

# Live preview while editing
quarto preview

# Force a cold rebuild, ignoring the freeze cache
rm -rf _freeze _site && quarto render
```

`execute: freeze: auto`, so a notebook only re-executes when its source changes. Delete `_freeze` to force
re-execution — worth doing when you change something in `R/`, since freeze does not track those files.

## Pipeline convention

Every notebook opens the same way: `## Libraries` (with `source(here::here("R", "..."))`), then
`## Global variables` defining `data_raw` / `data_interim` / `data_processed` / `figs` followed by
`fs::dir_create()` on the outputs, then `## Functions` for notebook-local helpers.

| Track | Notebooks | Notes |
|---|---|---|
| `01_bioscreen_strains` | `01_plot_analyze_growth` | Reads an import from hambiEvoEnvCoexist; no formatting step lives here |
| `02_biolog_ecoplates` | `01_format_growthcurves` → `02_growth_curve_stats` → `03_plot_analyze_growth` | Synergy H1; `03` joins `data/raw/carbon_compound_map.tsv` |
| `03_coculture_plates` | `01_format_growthcurves` → `02_growth_curve_stats` → `03_competition_outcomes` | Logphase 600 membrane plates; `03` emits the per-pair `competition_*_final.tsv` |
| `04_serial_transfer_coculture` | `01_pilot_20260202`, then `02_20260410_format_growthcurves` → `03_20260410_prune_and_plot` | One notebook per dated experiment; the pilot is standalone |
| `05_one_carbon_gcurves` | `01_carbon_source_platemap_format` → `02_read_format_gcurves` → `03_analyze_gcurves` | Ported from hambiEvoEnvCoexist `10_one_carbon_gcurves_logphase`; plates were run jointly with that project |

Stage 01 reads raw plate-reader files + samplesheet, joins by `well`, tags `plate_name`, adds a 5-point centered
rolling mean (`slider::slide_dbl`, `.before = 2, .after = 2`) as `OD600_rollmean`, and writes `gcurves_smooth.tsv`.
Stage 02 fits `growthrates::all_splines(OD600_rollmean ~ hours | id, spar = 0.5)` keyed on `id = "plate_name|well"`,
writes `spline_fits.rds` to `data/interim/`, and computes `DescTools::AUC` into `gcurve_auc_results.tsv`.

`R/gcurves.R` holds the plate reader IO and `plotplate()`. Both readers take a **full path**
(`read_logphase_xlsx(here::here(data_raw, file), sheet, skip)`); a single xlsx export usually holds several plates
on different sheets. `R/generic.R` holds strain labels, `straincols`/`histcols`, `theme_mct()`, and
`make_strainID()`.

## Gotchas

- `data/raw/` is committed and **never modified in place** — instrument output is the record. New acquisitions get
  their own `YYYYMMDD_description/` directory.
- The three `data/raw/*_biolog_ecoplate_*` dates come from file mtimes, not from the files themselves — the Gen5
  export template's embedded creation date is a 2011 stub. Same for `20240328_bioscreen_strains`, whose date is
  inferred from the sister project.
- `data/processed/`, `data/interim/`, `output/`, and `_site/` are committed build products. `_site/` **must** stay
  committed: the Pages workflow uploads it straight from the checkout without rendering.
- GitHub Pages must be set to deploy from **GitHub Actions**, not from a branch folder. The old `docs/` directory
  is gone.
- `scripts/04_.../03_prune_and_plot` writes `*_gcurves_pruned.tsv` rather than overwriting stage 02's
  `*_gcurves_smooth.tsv`, so partial re-renders cannot leave `data/processed/` in a half-pruned state.
- A `fig-cap` containing `: ` must be quoted or Quarto's YAML parser rejects the chunk.
- Streptomycin dosing in `05_one_carbon_gcurves` is **nested inside species**: HAMBI_1977 wells at 10 μg/ml,
  HAMBI_1287 at 512 μg/ml. No dose is shared, so species cannot be compared at matched streptomycin. Four carbon
  sources (D-mannitol, L-arginine, lactic acid, sodium acetate) were only run without streptomycin.
- The `Carbon_platemap_5` block transposition fixed in `05_one_carbon_gcurves/01` was inferred from the reads, not
  read off the spreadsheet, and it lands on the +Str columns that this project analyses (upstream only used the
  streptomycin-free columns, which the correction never touched). Re-check it against the bench notebook before
  publishing.
- `_notrack/` is a deliberate scratch area for untracked data and code (it still holds `_notrack/R/utils_cfus.R`,
  never wired into the pipeline); `wip/` likewise.
- New notebooks must be added to **both** the `render:` list and a sidebar `contents:` in `_quarto.yml`.
