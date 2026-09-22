# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

A research compendium (not an R package) for the MCTdownUnder project: modern coexistence theory applied to
streptomycin-resistance evolution in two bacteria. There is no `DESCRIPTION`, no `tests/`, and no `renv` —
`devtools::check()`/`devtools::test()`/`covr` do not apply here. Packages are expected to be available in the
user library: tidyverse, here, fs, readxl, lubridate, slider, ggforce, growthrates, DescTools, errors, scales,
glue, marquee, bayesboot.

The layout deliberately mirrors its sister project
[hambiEvoEnvCoexist](https://github.com/slhogle/hambiEvoEnvCoexist); when a convention here is ambiguous, that
repo is the reference.

Naming used everywhere in code, data, and prose:

- `1287` = *Citrobacter koseri* HAMBI_1287, `1977` = *Pseudomonas chlororaphis* HAMBI_1977
- `A`/`anc` = ancestral, streptomycin sensitive; `E`/`evo` = evolved, streptomycin resistant
- `strainID` = one-letter evolution prefix + species number, e.g. `A1287`, `E1977` (build it with
  `make_strainID()` from [R/generic.R](R/generic.R)). The `05` track is the exception; see Gotchas.
- `competition_pair` = the two strainIDs joined, e.g. `A1287_E1977`
- Wells are zero-padded (`A01`, not `A1`) so plate-reader output joins to samplesheets

## Core architecture: notebooks compute, R/ holds shared functions

- `scripts/**/*.qmd` are **fully executable**. Each reads from `data/raw/` (tracks `01` and `05` instead start from
  an imported table), writes analysis-ready tables to `data/processed/<notebook dir name>/`, caches expensive fits
  in `data/interim/<notebook dir name>/`, and saves publication figures to `output/figures/`. Rendering the site
  re-runs the analysis.
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
| `02_biolog_ecoplates` | `01_format_growthcurves` → `02_growth_curve_stats` → `03_plot_analyze_growth` | Synergy H1; `03` joins `data/raw/ecoplate_carbon_compound_map.tsv` |
| `03_coculture_plates` | `01_format_growthcurves` → `02_growth_curve_stats` → `03_competition_outcomes` | Logphase 600 membrane plates; `03` emits the per-pair `competition_*_final.tsv` |
| `04_serial_transfer_coculture` | `01_pilot_20260202`, then `02_20260410_format_growthcurves` → `03_20260410_prune_and_plot` | One notebook per dated experiment; the pilot is standalone |
| `05_one_carbon_gcurves` | `01_analyze_gcurves`, `02_analyze_gcurves_summary` (independent of each other) | No formatting step or raw files live here; both start from the imported `data/processed/05_one_carbon_gcurves/growth_curves_sm_bl_filt.tsv` (plates run jointly with hambiEvoEnvCoexist). `01` plots curves per carbon class; `02` writes `growth_summaries_carbon.tsv` and the K/AUC/μ figures |
| `06_cfu_competition` | `01_cfu_competition_outcomes` | Standalone; colony counts, not plate reader. Uses `R/cfus.R` maximum likelihood CFU estimators |

Stage 01 reads raw plate-reader files + samplesheet, joins by `well`, tags `plate_name`, adds a 5-point centered
rolling mean (`slider::slide_dbl`, `.before = 2, .after = 2`) as `OD600_rollmean`, and writes `gcurves_smooth.tsv`.
Stage 02 fits `growthrates::all_splines(OD600_rollmean ~ hours | id, spar = 0.5)` keyed on `id = "plate_name|well"`,
writes `spline_fits.rds` to `data/interim/`, and computes `DescTools::AUC` into `gcurve_auc_results.tsv`.

`R/cfus.R` holds the colony-forming-unit estimators of [Martini et al. 2024](https://doi.org/10.1128/spectrum.03946-23).
`find_estimators()` takes counts plus the dilution *fraction* (`1e-6`, not `-6`) and returns CFU in the volume
plated, so callers divide by the plating volume to get CFU/ml. Two traps:

- `find_MLE()` brackets its root on `c(1, 1e9)` and errors with "f() values at end points not of opposite sign"
  above that. So do not fold the plating volume into `dil` (use `V = 1` and divide afterwards), and guard against
  a series implying more than 1e9 CFU per plated volume — 100 colonies at 10^-7^ already breaks it.
- Pass **zero counts through** to the estimators. Both pool `sum(counts)/sum(dil)`, so an empty plate is what
  records the volume examined for nothing; filtering `count > 0` row-wise inflates the estimate (111 fold for one
  vial in `06_cfu_competition`). Withhold only series that are zero at *every* dilution — those have no
  root — and report them at a detection limit instead. Older sister-project scripts filter row-wise; they are wrong.

`R/gcurves.R` holds the plate reader IO and `plotplate()`. Both readers take a **full path**
(`read_logphase_xlsx(here::here(data_raw, file), sheet, skip)`); a single xlsx export usually holds several plates
on different sheets. It also holds the growth summaries used by `05`: `trap_auc()`, `dirty_mean_ci()` (mean ±
1.96 SE), and `tidy_bayes_boot()` (Bayesian bootstrap via `bayesboot`). `tidy_bayes_boot()` names its point
estimate `md` (a median), not `mn` like the others — except on empty input, where it returns `mn`. `R/generic.R`
holds strain labels, `straincols`/`histcols`, `theme_mct()`, and `make_strainID()`.

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
  HAMBI_1287 at 512 μg/ml. No dose is shared, so species cannot be compared at matched streptomycin. Every carbon
  source except cellobiose has a +Str arm. The +Str arms for D-mannitol, L-asparagine, lactic acid and sodium
  acetate came from a later run (28 Aug 2026), so older renders that call those four "never run with
  streptomycin" are stale.
- In the `05` input table, `strainID` is the **species** (`H1287`, `H1977`) and `evo_hist` is `ANC`/`EVO`, not the
  `A1287` form used elsewhere; `growth_summaries_carbon.tsv` renames them `species` and `genotype` (`A`/`E`). The
  imported table also holds `H0403` and `H1896` from the jointly run plates, so always filter on `strainID`.
- Wells labelled `L-arginine` in `05` outputs from before commit `1e5ebc0` are actually **L-asparagine**.
- `05` carbon sources were dosed at 10 mM, and citric, malonic, glycolic and salicylic acid went in as **free
  acids**. 10 mM citric acid lowers M9 by ~0.75 pH units, so zero growth on citrate (and on malonate for
  *C. koseri*, which is classically malonate-positive) may be a pH artifact. A rerun with neutralised stocks was
  pending as of 2026-09-17.
- `05` **lactic acid + streptomycin** shows no growth in any well, including the resistant EVO strains, so it is an
  artifact, not a drug effect. Exclude those wells until the rerun.
- `05_.../02_analyze_gcurves_summary` does **not** refit on render: its spline-fit and bootstrap chunks are
  `eval: false` and read the cached `data/interim/05_one_carbon_gcurves/one_carbon_spline_fits.rds` and
  `data/processed/05_one_carbon_gcurves/auc_k_mu_growth_summaries.rds`. Flip those chunks to `eval: true` to
  regenerate. `data/interim/05_one_carbon_gcurves/platemap_master_tidy.tsv` is left over from the removed
  formatting notebook; nothing reads it.
- Prefer **AUC** over K when comparing growth in `05`: streptomycin mostly delays the ancestors by 10–30 h and they
  then recover, which K scores as no effect.
- The `Carbon_platemap_5` block transposition was corrected by the old `05` formatting notebook, which is no longer
  in this repo. The fix was inferred from the reads, not read off the spreadsheet, and it lands on the +Str columns
  that this project analyses (upstream only used the streptomycin-free columns). Re-check it against the bench
  notebook before publishing.
- In `20260604_traditional_competition_cfu` the streptomycin in the **vials** (10 μg/ml for 1977, 512 μg/ml for
  1287) is **subinhibitory** — the sensitive ancestor still grows at it. The treatment asks whether ANC still
  outcompetes EVO under drug and whether they coexist; it is not a kill step, so the resistant fraction is not
  expected near 100%. Only the 5000 μg/ml in the **plating** agar is fully inhibitory, which is what makes a
  colony there a resistant cell.
- Ten selective platings in that dataset (nine of them 1287 at cycle 5) grew nothing at any dilution. A non-detect
  is reported at the detection limit, 3.3e6 CFU/ml, and flagged `detected == FALSE`; it cannot distinguish
  extinction from persistence below ~0.04–0.24% of the total. Treat those fractions as upper bounds.
- `transfer_cycle` in that file is the number of completed 48 h growth cycles: 1 = day 2, 5 = day 10. An earlier
  version of the file mislabelled cycle 1 as cycle 2.
- `_notrack/` is a deliberate scratch area for untracked data and code (it still holds `_notrack/R/utils_cfus.R`,
  never wired into the pipeline); `wip/` likewise.
- New notebooks must be added to **both** the `render:` list and a sidebar `contents:` in `_quarto.yml`.
