# MCTdownUnder

Etymology: 'MCT' = modern coexistence theory, 'downUnder' referring to Australian collaborators. Also because I am stupid

[Click here to view rendered notebooks of the analysis.](https://slhogle.github.io/MCTdownUnder/)

## Structure:

```
MCTdownUnder/
├── R/               # shared helper functions only (no analysis, no paths)
├── scripts/         # numbered, executable Quarto notebooks - all analysis lives here
├── data/
│   ├── raw/         # instrument exports, one dated directory per acquisition; never modified
│   ├── interim/     # cached model fits and other expensive intermediates
│   └── processed/   # cleaned, analysis-ready TSVs
├── output/
│   ├── figures/     # publication figures
│   └── tables/      # publication tables
└── _site/           # rendered project website, deployed to GitHub Pages by Actions
```

Every processed data file and figure is produced by a notebook in `scripts/`, so the analysis can be rebuilt
from `data/raw/` by rendering the project:

```bash
quarto render
```

## Manuscript:

### Published record

TBD

### Preprint

TBD

## Availability

Data and code in this GitHub repository (<https://github.com/slhogle/MCTdownUnder>) are provided under [GNU AGPL3](https://www.gnu.org/licenses/agpl-3.0.html).
- An archived release of the code here is available from Zenodo: <https://zenodo.org/records/EVENTUAL_ZENODO_RECORD>
