# phiX174 Recoding — Analysis Code and Data

Analysis code and supporting data for the study of codon-recoding (deoptimization)
of the bacteriophage phiX174, measuring the effects of synonymous recoding on
viral fitness in *E. coli*.

## What this repo contains

- **`code/`** — R scripts (01–08) in run order:
  - `01_pacbio_counts.R` — PacBio genotype counts (quality filter, trim, count)
  - `02_pacbio_fitness_craig.R` — Craig's PacBio fitness calculation
  - `03_compile_flask_data.R` — Flask assay data assembly (4 Excel workbooks)
  - `04_calc_divergence.R` — Codon-divergence tables per construct
  - `05_plot_data.R` — Plaque-size and figure plots
  - `06_fitness_stats.R` — Basic glm fitness statistics (exploratory)
  - `07_fitness_model.R` — **Primary mixed model** (lme4). Produces the
    headline statistics in `docs/MODEL_NUMBERS_R.txt`. R translation of the
    original Python analysis.
  - `08_relative_fitness.R` — **Relative-fitness analyses** (same-day WT
    control). Produces `docs/CORRECTED_NUMBERS_R.txt`. R translation of the
    original Python analysis.

- **`data/`** — Input data:
  - `data/counts/` — 18 PacBio count files
  - `data/fastas/all.fasta` — WT + recoded phiX174 sequences
  - `data/design/construct_table.csv` — Per-construct metadata (137 constructs)
  - `data/assay/` — 4 Excel assay workbooks
  - `data/craig/` — Craig's modeling inputs

- **`results/`** — Intermediate and compiled data:
  - `fitness_all.csv` — Per-assay rows (471 assays: 310 construct + 73 WT + folding)
  - `pacbio_may_counts.csv`, `pacbio_may_freqs.csv` — PacBio intermediates
  - `fit_summary.csv`, `fit_summary_techreps.csv` — Fitness summaries
  - Divergence and recovery CSVs

- **`docs/`** — Reference output files:
  - `MODEL_NUMBERS.txt` — Primary model output (Python reference)
  - `MODEL_NUMBERS_R.txt` — Primary model output (R reproduction)
  - `CORRECTED_NUMBERS.txt` — Relative-fitness output (Python reference)
  - `CORRECTED_NUMBERS_R.txt` — Relative-fitness output (R reproduction)

## Excluded from this repo

- **`data/mpileup/`** (~230 MB) — Raw pileup files are too large for a typical
  GitHub repository. These were used only for the initial PacBio count
  generation (step 01) and are not needed to reproduce the statistical
  analyses (steps 07–08).
- **Python scripts** — The statistical analyses were originally written in
  Python (statsmodels). They have been translated to R (scripts 07–08) for
  this repository. The original Python source is not included.
- **BLAST artifacts** (`.ndb`, `.nhr`, `.nin`, `.nsq`, `.nsd`) and other
  deprecated intermediate files.

## How to run

### Prerequisites

R ≥ 4.3 with the following packages:

```r
install.packages(c("lme4", "stringr", "ggplot2", "Biostrings",
                   "readxl", "cowplot", "reshape2", "seqinr", "sangerseqR"))
```

### Quick start (statistical analyses)

The two statistical scripts (07 and 08) only need `results/fitness_all.csv`
and `data/design/construct_table.csv`, which are both in this repo:

```bash
cd code
Rscript 07_fitness_model.R      # writes ../docs/MODEL_NUMBERS_R.txt
Rscript 08_relative_fitness.R   # writes ../docs/CORRECTED_NUMBERS_R.txt
```

Both scripts are run from the `code/` directory. They use relative paths
(`../results/`, `../data/`).

### Full pipeline

Scripts 01–06 are the data-compilation pipeline. They are provided as-written
and are **not guaranteed to be turnkey** — some have hardcoded paths that
were adjusted but may still need local tweaks, and one (`03`) has a
commented-out read that causes a downstream error. The statistical results
in `docs/` were produced by the original (Python) analysis; the R translations
(07–08) reproduce those values.

## Data notes

- `fitness_all.csv` has 471 rows. Types: random (125), shuffle (101),
  shuffle_Hi (43), shuffle_Lo (34), worst (95), wt (73). Genes: F (148),
  G (127), H (123), wt (73).
- `construct_table.csv` has 137 rows (all attempted constructs, including
  11 that were never recovered).
- Gene lengths (codons): F = 428, G = 176, H = 329.

## Citation

See `CITATION.cff`.
