# ============================================================================
# phiX174 synonymous recoding - analysis code
# ============================================================================
# These scripts are the "as-written" analysis code for the study. They are
# provided for transparency and reproducibility but are NOT a turnkey
# pipeline: they assume the data/ layout in this repo, and a few inputs
# (notably the per-sample PacBio SAM/depth files and the construction-failure
# sheet) are large or private and are not distributed here. See the repo
# README for the data that IS included.
#
# Numbered in intended run order:
#   01_pacbio_counts.R        - genotype read counts from PacBio alignments
#   02_pacbio_fitness_craig.R - per-variant competitive fitness + flask check (C.R.M.)
#   03_compile_flask_data.R   - parse flask-assay sheets -> results/fitness_all.csv
#   04_calc_divergence.R      - codon-divergence tables between replicates
#   05_plot_data.R            - figure / plaque-size plots
#   06_fitness_stats.R        - simple fitness GLM
#   07_fitness_model.R        - PRIMARY mixed model (headline statistics)
#   08_relative_fitness.R     - relative-fitness (same-day WT) sensitivity analyses
#
# Conventions
#   * Run each script from the code/ directory (so ../data and ../results
#     resolve). If you move things, adjust the relative paths.
#   * fitness is measured in doublings per hour.
#   * "EN"/"EV" = the two independent pooled competition experiments;
#     mix/T40/T80 = time 0 / 40 min / 80 min.
# ============================================================================
