# ============================================================================
# 06_fitness_stats.R
# Basic fitness statistics: a quick glm of fitness ~ gene + type + treat
# for an initial look at the compiled data (superseded by 07/08).
#
# Data dependencies (relative to repo root):
#   * results/fitness_all.csv
#   * results/fit_summary_techreps.csv
#
# NOTE: as-written exploratory script; not guaranteed turnkey. Paths are
# relative to the code/ directory.
# ============================================================================

library(ggplot2)
library(stringr)

fit.summary.tech <- read.csv('results/fit_summary_techreps.csv')
fit.all <- read.csv('results/fitness_all.csv')

mod <- glm(fitness~gene+type+treat, data=fit.all, )
summary(mod)
