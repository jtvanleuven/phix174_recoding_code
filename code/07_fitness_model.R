# ============================================================================
# 07_fitness_model.R
# PRIMARY model: raw fitness, wild-type assays included, assay date as a
# random effect. This is the model from which all headline statistics in the
# manuscript come (see docs/MODEL_NUMBERS.txt for the reference values).
#
# This is an R translation of the analysis originally written in Python
# (statsmodels mixedlm). It is provided for transparency. It requires:
#   * results/fitness_all.csv       (per-assay rows)
#   * data/design/construct_table.csv  (per-construct: codon_changes)
#
# Run from the code/ directory.
#
# Model (restricted maximum likelihood, mixed-effects):
#   fitness ~ day + lu + (gene intercepts) + (per-gene extent slope)
#             + strategy contrasts (shuffle = reference)
#             + (1 | date)
# where
#   day          = days since the first assay date (slow wild-type drift)
#   lu           = 1 for the "lu" technician, 0 for "emma"
#   g_g          = indicator for gene g in {F,G,H}  (gene intercept)
#   f_g          = fraction of gene g's codons changed  (extent slope)
#   s_random     = indicator: recoding strategy == "random"
#   s_worst      = indicator: recoding strategy == "worst"
#   (1 | date)   = per-assay-date random intercept
# ============================================================================

suppressPackageStartupMessages({
  library(lme4)        # mixed-effects models
  library(stringr)
})

# --- locate inputs (relative to code/) --------------------------------------
fit  <- read.csv("../results/fitness_all.csv", stringsAsFactors = FALSE)
con  <- read.csv("../data/design/construct_table.csv", stringsAsFactors = FALSE)

# gene lengths in codons (F=428, G=176, H=329)
gene_len <- c(F = 428, G = 176, H = 329)
con$gene_len <- gene_len[con$gene]
con$key      <- paste(con$gene, con$strategy, con$level, con$rep, sep = "|")

# --- build the model data frame ---------------------------------------------
# Map each assay to its construct's codon-change count. Key = gene|type|level|rep.
fit$repc <- str_replace(fit$rep, "rep", "")
fit$key  <- paste(fit$gene, fit$type, fit$treat_numeric, fit$repc, sep = "|")

m <- merge(fit, con[, c("key", "codon_changes", "gene_len")],
           by = "key", all.x = TRUE)

# Wild-type assays have no construct (gene == "wt"); give them zero extent.
m$codon_changes[is.na(m$codon_changes)] <- 0
m$gene_len[is.na(m$gene_len)] <- NA

# Exclude construction failures (fitness == 0 placeholder) and the two
# folding-selected groups that the headline model is not fitted to.
m <- m[m$fitness > 0, ]
m <- m[!m$type %in% c("shuffle_Hi", "shuffle_Lo"), ]

# Numeric covariates (wt assays: frac = 0, gene_len NA -> frac NA -> coerce to 0)
m$frac   <- m$codon_changes / (m$gene_len / 3)
m$frac   <- ifelse(is.na(m$frac), 0, m$frac)
m$day    <- as.numeric(as.Date(m$date) - min(as.Date(m$date), na.rm = TRUE))
m$lu     <- as.numeric(m$sci == "lu")
m$s_random <- as.numeric(m$type == "random")
m$s_worst  <- as.numeric(m$type == "worst")

# Per-gene extent terms: fraction of that gene's codons changed
m$f_F <- m$frac * (m$gene == "F")
m$f_G <- m$frac * (m$gene == "G")
m$f_H <- m$frac * (m$gene == "H")

# Per-gene "number of codons changed per 100" terms (for the per-100-codon view)
m$c_F <- (m$codon_changes / 100) * (m$gene == "F")
m$c_G <- (m$codon_changes / 100) * (m$gene == "G")
m$c_H <- (m$codon_changes / 100) * (m$gene == "H")

# Gene intercept terms (wild-type rows have gene == "wt", so these are 0)
m$g_F <- as.numeric(m$gene == "F")
m$g_G <- as.numeric(m$gene == "G")
m$g_H <- as.numeric(m$gene == "H")

# The model is fitted to the three recoding strategies together with the
# wild-type assays (wt rows act as the gene=0, extent=0 reference).
core <- m[m$type %in% c("shuffle", "random", "worst", "wt"), ]
wt   <- core[core$gene == "wt", ]

report <- c()
emit <- function(...) report <<- c(report, paste0(...))

emit("PRIMARY MODEL - raw fitness, wild-type assays included, assay date as a random effect")
emit(strrep("=", 84)); emit("")
emit(sprintf("  %d construct assays + %d wild-type assays, %d dates",
             sum(core$gene != "wt"), nrow(wt), length(unique(core$date))))

# In R, `formula + (1|g)` is not a valid lmer call, so we paste the RE in.
# Collapse deparse to one line so a long formula does not wrap and corrupt.
with_re <- function(f) as.formula(paste(paste(deparse(f), collapse = " "),
                                        "(1 | date)", sep = " + "))

# Fixed-effect formula (gene intercepts + per-gene extent slopes + strategy)
fmla   <- fitness ~ day + lu + g_F + g_G + g_H + f_F + f_G + f_H + s_random + s_worst
fmla_c <- fitness ~ day + lu + g_F + g_G + g_H + c_F + c_G + c_H + s_random + s_worst

# --- fit the primary model (RE on date) -------------------------------------
mm <- lmer(with_re(fmla),  data = core, REML = TRUE)

# helper to pull a coefficient, Wald P, and 95% CI.
# lmer summary() gives Estimate / Std. Error / t value (no P column),
# so the two-sided Wald P is computed from the t value.
co <- function(fit_obj, term) {
  b    <- summary(fit_obj)$coefficients
  est  <- unname(b[term, "Estimate"])
  se   <- unname(b[term, "Std. Error"])
  tval <- unname(b[term, "t value"])
  data.frame(Estimate = est, SE = se, P = 2 * pt(-abs(tval), df = Inf))
}
ci95 <- function(x) x$Estimate + c(-1, 1) * qnorm(0.975) * x$SE

emit(""); emit("[1] Technician")
r <- co(mm, "lu"); lo <- ci95(r)
emit(sprintf("  Lu - Emma: %+.2f [%+.2f, %+.2f], P=%.3f",
             r$Estimate, lo[1], lo[2], r$P))

# technician x construct interaction (the test for a strategy-confounded effect)
core$luC <- core$lu * (core$gene != "wt")
m5 <- lmer(with_re(update(fmla, . ~ . + luC)), data = core, REML = TRUE)
r <- co(m5, "luC"); lo <- ci95(r)
emit(sprintf("  technician x construct: %+.2f [%+.2f, %+.2f], P=%.3f",
             r$Estimate, lo[1], lo[2], r$P))

emit(""); emit("[2] Recoding strategy (reference = shuffle)")
for (nm in c("s_worst", "s_random")) {
  lab <- c(s_worst = "worst", s_random = "random")[[nm]]
  r <- co(mm, nm); lo <- ci95(r)
  emit(sprintf("  %-7s vs shuffle: %+.2f [%+.2f, %+.2f], P=%.3f",
               lab, r$Estimate, lo[1], lo[2], r$P))
}
fmla_nolu <- fitness ~ day + g_F + g_G + g_H + f_F + f_G + f_H + s_random + s_worst
m3 <- lmer(with_re(fmla_nolu), data = core, REML = TRUE)
mw <- co(m3, "s_worst")
emit(sprintf("  (technician omitted: worst %+.2f, P=%.3f)", mw$Estimate, mw$P))

emit(""); emit("[3] Extent of recoding, per gene")
emit("  per whole gene:")
for (nm in c("f_F", "f_G", "f_H")) {
  lab <- sub("^f_", "", nm)
  r <- co(mm, nm); lo <- ci95(r)
  emit(sprintf("    %s: %+.2f [%+.2f, %+.2f], P=%.2g",
               lab, r$Estimate, lo[1], lo[2], r$P))
}
mc <- lmer(with_re(fmla_c), data = core, REML = TRUE)
emit("  per 100 codons changed:")
for (nm in c("c_F", "c_G", "c_H")) {
  lab <- sub("^c_", "", nm)
  r <- co(mc, nm); lo <- ci95(r)
  emit(sprintf("    %s: %+.2f [%+.2f, %+.2f], P=%.2g",
               lab, r$Estimate, lo[1], lo[2], r$P))
}

# gene x extent interaction: likelihood-ratio test vs a single pooled slope
core$cc <- core$codon_changes / 100
m0 <- lmer(with_re(fitness ~ day + lu + g_F + g_G + g_H + cc + s_random + s_worst),
           data = core, REML = FALSE)
lr <- 2 * (logLik(mc) - logLik(m0))
p  <- pchisq(lr, df = 2, lower.tail = FALSE)
emit(sprintf("  gene x extent: likelihood ratio %.2f on 2 df, P=%.2g", lr, p))

emit(""); emit("[4] Construct-level values for tables (model-adjusted: fitness")
emit("    minus the fitted date effect)")
# With a per-date RE the date random effects are the BLUPs; on a singular
# fit the date variance is ~0 so the BLUPs are ~0. Fall back to 0 if absent.
# On a singular fit the date variance is ~0 and lme4 returns only the
# variance component (no per-date BLUPs), so the date adjustment is ~0.
if ("vcORE" %in% colnames(as.data.frame(VarCorr(mm)$date[, , drop = TRUE]))) {
  re  <- as.data.frame(VarCorr(mm)$date[, , drop = TRUE])
  core$dateRE <- re$vcORE[match(core$date, rownames(re))]
} else {
  core$dateRE <- 0
}
core$adj <- core$fitness - core$dateRE - summary(mm)$coefficients["day", 1] * core$day
           - summary(mm)$coefficients["lu",   1] * core$lu
wtadj <- mean(core$adj[core$gene == "wt"], na.rm = TRUE)
emit(sprintf("  wild-type baseline after adjustment: %.2f doublings/hr", wtadj))

# --- write and print ---------------------------------------------------------
# Writes the R-reproduced numbers next to the original Python reference
# (docs/MODEL_NUMBERS.txt). The two should agree on all headline statistics.
out <- paste(report, collapse = "\n")
writeLines(out, "../docs/MODEL_NUMBERS_R.txt")
cat(out, sep = "\n")
