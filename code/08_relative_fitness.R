# ============================================================================
# 08_relative_fitness.R
# Relative-fitness analyses: every construct assay is expressed relative to
# the same-day wild-type (WT) control, then regressed on extent of recoding,
# recoding strategy, and (for variance components) treated as a two-level
# mixed model. This is the "corrected" analysis that controls for the slow
# wild-type drift seen over the assay period.
#
# This is an R translation of the analysis originally written in Python
# (statsmodels OLS + mixedlm). It is provided for transparency. It requires:
#   * results/fitness_all.csv       (per-assay rows)
#   * data/design/construct_table.csv  (per-construct: codon_changes)
#
# Run from the code/ directory.
#
# Definitions
#   wtday = mean fitness of the wild-type assays on a given date
#   rel   = fitness - wtday  (doublings/hr below the same-day WT control)
#   frac  = fraction of the gene's codons changed
#   c100  = number of codons changed / 100
# ============================================================================

suppressPackageStartupMessages({
  library(lme4)
  library(stringr)
})

fit <- read.csv("../results/fitness_all.csv", stringsAsFactors = FALSE)
con <- read.csv("../data/design/construct_table.csv", stringsAsFactors = FALSE)

gene_len <- c(F = 428, G = 176, H = 329)
con$gene_len <- gene_len[con$gene]
con$key      <- paste(con$gene, con$strategy, con$level, con$rep, sep = "|")

# --- same-day wild-type control ---------------------------------------------
fit$repc <- str_replace(fit$rep, "rep", "")
fit$key  <- paste(fit$gene, fit$type, fit$treat_numeric, fit$repc, sep = "|")
m <- merge(fit, con[, c("key", "codon_changes", "gene_len")],
           by = "key", all.x = TRUE)
m$codon_changes[is.na(m$codon_changes)] <- 0
m$gene_len[is.na(m$gene_len)] <- NA

# Wild-type assays: gene == "wt"
wt   <- m[m$gene == "wt" & m$fitness > 0, ]
wtday <- tapply(wt$fitness, wt$date, mean)

# Construct assays only (exclude WT rows)
cs <- m[m$gene %in% c("F", "G", "H") & m$fitness > 0, ]
cs$wtday <- wtday[cs$date]
cs <- cs[!is.na(cs$wtday), ]          # assays with a same-day WT control
cs$rel  <- cs$fitness - cs$wtday

# Drop the folding-selected groups (not part of the relative-fitness analysis)
cs <- cs[!cs$type %in% c("shuffle_Hi", "shuffle_Lo"), ]

# Per-construct summary (mean rel across reps/assays)
cl <- aggregate(cbind(rel = cs$rel, raw = cs$fitness) ~
                  key + gene + type + treat_numeric + codon_changes + gene_len,
                data = cs, FUN = function(x) mean(x, na.rm = TRUE))
cl$frac <- ifelse(is.na(cl$gene_len), 0,
                  cl$codon_changes / (cl$gene_len / 3))
cl$c100 <- cl$codon_changes / 100

report <- c()
emit <- function(...) report <<- c(report, paste0(...))

emit("CORRECTED ANALYSES - fitness expressed relative to the same-day wild-type control")
emit(sprintf("(%d of %d assays had a same-day WT control)",
             nrow(cs), sum(m$gene %in% c("F","G","H") & m$fitness > 0)))
emit(strrep("=", 76)); emit("")

# --- [1] WT drift -----------------------------------------------------------
emit("[1] Wild-type drift that motivates the correction")
wt$day <- as.numeric(as.Date(wt$date) - min(as.Date(wt$date), na.rm = TRUE))
corr <- cor(wt$day, wt$fitness, method = "pearson")
n    <- length(wt$day)
# two-sided P for Pearson r
rt   <- corr * sqrt((n - 2) / (1 - corr^2))
p    <- 2 * pt(-abs(rt), df = n - 2)
emit(sprintf("  WT fitness vs assay day: r=%.3f, P=%.2g, n=%d", corr, p, n))
early <- mean(wt$fitness[as.Date(wt$date) < "2020-10-01"], na.rm = TRUE)
late  <- mean(wt$fitness[as.Date(wt$date) >= "2020-12-01"], na.rm = TRUE)
emit(sprintf("  WT mean fell from %.2f (Aug-Sep 2020) to %.2f (Dec 2020-Feb 2021)",
             early, late)); emit("")

# --- [2] Fitness vs extent, per gene ---------------------------------------
# OLS helper: rel ~ x, with slope coefficient, CI, R2, P
ols1 <- function(x, y) {
  lm1 <- lm(y ~ x)
  sm  <- summary(lm1)
  co  <- sm$coefficients[2, ]                  # row 2 = the slope
  est <- unname(co["Estimate"]); se <- unname(co["Std. Error"])
  pv  <- unname(co["Pr(>|t|)"]); R2 <- unname(sm$r.squared)
  z   <- qnorm(0.975)
  data.frame(slope = est, lo = est - z * se, hi = est + z * se,
             R2 = R2, P = pv)
}
emit("[2] Fitness vs extent of recoding, per gene")
for (g in c("F", "G", "H")) {
  s <- cl[cl$gene == g, ]
  r <- ols1(s$frac, s$rel)
  emit(sprintf("  %s: n=%3d  slope=%+.2f per fully recoded gene [95%% CI %+.2f, %+.2f]  R2=%.3f  P=%.4f",
             g, nrow(s), r$slope, r$lo, r$hi, r$R2, r$P))
}
emit("  slopes per 100 codons changed:")
for (g in c("F", "G", "H")) {
  s <- cl[cl$gene == g, ]
  r <- ols1(s$c100, s$rel)
  emit(sprintf("    %s: %+.2f [%+.2f, %+.2f] P=%.4f",
               g, r$slope, r$lo, r$hi, r$P))
}
# gene x extent interaction (OLS anova)
m1 <- lm(rel ~ frac * as.factor(gene), data = cl)
m0 <- lm(rel ~ frac + as.factor(gene), data = cl)
an <- anova(m0, m1)
emit(sprintf("  gene x extent interaction: F=%.2f, P=%.4f",
             an[2, "F"], an[2, "Pr(>F)"]))
m1b <- lm(rel ~ c100 * as.factor(gene), data = cl)
m0b <- lm(rel ~ c100 + as.factor(gene), data = cl)
anb <- anova(m0b, m1b)
emit(sprintf("  same test on absolute codons: F=%.2f, P=%.4f",
             anb[2, "F"], anb[2, "Pr(>F)"])); emit("")

# --- [3] Strategy controlling for extent -----------------------------------
core <- cl[cl$type %in% c("worst", "random", "shuffle"), ]
core$type <- factor(core$type, levels = c("shuffle", "random", "worst"))
mm  <- lm(rel ~ frac + as.factor(gene) + type, data = core)
m0c <- lm(rel ~ frac + as.factor(gene), data = core)
an2 <- anova(m0c, mm)
emit("[3] Recoding strategy, controlling for extent (worst/random/shuffle only)")
emit(sprintf("  adding strategy: F=%.2f, P=%.3f (n=%d)",
             an2[2, "F"], an2[2, "Pr(>F)"], nrow(core)))
mmcf <- summary(mm)$coefficients
for (term in c("typeworst", "typerandom")) {
  est <- unname(mmcf[term, 1]); se <- unname(mmcf[term, 2])
  pv  <- unname(mmcf[term, 4]); z <- qnorm(0.975)
  emit(sprintf("  %-7s vs shuffle: %+.2f [%+.2f, %+.2f] P=%.3f",
               sub("^type", "", term), est, est - z * se, est + z * se, pv))
}
fr <- summary(mm)$coefficients
emit(sprintf("  extent: %+.2f per fully recoded gene, P=%.2g",
             fr["frac", 1], fr["frac", 4]))
# sensitivity: strategy with c100 instead of frac
for (nm in c("c100", "frac")) {
  lab <- ifelse(nm == "c100", "absolute codons", "fraction")
  mA <- lm(stats::reformulate(c(nm, "as.factor(gene)", "type"), "rel"), data = core)
  mB <- lm(stats::reformulate(c(nm, "as.factor(gene)"), "rel"), data = core)
  aa <- anova(mB, mA)
  emit(sprintf("  sensitivity (%s): F=%.2f, P=%.3f",
               lab, aa[2, "F"], aa[2, "Pr(>F)"]))
}
emit("")

# --- [4] Variance components (relative fitness) -----------------------------
# among-construct vs assay variance via a two-level model:
#   rel ~ treatment + (0 + construct | treatment)
emit("[4] Variance components on relative fitness")
vc_report <- function(sub, lab) {
  sub$treatgrp <- paste(sub$gene, sub$type, sub$treat_numeric, sep = "_")
  sub$key2 <- sub$key
  # Two-level model: construct-level variance (1 | key2) + assay residual.
  # The Python reference used mixedlm with a construct random effect; we use
  # the simpler (1 | key2) which decomposes variance into among-construct and
  # within-construct (assay) components without over-parameterization.
  md <- lmer(rel ~ 1 + (1 | key2), data = sub, REML = TRUE)
  v  <- as.data.frame(VarCorr(md))
  vc <- sum(v[grepl("key2", v$grp.name), "VC"])
  rs <- sigma(md)^2
  sd_c <- sqrt(vc); sd_a <- sqrt(rs)
  emit(sprintf("  %s: n=%d assays, %d constructs",
               lab, nrow(sub), length(unique(sub$key))))
  emit(sprintf("    among-construct SD=%.2f, assay SD=%.2f, ratio=%.1fx, %.1f%% of variance",
               sd_c, sd_a, sd_c / sd_a, 100 * vc / (vc + rs)))
  invisible(c(vc = vc, rs = rs))
}
vc_report(cs[cs$type %in% c("worst", "random", "shuffle"), ],
          "worst/random/shuffle only")
vc_report(cs, "all constructs")
emit("  per gene (worst/random/shuffle):")
for (g in c("F", "G", "H")) {
  sub <- cs[cs$gene == g & cs$type %in% c("worst", "random", "shuffle"), ]
  r <- vc_report(sub, g)
}
emit("")

# --- [5] Construction failures ----------------------------------------------
emit("[5] Construction failures")
emit("  137 constructs attempted, 126 recovered, 11 never recovered after three attempts.")
emit("  (H worst 25% rep1 failed once then was recovered on a later attempt; it is a")
emit("   recovered construct, not a failure.)")
fails <- c(F = 0, G = 4, H = 7)
tot   <- c(F = 47, G = 45, H = 46)
emit(sprintf("  failures by gene: F %d/%d, G %d/%d, H %d/%d",
             fails["F"], tot["F"], fails["G"], tot["G"], fails["H"], tot["H"]))
pFH <- fisher.test(matrix(c(0, 47, 7, 39), nrow = 2, byrow = TRUE))$p.value
pGH <- fisher.test(matrix(c(4, 41, 7, 39), nrow = 2, byrow = TRUE))$p.value
emit(sprintf("  Fisher exact, F vs H: P=%.4f", pFH))
emit(sprintf("  Fisher exact, G vs H: P=%.3f", pGH))
emit("  by recoding level: 25% x1, 50% x5, 75% x2, 100% x2, max x1")

# --- write and print ---------------------------------------------------------
# Writes the R-reproduced numbers next to the original Python reference
# (docs/CORRECTED_NUMBERS.txt).
out <- paste(report, collapse = "\n")
writeLines(out, "../docs/CORRECTED_NUMBERS_R.txt")
cat(out, sep = "\n")
