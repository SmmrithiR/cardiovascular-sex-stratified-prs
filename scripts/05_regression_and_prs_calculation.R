#!/usr/bin/env Rscript
# Script: 05_prs_calculation_corrected.R
# Author: Smmrithi Ravindran
# Description: Calculate PRS using frequency-difference weights and classify
#              risk using OVERALL percentile thresholds (one shared distribution).
#              Replaces the previous version which had a bug in sex-specific
#              percentile assignment using ifelse() on a data.table, and used
#              incorrect 80th/20th cutoffs instead of 95th/75th/25th/5th.

library(data.table)

setwd("/ix1/vgopalakrishnan/thesis_genomics_analysis")

# ── Load data ──────────────────────────────────────────────────────────────────
geno_data <- fread("target_snps_raw_5snps.raw")

snp_cols <- c("rs2501352_A", "rs1426810_A", "rs10965183_A",
              "rs111339851_A", "rs7258841_A")

# ── Remove individuals with missing genotype data ──────────────────────────────
complete_data <- geno_data[complete.cases(geno_data[, ..snp_cols])]

cat("=== SAMPLE FLOW ===\n")
cat("Raw total:            ", nrow(geno_data),     "(should be 2524)\n")
cat("After SNP exclusions: ", nrow(complete_data), "(should be 2509)\n")
cat("  Male:               ", sum(complete_data$SEX == 1), "(should be 1237)\n")
cat("  Female:             ", sum(complete_data$SEX == 2), "(should be 1272)\n\n")

# ── Calculate frequency-difference weights from full dataset ───────────────────
cat("=== FREQUENCY-DIFFERENCE WEIGHTS ===\n")
weights <- c()
for (snp in snp_cols) {
  male_freq   <- mean(complete_data[SEX == 1, get(snp)], na.rm = TRUE) / 2
  female_freq <- mean(complete_data[SEX == 2, get(snp)], na.rm = TRUE) / 2
  freq_diff   <- female_freq - male_freq
  weights     <- c(weights, freq_diff)
  cat(sprintf("%s: male=%.4f, female=%.4f, diff=%+.4f\n",
              snp, male_freq, female_freq, freq_diff))
}
names(weights) <- snp_cols

# ── Calculate PRS ──────────────────────────────────────────────────────────────
score_matrix <- as.matrix(complete_data[, ..snp_cols])
complete_data[, PRS := as.vector(score_matrix %*% weights)]

# ── PRS summary by sex ─────────────────────────────────────────────────────────
cat("\n=== PRS SUMMARY BY SEX ===\n")
male_prs   <- complete_data[SEX == 1, PRS]
female_prs <- complete_data[SEX == 2, PRS]

male_mean   <- mean(male_prs)
female_mean <- mean(female_prs)
pct_diff    <- (female_mean - male_mean) / male_mean * 100

male_se   <- sd(male_prs)   / sqrt(length(male_prs))
female_se <- sd(female_prs) / sqrt(length(female_prs))

cat(sprintf("Male mean PRS:   %.4f (95%% CI: %.3f-%.3f, SD=%.3f)\n",
            male_mean,
            male_mean - 1.96 * male_se,
            male_mean + 1.96 * male_se,
            sd(male_prs)))
cat(sprintf("Female mean PRS: %.4f (95%% CI: %.3f-%.3f, SD=%.3f)\n",
            female_mean,
            female_mean - 1.96 * female_se,
            female_mean + 1.96 * female_se,
            sd(female_prs)))
cat(sprintf("Sex difference:  %.1f%%\n\n", pct_diff))

# ── t-test ─────────────────────────────────────────────────────────────────────
t_result <- t.test(male_prs, female_prs)
pooled_sd <- sqrt(((length(male_prs)-1)*var(male_prs) +
                   (length(female_prs)-1)*var(female_prs)) /
                  (length(male_prs) + length(female_prs) - 2))
cohens_d  <- (female_mean - male_mean) / pooled_sd

cat(sprintf("t = %.2f, df = %d, p = %.2e, Cohen's d = %.2f\n\n",
            t_result$statistic,
            round(t_result$parameter),
            t_result$p.value,
            cohens_d))

# ── OVERALL percentile thresholds (one shared distribution) ───────────────────
# NOTE: We use OVERALL percentiles — a single threshold applied to all
# individuals regardless of sex. This is the correct approach for asking
# "are women more likely to be flagged as high risk by a shared clinical
# threshold?" and is what is described in the manuscript Methods section.
# The previous script incorrectly used sex-specific percentiles via a buggy
# ifelse() call, and also used non-standard 80th/20th cutoffs.

cat("=== RISK STRATIFICATION (OVERALL PERCENTILES) ===\n")

p95 <- quantile(complete_data$PRS, 0.95)
p75 <- quantile(complete_data$PRS, 0.75)
p25 <- quantile(complete_data$PRS, 0.25)
p05 <- quantile(complete_data$PRS, 0.05)

cat(sprintf("Thresholds: <5th=%.4f, 25th=%.4f, 75th=%.4f, 95th=%.4f\n\n",
            p05, p25, p75, p95))

complete_data[, Risk_Category := fcase(
  PRS >= p95, "Very High Risk",
  PRS >= p75, "High Risk",
  PRS >= p25, "Intermediate Risk",
  PRS >= p05, "Low Risk",
  default =   "Very Low Risk"
)]

# ── Table 4 ───────────────────────────────────────────────────────────────────
table4      <- complete_data[, .N, by = .(SEX, Risk_Category)]
table4_wide <- dcast(table4, Risk_Category ~ SEX, value.var = "N", fill = 0)
setnames(table4_wide, c("1", "2"), c("Male_N", "Female_N"))

male_total   <- sum(complete_data$SEX == 1)
female_total <- sum(complete_data$SEX == 2)

table4_wide[, Male_Pct   := round(Male_N   / male_total   * 100, 1)]
table4_wide[, Female_Pct := round(Female_N / female_total * 100, 1)]

# Relative difference = (female rate - male rate) / male rate * 100
table4_wide[, Rel_Diff := round(
  (Female_N / female_total - Male_N / male_total) /
  (Male_N / male_total) * 100, 1)]

cat_order   <- c("Very High Risk", "High Risk", "Intermediate Risk",
                 "Low Risk", "Very Low Risk")
table4_wide <- table4_wide[match(cat_order, Risk_Category)]

cat("=== TABLE 4 ===\n")
print(table4_wide)

# High Risk Total row
high_m <- sum(table4_wide[Risk_Category %in% c("Very High Risk","High Risk"), Male_N])
high_f <- sum(table4_wide[Risk_Category %in% c("Very High Risk","High Risk"), Female_N])
cat(sprintf("\nHigh Risk Total (>=75th): Male %d (%.1f%%), Female %d (%.1f%%), Rel Diff %.1f%%\n",
            high_m, high_m/male_total*100,
            high_f, high_f/female_total*100,
            (high_f/female_total - high_m/male_total)/(high_m/male_total)*100))

# Verify totals
cat(sprintf("\nMale total:   %d (should be %d) %s\n",
            sum(table4_wide$Male_N), male_total,
            ifelse(sum(table4_wide$Male_N)==male_total, "✓", "ERROR")))
cat(sprintf("Female total: %d (should be %d) %s\n",
            sum(table4_wide$Female_N), female_total,
            ifelse(sum(table4_wide$Female_N)==female_total, "✓", "ERROR")))

# ── Key manuscript figures verification ───────────────────────────────────────
cat("\n=== MANUSCRIPT FIGURES VERIFICATION ===\n")
vh_m <- table4_wide[Risk_Category == "Very High Risk", Male_N]
vh_f <- table4_wide[Risk_Category == "Very High Risk", Female_N]
vh_rel_diff <- (vh_f/female_total - vh_m/male_total) / (vh_m/male_total) * 100
per_1000 <- (vh_f/female_total - vh_m/male_total) * 1000

cat(sprintf("Very High Risk males:   %d (%.1f%%)\n", vh_m, vh_m/male_total*100))
cat(sprintf("Very High Risk females: %d (%.1f%%)\n", vh_f, vh_f/female_total*100))
cat(sprintf("Relative difference:    %.1f%% (manuscript reports 71%%)\n", vh_rel_diff))
cat(sprintf("Per 1,000 screened:     %.0f additional females\n", per_1000))

# ── Save outputs ───────────────────────────────────────────────────────────────
fwrite(complete_data, "results/prs_scores_final.csv")
fwrite(table4_wide,   "results/table4_corrected.csv")

cat("\n=== SAVED ===\n")
cat("results/prs_scores_final.csv\n")
cat("results/table4_corrected.csv\n")
cat("\nAnalysis complete.\n")
