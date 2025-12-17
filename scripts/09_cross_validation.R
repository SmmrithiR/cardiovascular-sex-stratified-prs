#!/usr/bin/env Rscript
# Script: 09_cross_validation.R
# Author: Smmrithi Ravindran
# Date: 2025
# Description: 10-fold cross-validation with sex-balanced fold assignment

library(data.table)
library(dplyr)

setwd("/ix1/vgopalakrishnan/thesis_genomics_analysis")

# Load data
geno_data <- fread("target_snps_raw_5snps.raw")

snp_cols <- c("rs2501352_A", "rs1426810_A", "rs10965183_A", "rs111339851_A", "rs7258841_A")
complete_data <- geno_data[complete.cases(geno_data[, ..snp_cols])]

cat(" 10-FOLD CROSS-VALIDATION \n\n")
cat("Total sample:", nrow(complete_data), "\n")
cat("Males:", sum(complete_data$SEX == 1), "\n")
cat("Females:", sum(complete_data$SEX == 2), "\n\n")

# CREATE 10 SEX-BALANCED FOLDS
set.seed(12345)

# Assign fold numbers separately by sex
males <- complete_data[SEX == 1]
females <- complete_data[SEX == 2]

males[, fold := sample(rep(1:10, length.out = .N))]
females[, fold := sample(rep(1:10, length.out = .N))]

# Combine
geno_data_folds <- rbind(males, females)

# Verify balance
cat(" FOLD BALANCE \n")
fold_summary <- geno_data_folds[, .(
  N = .N,
  N_Male = sum(SEX == 1),
  N_Female = sum(SEX == 2),
  Pct_Male = round(sum(SEX == 1)/.N * 100, 1)
), by = fold]
print(fold_summary)
cat("\n")

# RUN 10-FOLD CROSS-VALIDATION
cv_results <- list()

cat(" RUNNING CROSS-VALIDATION \n\n")

for(i in 1:10) {
  cat(sprintf("Fold %d/10...\n", i))
  
  # Split: fold i is test, rest are training
  train_data <- geno_data_folds[fold != i]
  test_data <- geno_data_folds[fold == i]
  
  # Calculate weights from training data
  weights <- c()
  for(snp in snp_cols) {
    male_freq <- mean(train_data[SEX == 1, get(snp)]) / 2
    female_freq <- mean(train_data[SEX == 2, get(snp)]) / 2
    weights <- c(weights, female_freq - male_freq)
  }
  names(weights) <- snp_cols
  
  # Calculate scores on test data
  calculate_score <- function(data, weights, cols) {
    score_matrix <- as.matrix(data[, ..cols])
    return(as.vector(score_matrix %*% weights))
  }
  
  test_data[, Score := calculate_score(.SD, weights, snp_cols)]
  
  # Statistics
  male_scores <- test_data[SEX == 1, Score]
  female_scores <- test_data[SEX == 2, Score]
  
  male_mean <- mean(male_scores)
  female_mean <- mean(female_scores)
  diff <- female_mean - male_mean
  pct_diff <- (diff / male_mean) * 100
  
  t_result <- t.test(male_scores, female_scores)
  p_value <- t_result$p.value
  
  pooled_sd <- sqrt((sd(male_scores)^2 + sd(female_scores)^2) / 2)
  cohens_d <- (female_mean - male_mean) / pooled_sd
  
  # Store results
  cv_results[[i]] <- data.table(
    Fold = i,
    N_Train = nrow(train_data),
    N_Test = nrow(test_data),
    Male_Mean = male_mean,
    Female_Mean = female_mean,
    Difference = diff,
    Pct_Difference = pct_diff,
    P_value = p_value,
    Cohens_d = cohens_d,
    Significant = ifelse(p_value < 0.05, "Yes", "No")
  )
  
  cat(sprintf("  Female %.1f%% higher (p=%.3f)\n", pct_diff, p_value))
}

cat("\n")

# AGGREGATE RESULTS
cv_results_df <- rbindlist(cv_results)

cat(" INDIVIDUAL FOLD RESULTS \n")
print(cv_results_df[, .(Fold, N_Test, Pct_Difference, P_value, Significant)])

cat("\n SUMMARY ACROSS ALL FOLDS \n")

mean_diff <- mean(cv_results_df$Pct_Difference)
sd_diff <- sd(cv_results_df$Pct_Difference)
se_diff <- sd_diff / sqrt(10)
ci_lower <- mean_diff - 1.96 * se_diff
ci_upper <- mean_diff + 1.96 * se_diff

cat(sprintf("Mean %% difference: %.2f%%\n", mean_diff))
cat(sprintf("SD: %.2f%%\n", sd_diff))
cat(sprintf("95%% CI: [%.2f%%, %.2f%%]\n", ci_lower, ci_upper))
cat(sprintf("Folds with p<0.05: %d/10\n", sum(cv_results_df$Significant == "Yes")))
cat(sprintf("Folds showing females higher: %d/10\n", sum(cv_results_df$Pct_Difference > 0)))
cat(sprintf("Mean Cohen's d: %.3f\n", mean(cv_results_df$Cohens_d)))

# SAVE RESULTS
fwrite(cv_results_df, "cv_10fold_detailed_results.csv")

summary_table <- data.table(
  Metric = c("Mean % Difference", "95% CI Lower", "95% CI Upper", 
             "Folds Significant", "Consistent Direction"),
  Value = c(sprintf("%.2f%%", mean_diff),
            sprintf("%.2f%%", ci_lower),
            sprintf("%.2f%%", ci_upper),
            sprintf("%d/10", sum(cv_results_df$Significant == "Yes")),
            sprintf("%d/10", sum(cv_results_df$Pct_Difference > 0)))
)

fwrite(summary_table, "cv_10fold_summary.csv")

cat("\n FILES SAVED \n")
cat("1. cv_10fold_detailed_results.csv\n")
cat("2. cv_10fold_summary.csv\n")

cat("\n CROSS-VALIDATION COMPLETE \n")
