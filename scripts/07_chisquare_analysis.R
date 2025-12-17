#!/usr/bin/env Rscript
# Script: 07_chisquare_analysis.R
# Author: Smmrithi Ravindran
# Date: 2025
# Description: Chi-square tests for sex-differential variants (replacing linear regression)

library(data.table)
library(dplyr)

setwd("/ix1/vgopalakrishnan/thesis_genomics_analysis")

# Load the 5-SNP data
geno_data <- fread("target_snps_raw_5snps.raw")

# Remove individuals with missing data
snp_cols <- c("rs2501352_A", "rs1426810_A", "rs10965183_A", "rs111339851_A", "rs7258841_A")
complete_data <- geno_data[complete.cases(geno_data[, ..snp_cols])]

cat(" CHI-SQUARE ANALYSIS FOR SEX DIFFERENCES \n\n")
cat("Sample size:", nrow(complete_data), "\n")
cat("Males:", sum(complete_data$SEX == 1), "\n")
cat("Females:", sum(complete_data$SEX == 2), "\n\n")

# FULL DATASET ANALYSIS
cat(" FULL DATASET: CHI-SQUARE TESTS \n\n")

results_full <- data.table()

for(snp in snp_cols) {
  
  # Calculate allele frequencies
  male_genotypes <- complete_data[SEX == 1, get(snp)]
  female_genotypes <- complete_data[SEX == 2, get(snp)]
  
  # Convert genotypes to allele counts
  male_allele_count <- sum(male_genotypes, na.rm = TRUE)
  male_total_alleles <- length(male_genotypes) * 2
  male_freq <- male_allele_count / male_total_alleles
  
  female_allele_count <- sum(female_genotypes, na.rm = TRUE)
  female_total_alleles <- length(female_genotypes) * 2
  female_freq <- female_allele_count / female_total_alleles
  
  # Create 2x2 contingency table
  contingency_table <- matrix(
    c(male_total_alleles - male_allele_count, male_allele_count,
      female_total_alleles - female_allele_count, female_allele_count),
    nrow = 2,
    byrow = TRUE
  )
  
  # Chi-square test
  chisq_result <- chisq.test(contingency_table, correct = FALSE)
  
  # Calculate odds ratio
  a <- contingency_table[1,1]
  b <- contingency_table[1,2]
  c <- contingency_table[2,1]
  d <- contingency_table[2,2]
  
  odds_ratio <- (a * d) / (b * c)
  log_or <- log(odds_ratio)
  
  # Frequency difference
  freq_diff <- female_freq - male_freq
  abs_freq_diff <- abs(freq_diff)
  
  # Store results
  results_full <- rbind(results_full, data.table(
    SNP = snp,
    Male_Freq = male_freq,
    Female_Freq = female_freq,
    Freq_Difference = freq_diff,
    Abs_Freq_Diff = abs_freq_diff,
    Chi_Square = chisq_result$statistic,
    P_Value = chisq_result$p.value,
    Odds_Ratio = odds_ratio,
    Log_Odds_Ratio = log_or,
    Direction = ifelse(freq_diff > 0, "Higher in Females", "Higher in Males")
  ))
  
  # Print detailed results
  cat(sprintf("%s:\n", snp))
  cat(sprintf("  Male frequency: %.4f (%.1f%%)\n", male_freq, male_freq*100))
  cat(sprintf("  Female frequency: %.4f (%.1f%%)\n", female_freq, female_freq*100))
  cat(sprintf("  Difference: %.4f (%.1f%%)\n", freq_diff, freq_diff*100))
  cat(sprintf("  Chi-square: %.4f\n", chisq_result$statistic))
  cat(sprintf("  P-value: %.4e\n", chisq_result$p.value))
  cat(sprintf("  Odds Ratio: %.4f\n", odds_ratio))
  cat(sprintf("  Direction: %s\n\n", ifelse(freq_diff > 0, "Higher in Females", "Higher in Males")))
}

# Print summary table
cat("\n SUMMARY TABLE (FULL DATASET) \n")
print(results_full[, .(SNP, Male_Freq, Female_Freq, Freq_Difference, 
                       Chi_Square, P_Value, Odds_Ratio, Direction)])

# Save results
fwrite(results_full, "chisquare_results_full_dataset.csv")

# TRAIN/TEST SPLIT WITH CHI-SQUARE
cat("\n\n TRAIN/TEST SPLIT ANALYSIS \n\n")

set.seed(12345)

# Split by sex
males <- complete_data[SEX == 1]
females <- complete_data[SEX == 2]

n_males_train <- floor(0.7 * nrow(males))
male_train_idx <- sample(1:nrow(males), n_males_train)
males_train <- males[male_train_idx]
males_test <- males[-male_train_idx]

n_females_train <- floor(0.7 * nrow(females))
female_train_idx <- sample(1:nrow(females), n_females_train)
females_train <- females[female_train_idx]
females_test <- females[-female_train_idx]

train_data <- rbind(males_train, females_train)
test_data <- rbind(males_test, females_test)

cat("Training set:", nrow(train_data), "\n")
cat("Testing set:", nrow(test_data), "\n\n")

# Calculate statistics on TRAINING data
cat(" TRAINING DATA: CHI-SQUARE TESTS \n\n")

results_train <- data.table()

for(snp in snp_cols) {
  
  male_genotypes <- train_data[SEX == 1, get(snp)]
  female_genotypes <- train_data[SEX == 2, get(snp)]
  
  male_allele_count <- sum(male_genotypes, na.rm = TRUE)
  male_total_alleles <- length(male_genotypes) * 2
  male_freq <- male_allele_count / male_total_alleles
  
  female_allele_count <- sum(female_genotypes, na.rm = TRUE)
  female_total_alleles <- length(female_genotypes) * 2
  female_freq <- female_allele_count / female_total_alleles
  
  contingency_table <- matrix(
    c(male_total_alleles - male_allele_count, male_allele_count,
      female_total_alleles - female_allele_count, female_allele_count),
    nrow = 2,
    byrow = TRUE
  )
  
  chisq_result <- chisq.test(contingency_table, correct = FALSE)
  
  a <- contingency_table[1,1]
  b <- contingency_table[1,2]
  c <- contingency_table[2,1]
  d <- contingency_table[2,2]
  
  odds_ratio <- (a * d) / (b * c)
  
  freq_diff <- female_freq - male_freq
  
  results_train <- rbind(results_train, data.table(
    SNP = snp,
    Male_Freq = male_freq,
    Female_Freq = female_freq,
    Freq_Difference = freq_diff,
    Chi_Square = chisq_result$statistic,
    P_Value = chisq_result$p.value,
    Odds_Ratio = odds_ratio
  ))
  
  cat(sprintf("%s: Δf=%.4f, χ²=%.2f, p=%.4e, OR=%.4f\n", 
              snp, freq_diff, chisq_result$statistic, 
              chisq_result$p.value, odds_ratio))
}

cat("\n")

# Calculate genetic scores using frequency differences as weights
cat(" CALCULATING GENETIC SCORES \n\n")

# Use frequency differences as weights
weights <- results_train$Freq_Difference
names(weights) <- snp_cols

# Calculate scores on test data
calculate_score <- function(data, weights, cols) {
  score_matrix <- as.matrix(data[, ..cols])
  scores <- as.vector(score_matrix %*% weights)
  return(scores)
}

test_data[, Genetic_Score := calculate_score(.SD, weights, snp_cols)]

# Test for differences
male_scores <- test_data[SEX == 1, Genetic_Score]
female_scores <- test_data[SEX == 2, Genetic_Score]

cat("TEST DATA RESULTS \n")
cat("Male mean score:", round(mean(male_scores), 4), "\n")
cat("Female mean score:", round(mean(female_scores), 4), "\n")
cat("Difference:", round(mean(female_scores) - mean(male_scores), 4), "\n")
cat("% Difference:", round((mean(female_scores) - mean(male_scores))/mean(male_scores)*100, 1), "%\n")

t_result <- t.test(male_scores, female_scores)
cat("T-test p-value:", format(t_result$p.value, scientific = TRUE, digits = 3), "\n")

pooled_sd <- sqrt((sd(male_scores)^2 + sd(female_scores)^2) / 2)
cohens_d <- (mean(female_scores) - mean(male_scores)) / pooled_sd
cat("Cohen's d:", round(cohens_d, 3), "\n")

# Save results
fwrite(results_train, "chisquare_results_training.csv")
fwrite(test_data, "test_data_with_scores.csv")

# CREATE PUBLICATION TABLE
cat("\n CREATING PUBLICATION TABLE \n")

pub_table <- results_full[, .(
  SNP = SNP,
  Male_Frequency = sprintf("%.3f", Male_Freq),
  Female_Frequency = sprintf("%.3f", Female_Freq),
  Difference = sprintf("%+.3f", Freq_Difference),
  Chi_Square = sprintf("%.2f", Chi_Square),
  P_Value = ifelse(P_Value < 0.001, 
                   sprintf("%.2e", P_Value),
                   sprintf("%.3f", P_Value)),
  Odds_Ratio = sprintf("%.3f", Odds_Ratio),
  Direction = Direction
)]

print(pub_table)
fwrite(pub_table, "chisquare_publication_table.csv")

cat("\n FILES SAVED \n")
cat("1. chisquare_results_full_dataset.csv\n")
cat("2. chisquare_results_training.csv\n")
cat("3. test_data_with_scores.csv\n")
cat("4. chisquare_publication_table.csv\n")

cat("\n ANALYSIS COMPLETE \n")
