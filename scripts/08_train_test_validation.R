#!/usr/bin/env Rscript
# Script: 08_train_test_validation.R
# Author: Smmrithi Ravindran
# Date: 2025
# Description: 70/30 train-test split validation with sex-balanced allocation

library(data.table)
library(dplyr)

setwd("/ix1/vgopalakrishnan/thesis_genomics_analysis")

# Load the 5-SNP data
geno_data <- fread("target_snps_raw_5snps.raw")

# Remove individuals with missing data
snp_cols <- c("rs2501352_A", "rs1426810_A", "rs10965183_A", "rs111339851_A", "rs7258841_A")
complete_data <- geno_data[complete.cases(geno_data[, ..snp_cols])]

cat(" TRAIN/TEST SPLIT VALIDATION \n\n")
cat("Total sample:", nrow(complete_data), "\n")
cat("Males:", sum(complete_data$SEX == 1), "\n")
cat("Females:", sum(complete_data$SEX == 2), "\n\n")

# 70/30 SPLIT WITH SEX-BALANCED ALLOCATION
set.seed(12345)

# Split males
males <- complete_data[SEX == 1]
n_males_train <- floor(0.7 * nrow(males))
male_train_idx <- sample(1:nrow(males), n_males_train)
males_train <- males[male_train_idx]
males_test <- males[-male_train_idx]

# Split females
females <- complete_data[SEX == 2]
n_females_train <- floor(0.7 * nrow(females))
female_train_idx <- sample(1:nrow(females), n_females_train)
females_train <- females[female_train_idx]
females_test <- females[-female_train_idx]

# Combine
train_data <- rbind(males_train, females_train)
test_data <- rbind(males_test, females_test)

cat(" SPLIT VERIFICATION \n")
cat("Training set:", nrow(train_data), 
    sprintf("(%.1f%% male)\n", sum(train_data$SEX==1)/nrow(train_data)*100))
cat("Testing set:", nrow(test_data), 
    sprintf("(%.1f%% male)\n\n", sum(test_data$SEX==1)/nrow(test_data)*100))

# CALCULATE EFFECT SIZES ON TRAINING DATA ONLY
cat(" TRAINING DATA: CALCULATING WEIGHTS \n\n")

weights <- c()
for(snp in snp_cols) {
  male_freq <- mean(train_data[SEX == 1, get(snp)], na.rm=TRUE) / 2
  female_freq <- mean(train_data[SEX == 2, get(snp)], na.rm=TRUE) / 2
  freq_diff <- female_freq - male_freq
  
  weights <- c(weights, freq_diff)
  
  cat(sprintf("%s: Male=%.4f, Female=%.4f, Diff=%+.4f\n", 
              snp, male_freq, female_freq, freq_diff))
}
names(weights) <- snp_cols

cat("\n")

# APPLY TO TEST DATA
cat(" APPLYING TO TEST DATA \n\n")

calculate_score <- function(data, weights, cols) {
  score_matrix <- as.matrix(data[, ..cols])
  return(as.vector(score_matrix %*% weights))
}

# Calculate scores
train_data[, Score := calculate_score(.SD, weights, snp_cols)]
test_data[, Score := calculate_score(.SD, weights, snp_cols)]

# Training results
cat("TRAINING SET:\n")
cat("Male mean:", round(mean(train_data[SEX==1, Score]), 4), "\n")
cat("Female mean:", round(mean(train_data[SEX==2, Score]), 4), "\n")
cat("Difference:", round(mean(train_data[SEX==2, Score]) - mean(train_data[SEX==1, Score]), 4), "\n")
cat("% Diff:", round((mean(train_data[SEX==2, Score]) - mean(train_data[SEX==1, Score])) / 
                      mean(train_data[SEX==1, Score]) * 100, 1), "%\n\n")

# Testing results
cat("TESTING SET:\n")
male_scores <- test_data[SEX == 1, Score]
female_scores <- test_data[SEX == 2, Score]

cat("Male mean:", round(mean(male_scores), 4), "\n")
cat("Female mean:", round(mean(female_scores), 4), "\n")
cat("Difference:", round(mean(female_scores) - mean(male_scores), 4), "\n")
cat("% Diff:", round((mean(female_scores) - mean(male_scores))/mean(male_scores)*100, 1), "%\n")

# Statistical test
t_result <- t.test(male_scores, female_scores)
cat("P-value:", format(t_result$p.value, scientific=TRUE, digits=3), "\n")

pooled_sd <- sqrt((sd(male_scores)^2 + sd(female_scores)^2) / 2)
cohens_d <- (mean(female_scores) - mean(male_scores)) / pooled_sd
cat("Cohen's d:", round(cohens_d, 3), "\n\n")

# SAVE RESULTS
fwrite(train_data, "train_data_70percent.csv")
fwrite(test_data, "test_data_30percent.csv")

summary_results <- data.table(
  Dataset = c("Training", "Testing"),
  N = c(nrow(train_data), nrow(test_data)),
  Male_Mean = c(mean(train_data[SEX==1, Score]), mean(male_scores)),
  Female_Mean = c(mean(train_data[SEX==2, Score]), mean(female_scores)),
  Difference = c(mean(train_data[SEX==2, Score]) - mean(train_data[SEX==1, Score]),
                 mean(female_scores) - mean(male_scores)),
  Pct_Diff = c((mean(train_data[SEX==2, Score]) - mean(train_data[SEX==1, Score])) / 
                mean(train_data[SEX==1, Score]) * 100,
               (mean(female_scores) - mean(male_scores))/mean(male_scores)*100)
)

fwrite(summary_results, "train_test_summary.csv")

cat("FILES SAVED \n")
cat("1. train_data_70percent.csv\n")
cat("2. test_data_30percent.csv\n")
cat("3. train_test_summary.csv\n")

cat("\n VALIDATION COMPLETE \n")
