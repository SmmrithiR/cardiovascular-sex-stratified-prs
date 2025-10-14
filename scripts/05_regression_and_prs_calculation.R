#!/usr/bin/env Rscript
# Script: 05_regression_and_prs_calculation.R
# Author: Smmrithi Ravindran
# Description: Calculate proper effect sizes using regression and compute PRS

# Load required libraries
library(data.table)
library(dplyr)

# Set working directory
setwd("/ix1/vgopalakrishnan/thesis_genomics_analysis")

# Load the genotype data
geno_data <- fread("target_snps_raw.raw")

# Check dimensions
cat("Dataset dimensions:", dim(geno_data), "\n")
cat("Number of individuals:", nrow(geno_data), "\n")

# Extract SNP columns
snp_cols <- grep("rs.*_A", colnames(geno_data), value = TRUE)
cat("\nSNPs analyzed:", snp_cols, "\n\n")

# STEP 1: Calculate allele frequencies and t-tests 
cat("Sex-Stratified Allele Frequency Analysis \n\n")

for(snp in snp_cols) {
  cat(" Analysis for", snp, "\n")
  
  # Get genotype counts by sex
  male_geno <- geno_data[SEX == 1, get(snp)]
  female_geno <- geno_data[SEX == 2, get(snp)]
  
  # Calculate allele frequencies (divide by 2 because genotypes are 0/1/2)
  male_freq <- mean(male_geno, na.rm = TRUE) / 2
  female_freq <- mean(female_geno, na.rm = TRUE) / 2
  
  cat("Male allele frequency:", round(male_freq, 4), "\n")
  cat("Female allele frequency:", round(female_freq, 4), "\n")
  cat("Difference:", round(female_freq - male_freq, 4), "\n")
  
  # Statistical test for difference
  t_test <- t.test(male_geno, female_geno)
  cat("T-test p-value:", t_test$p.value, "\n\n")
}

# STEP 2: Calculate proper effect sizes using linear regression 
cat("\n Linear Regression Analysis for Effect Sizes \n\n")

effect_sizes <- c()

for(snp in snp_cols) {
  # Linear regression: genotype count ~ sex
  # SEX coding: 1=male, 2=female
  model <- glm(get(snp) ~ SEX, data = geno_data, family = gaussian())
  
  # Extract beta coefficient for SEX
  beta_coef <- coef(model)[2]
  p_value <- summary(model)$coefficients[2, 4]
  
  effect_sizes[snp] <- beta_coef
  
  cat("SNP:", snp, "\n")
  cat("  Beta coefficient:", round(beta_coef, 4), "\n")
  cat("  P-value:", p_value, "\n\n")
}

#  STEP 3: Calculate PRS using regression-derived weights 
cat("\n Polygenic Risk Score Calculation \n")

# Function to calculate PRS with proper beta weights
calculate_prs <- function(genotype_data) {
  # Extract SNP genotypes as matrix
  snp_matrix <- as.matrix(genotype_data[, .(rs2501352_A, rs1426810_A, 
                                             rs10965183_A, rs7258841_A)])
  
  # Use regression-derived beta coefficients as weights
  weights <- c(
    rs2501352_A = 0.0593,   # CRP
    rs1426810_A = -0.0540,  # ADIPOQ (negative = higher in males)
    rs10965183_A = 0.0638,  # 9p21
    rs7258841_A = 0.0947    # APOE
  )
  
  # Calculate weighted PRS for each individual
  prs_scores <- as.vector(snp_matrix %*% weights)
  
  return(prs_scores)
}

# Calculate PRS for all individuals
geno_data[, PRS := calculate_prs(.SD)]

# STEP 4: Calculate summary statistics by sex
cat("\n PRS Summary by Sex\n")

summary_stats <- geno_data[, .(
  N = .N,
  Mean_PRS = mean(PRS),
  SD_PRS = sd(PRS),
  Min_PRS = min(PRS),
  Max_PRS = max(PRS)
), by = SEX]

print(summary_stats)

# Calculate 95% confidence intervals
calculate_ci <- function(data, conf_level = 0.95) {
  n <- length(data)
  mean_val <- mean(data)
  se <- sd(data) / sqrt(n)
  alpha <- 1 - conf_level
  t_val <- qt(1 - alpha/2, df = n-1)
  
  return(list(
    mean = mean_val,
    lower_ci = mean_val - t_val * se,
    upper_ci = mean_val + t_val * se
  ))
}

male_ci <- calculate_ci(geno_data[SEX == 1, PRS])
female_ci <- calculate_ci(geno_data[SEX == 2, PRS])

cat("\n Statistical Validation \n")
cat("Males - Mean PRS:", round(male_ci$mean, 4), 
    " (95% CI:", round(male_ci$lower_ci, 4), "-", round(male_ci$upper_ci, 4), ")\n")
cat("Females - Mean PRS:", round(female_ci$mean, 4), 
    " (95% CI:", round(female_ci$lower_ci, 4), "-", round(female_ci$upper_ci, 4), ")\n")

# Sex difference test
prs_comparison <- t.test(geno_data[SEX == 1, PRS], geno_data[SEX == 2, PRS])
cat("Sex difference p-value:", format(prs_comparison$p.value, scientific = TRUE), "\n")

# Calculate percent difference
percent_diff <- ((female_ci$mean - male_ci$mean) / male_ci$mean) * 100
cat("Percent difference: Females", round(percent_diff, 1), "% higher than males\n")

# STEP 5: Risk stratification 
cat("\n Risk Stratification \n")

# Add sex-specific percentile rankings
geno_data[, PRS_percentile := ifelse(SEX == 1, 
                                      rank(PRS[SEX == 1]) / sum(SEX == 1) * 100,
                                      rank(PRS[SEX == 2]) / sum(SEX == 2) * 100)]

# Create risk categories
geno_data[, Risk_Category := ifelse(PRS_percentile >= 95, "Very High Risk",
                              ifelse(PRS_percentile >= 80, "High Risk",
                              ifelse(PRS_percentile <= 20, "Low Risk", 
                                     "Intermediate Risk")))]

# Risk summary by sex
risk_summary <- geno_data[, .N, by = .(SEX, Risk_Category)]
risk_summary_wide <- dcast(risk_summary, Risk_Category ~ SEX, value.var = "N", fill = 0)
setnames(risk_summary_wide, c("1", "2"), c("Males", "Females"))

# Add percentages
risk_summary_wide[, Male_Percent := round(Males / sum(Males) * 100, 1)]
risk_summary_wide[, Female_Percent := round(Females / sum(Females) * 100, 1)]

print(risk_summary_wide)

# Save results
fwrite(geno_data, "results/prs_scores_final.csv")

cat("\nAnalysis complete! Results saved to results/prs_scores_final.csv\n")
