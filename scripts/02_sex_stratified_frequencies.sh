#!/bin/bash
#SBATCH --job-name=sex_frequencies
#SBATCH --account=vgopalakrishnan
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

# Script: 02_sex_stratified_frequencies.sh
# Author: Smmrithi Ravindran
# Description: Calculate sex-stratified allele frequencies for each gene region

module load plink/1.90b7.7
cd /ix1/vgopalakrishnan/thesis_genomics_analysis

# Array of gene regions
GENES=("chr19_data" "crp_region" "9p21_region" "adipoq_region" "il6_region")

# Calculate frequencies for each gene region, stratified by sex
for GENE in "${GENES[@]}"; do
    echo "Processing ${GENE}..."
    
    # Male frequencies
    plink --bfile results/${GENE} \
          --keep results/males.txt \
          --freq \
          --out results/${GENE}_males_freq
    
    # Female frequencies
    plink --bfile results/${GENE} \
          --keep results/females.txt \
          --freq \
          --out results/${GENE}_females_freq
done

echo "Frequency calculation complete!"
