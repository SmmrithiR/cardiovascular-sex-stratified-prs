#!/bin/bash
#SBATCH --job-name=cv_genetics_extraction
#SBATCH --account=vgopalakrishnan
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

# Script: 01_data_extraction_and_setup.sh
# Author: Smmrithi Ravindran
# Date: 2025
# Description: Extract gene regions and create sex-stratified sample lists

# Load required modules
module load plink/1.90b7.7

# Set working directory
cd /ix1/vgopalakrishnan/thesis_genomics_analysis

# Create results directory
mkdir -p results

# Create male and female sample lists from family file
echo "Creating sex-specific sample lists..."
awk '$5==1 {print $1, $2}' data/GENOTYPE_DATA.fam > results/males.txt
awk '$5==2 {print $1, $2}' data/GENOTYPE_DATA.fam > results/females.txt

# Verify sample counts
echo "Males: $(wc -l < results/males.txt)"
echo "Females: $(wc -l < results/females.txt)"

# Extract chromosome 19 (APOE region)
echo "Extracting chromosome 19 (APOE)..."
plink --bfile data/GENOTYPE_DATA \
      --chr 19 \
      --make-bed \
      --out results/chr19_data

# Extract CRP region (chromosome 1, bp 159000000-160000000)
echo "Extracting CRP region..."
plink --bfile data/GENOTYPE_DATA \
      --chr 1 \
      --from-bp 159000000 \
      --to-bp 160000000 \
      --make-bed \
      --out results/crp_region

# Extract 9p21 region (chromosome 9, bp 21000000-22500000)
echo "Extracting 9p21 region..."
plink --bfile data/GENOTYPE_DATA \
      --chr 9 \
      --from-bp 21000000 \
      --to-bp 22500000 \
      --make-bed \
      --out results/9p21_region

# Extract ADIPOQ region (chromosome 3, bp 186500000-186600000)
echo "Extracting ADIPOQ region..."
plink --bfile data/GENOTYPE_DATA \
      --chr 3 \
      --from-bp 186500000 \
      --to-bp 186600000 \
      --make-bed \
      --out results/adipoq_region

# Extract IL-6 region (chromosome 7, bp 22700000-22800000)
echo "Extracting IL-6 region..."
plink --bfile data/GENOTYPE_DATA \
      --chr 7 \
      --from-bp 22700000 \
      --to-bp 22800000 \
      --make-bed \
      --out results/il6_region

echo "Data extraction complete!"
