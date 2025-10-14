#!/bin/bash
#SBATCH --job-name=prs_analysis
#SBATCH --account=vgopalakrishnan
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

# Script: 04_extract_target_snps.sh
# Author: Smmrithi Ravindran
# Description: Extract the 4 target SNPs for PRS calculation

module load plink/1.90b7.7
cd /ix1/vgopalakrishnan/thesis_genomics_analysis

# Create file with target SNPs
cat > target_snps.txt << EOF
rs7258841
rs10965183
rs2501352
rs1426810
EOF

# Extract the 4 target SNPs and convert to raw format for analysis
echo "Extracting target SNPs..."
plink --bfile data/GENOTYPE_DATA \
      --extract target_snps.txt \
      --recode A \
      --out target_snps_raw

echo "Target SNP extraction complete!"
echo "Output file: target_snps_raw.raw"
