#!/bin/bash
# Script: 03_identify_sex_differences.sh
# Author: Smmrithi Ravindran
# Description: Identify SNPs with >2% sex difference in allele frequency

cd /ix1/vgopalakrishnan/thesis_genomics_analysis

# Array of gene regions
GENES=("chr19_data" "crp_region" "9p21_region" "adipoq_region" "il6_region")

for GENE in "${GENES[@]}"; do
    echo "Finding sex differences for ${GENE}..."
    
    # Find SNPs with >2% sex difference
    awk 'NR==FNR{male[$2]=$5; next} 
         NR>1 && $2 in male {
             diff=$5-male[$2]; 
             if(diff>0.02 || diff<-0.02) 
                 print $2, male[$2], $5, diff
         }' \
         results/${GENE}_males_freq.frq \
         results/${GENE}_females_freq.frq \
         | sort -k4 -nr \
         > results/${GENE}_gender_differences.txt
    
    # Count how many SNPs found
    COUNT=$(wc -l < results/${GENE}_gender_differences.txt)
    echo "${GENE}: ${COUNT} SNPs with >2% sex difference"
done

echo "Sex difference identification complete!"
