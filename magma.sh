#!/bin/bash
#unique ID (chr:pos:allele1:allele2) is used
#GWAS summary data is in format of SNP CHR BP A1 A2 FRQ BETA SE P N

# ==================== CONFIGURATION ====================
# Set base directories 
BASE_DIR="/path/to/the/project"
DATA_DIR="${BASE_DIR}/data"
REF_DIR="${BASE_DIR}/reference"
OUTPUT_DIR="${BASE_DIR}/results"

# ==================== MAIN PIPELINE ====================
# Gene Annotation
magma --annotate window=35,10 \
      --snp-loc "${DATA_DIR}/GWAS_mdd2025_eur.txt" \
      --gene-loc "${REF_DIR}/NCBI37.3.gene.loc.MHCexcluded" \
      --out "${OUTPUT_DIR}/Gannotation"


# Gene Analysis
magma --bfile "${REF_DIR}/g1000_eur" \
      --pval "${DATA_DIR}/GWAS_mdd2025_eur.txt" ncol=N \
      --gene-annot "${OUTPUT_DIR}/Gannotation.genes.annot" \
      --out "${OUTPUT_DIR}/Ganalysis"


# Gene property analysis for tissue specificity 
magma --gene-results "${OUTPUT_DIR}/Ganalysis.genes.raw" \
      --gene-covar "${REF_DIR}/gtex_v8_26_avg_log2TPM.txt" \
      --model direction-covar=positive condition-hide=Average \
      --out "${OUTPUT_DIR}/Gproperty"
    

#Gene binary analysis for tissue specificity
magma  --gene-results "${OUTPUT_DIR}/Ganalysis.genes.raw" \
          --set-annot "${REF_DIR}/gtexset_genes.txt" col=2,1 \
          --out "${OUTPUT_DIR}/Gbinary"