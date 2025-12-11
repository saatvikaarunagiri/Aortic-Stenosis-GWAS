#!/bin/bash
#PLINK Quality Control Pipeline

#Quality Control Parameters
#Call rate: 95%
#Hardy-Weinberg equilibrium: p > 0.001
#Minor allele frequency: No filter (report all)

#Step 1: Calculate call rates and generate summary statistics
plink --bfile input_genotypes --missing --out qc_missing

#Step 2: Filter by call rate (95%)
plink --bfile input_genotypes --geno 0.05 --mind 0.05 --make-bed --out qc_callrate

#Step 3: Hardy-Weinberg Equilibrium test
plink --bfile qc_callrate --hardy --out qc_hwe

#Step 4: Remove SNPs failing HWE (p < 0.001)
plink --bfile qc_callrate --hwe 0.001 --make-bed --out qc_hwe_filtered

#Step 5: Calculate allele frequencies
plink --bfile qc_hwe_filtered --freq --out qc_frequencies

#Step 6: Generate final QC report
plink --bfile qc_hwe_filtered --missing --hardy --freq --out final_qc_report
