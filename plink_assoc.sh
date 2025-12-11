#!/bin/bash
#PLINK Association Testing Pipeline
#Association Testing Parameters-
#Test: Logistic regression (binary trait: case/control)
#Adjustment: None (primary analysis)
#Output: Odds ratios, p-values, confidence intervals

#Step 1: Case-control association test (aortic stenosis)
plink --bfile qc_hwe_filtered --logistic --ci 0.95 --out aortic_stenosis_association

#Step 2: Test specific SNP rs11166276 (novel finding)
plink --bfile qc_hwe_filtered --snp rs11166276 --logistic --ci 0.95 --out rs11166276_detailed

#Step 3: Quantitative trait analysis (serum phosphate)
#Note- Requires phosphate phenotype file
plink --bfile qc_hwe_filtered --pheno serum_phosphate.txt --linear --ci 0.95 --out serum_phosphate_association

#Step 4: Test all 15 European discovery SNPs
#Create SNP list file containing rs IDs

#Test European SNPs specifically
plink --bfile qc_hwe_filtered --extract european_snps.txt --logistic --ci 0.95 --out european_snps_replication

#Step 5: Calculate lambda (genomic inflation factor)
#Used R for calculation
#{ Lambda calculation: median(qchisq(1-p,1))/qchisq(0.5,1) }

#Step 6: Generate Manhattan plot data
plink --bfile qc_hwe_filtered --assoc --adjust --out manhattan_plot_data

#Step 7: Family-based analysis
plink --bfile qc_hwe_filtered --tdt --out family_based_test
