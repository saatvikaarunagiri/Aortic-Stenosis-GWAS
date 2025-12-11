#Statistical Genetics Analysis Pipeline
#Project: Aortic Stenosis Genetic Association Study in South Asians
#Course: BS858 Statistical Genetics I

#Description: Complete analysis pipeline for genetic association study of aortic stenosis
#and serum phosphate in South Asian population. 
#Replication of Trenkwalder et al. European GWAS findings with trans-ancestry comparison.

#Load required libraries
required_packages <- c(
  "data.table",
  "dplyr",
  "ggplot2",
  "epitools",
  "pwr",
  "gridExtra",
  "scales",
  "RColorBrewer"
)

# Install missing packages
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(epitools)
library(pwr)
library(gridExtra)
library(scales)
library(RColorBrewer)

#Set working directory
setwd("/path/to/my/data")

#Set random seed for reproducibility
set.seed(42)

#Define output directory
output_dir <- "statistical_genetics_output"
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

#QUESTION 1: HERITABILITY & FAMILIALITY ANALYSIS
#Load data
data <- fread("project_sib_pheno_and_RV_data.csv")

cat("\nData loaded successfully!\n")
cat("Sample size:", nrow(data), "sibling pairs\n")
cat("Variables:", ncol(data), "\n\n")

#1a. HERITABILITY OF SERUM PHOSPHATE (Continuous Trait)
#Calculate sibling correlation for STSPP
cor_spp <- cor(data$STSPP1, data$STSPP2, use = "complete.obs")
cat("Sibling correlation for STSPP:", round(cor_spp, 3), "\n")

#Estimate heritability
h2_spp <- 2 * cor_spp
cat("Estimated heritability (h²) for STSPP:", round(h2_spp, 3), 
    paste0("(", round(h2_spp*100, 1), "%)"), "\n\n")

#Statistical test
cor_test_spp <- cor.test(data$STSPP1, data$STSPP2)
print(cor_test_spp)

#Create results dataframe
spp_heritability_results <- data.frame(
  Statistic = c("Sibling correlation (r)", "Estimated heritability (h²)", 
                "Test statistic", "P-value", "95% CI (lower)", "95% CI (upper)"),
  Value = c(
    round(cor_spp, 3),
    paste0(round(h2_spp, 3), " (", round(h2_spp*100, 1), "%)"),
    paste0("t = ", round(cor_test_spp$statistic, 3), 
           ", df = ", cor_test_spp$parameter),
    format(cor_test_spp$p.value, scientific = TRUE, digits = 3),
    round(cor_test_spp$conf.int[1], 3),
    round(cor_test_spp$conf.int[2], 3)
  )
)

print(spp_heritability_results)

#Save results
write.csv(spp_heritability_results, 
          file.path(output_dir, "Q1a_spp_heritability.csv"), 
          row.names = FALSE)

#Visualization: Scatter plot of sibling correlations
png(file.path(output_dir, "Q1a_spp_sibling_correlation.png"), 
    width = 10, height = 8, units = "in", res = 300)
ggplot(data, aes(x = STSPP1, y = STSPP2)) +
  geom_point(alpha = 0.3, color = "#3498db") +
  geom_smooth(method = "lm", color = "#e74c3c", size = 1.5) +
  labs(
    title = "Serum Phosphate Sibling Correlation",
    subtitle = paste0("r = ", round(cor_spp, 3), 
                      ", h² = ", round(h2_spp, 3),
                      ", p < 2.2e-16"),
    x = "Sibling 1 Serum Phosphate (mg/dL)",
    y = "Sibling 2 Serum Phosphate (mg/dL)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12, color = "gray30")
  )
dev.off()

cat("\n SPP heritability analysis complete!\n")
cat("  h2 =", round(h2_spp, 3), "(65.6% of variance explained by genetics)\n\n")

#1b. FAMILIALITY OF AORTIC STENOSIS (Binary Trait)
cat("1b. AORTIC STENOSIS FAMILIALITY\n")

#Create contingency table for AS status
table_AS <- table(data$AS1, data$AS2)
rownames(table_AS) <- c("Proband Unaffected", "Proband Affected")
colnames(table_AS) <- c("Sibling Unaffected", "Sibling Affected")

cat("AS Status Contingency Table:\n")
print(table_AS)
cat("\n")

#Calculate prevalence
prev_proband <- sum(data$AS1 == 2) / nrow(data)
prev_sibling <- sum(data$AS2 == 2) / nrow(data)

#Calculate risk to sibling given proband affected
risk_sib_if_proband_affected <- sum(data$AS1 == 2 & data$AS2 == 2) / sum(data$AS1 == 2)

#Recurrence risk ratio
lambda_s <- risk_sib_if_proband_affected / prev_sibling

cat("Prevalence in probands:", round(prev_proband, 3), "\n")
cat("Prevalence in siblings:", round(prev_sibling, 3), "\n")
cat("Risk to sibling if proband affected:", round(risk_sib_if_proband_affected, 3), "\n")
cat("Recurrence risk ratio (λs):", round(lambda_s, 2), "\n\n")

#Chi-square test for independence
chisq_test_AS <- chisq.test(table_AS)
print(chisq_test_AS)
cat("\n")

#Calculate odds ratio
or_result <- oddsratio(table_AS)
print(or_result)

#Create results dataframe
as_familiality_results <- data.frame(
  Statistic = c("Total sibling pairs", "Proband prevalence", "Sibling prevalence",
                "Risk to sibling (proband affected)", "Recurrence risk ratio (λs)",
                "Chi-square statistic", "Chi-square p-value",
                "Odds ratio", "OR 95% CI lower", "OR 95% CI upper"),
  Value = c(
    nrow(data),
    round(prev_proband, 3),
    round(prev_sibling, 3),
    round(risk_sib_if_proband_affected, 3),
    round(lambda_s, 2),
    round(chisq_test_AS$statistic, 3),
    round(chisq_test_AS$p.value, 3),
    round(or_result$measure[2, "estimate"], 3),
    round(or_result$measure[2, "lower"], 3),
    round(or_result$measure[2, "upper"], 3)
  )
)

print(as_familiality_results)

#Save results
write.csv(as_familiality_results, 
          file.path(output_dir, "Q1b_as_familiality.csv"), 
          row.names = FALSE)

#Visualization: Mosaic plot
png(file.path(output_dir, "Q1b_as_contingency_mosaic.png"), 
    width = 10, height = 8, units = "in", res = 300)
mosaicplot(table_AS, 
           main = "Aortic Stenosis Concordance in Sibling Pairs",
           xlab = "Proband AS Status",
           ylab = "Sibling AS Status",
           color = c("#2ecc71", "#e74c3c"),
           las = 1)
dev.off()

cat("\n AS familiality analysis complete!\n")
cat("  lambda_s =", round(lambda_s, 2), "(modest familial aggregation)\n")
cat("  Not significant (p = 0.44)\n\n")

#HELPER FUNCTIONS FOR SUBSEQUENT ANALYSES
#Function to calculate power for quantitative traits
calculate_power_quantitative <- function(n, maf, beta, alpha = 0.05) {
  #Variance explained by SNP
  var_explained <- 2 * maf * (1 - maf) * beta^2
  
  #Non-centrality parameter
  ncp <- n * var_explained / (1 - var_explained)
  
  #Power calculation
  critical_value <- qchisq(1 - alpha, df = 1)
  power <- 1 - pchisq(critical_value, df = 1, ncp = ncp)
  
  return(list(
    var_explained = var_explained,
    ncp = ncp,
    power = power
  ))
}

#Function to format p-values
format_pval <- function(p) {
  if (p < 0.001) {
    return(format(p, scientific = TRUE, digits = 2))
  } else {
    return(round(p, 3))
  }
}

#Function to create Manhattan-style plot
manhattan_plot <- function(results_df, sig_threshold, title) {
  results_df$log10p <- -log10(results_df$P)
  results_df$chr_factor <- factor(results_df$CHR)
  
  #Color by chromosome
  chr_colors <- rep(c("#3498db", "#95a5a6"), length.out = length(unique(results_df$CHR)))
  
  p <- ggplot(results_df, aes(x = seq_along(SNP), y = log10p, color = chr_factor)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_hline(yintercept = -log10(sig_threshold), linetype = "dashed", 
               color = "#e74c3c", size = 1) +
    scale_color_manual(values = chr_colors) +
    labs(
      title = title,
      subtitle = paste0("Significance threshold: α = ", sig_threshold),
      x = "SNP",
      y = expression(-log[10](P))
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 16),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    scale_x_continuous(breaks = seq_along(results_df$SNP), 
                       labels = results_df$SNP)
  
  return(p)
}

cat("\n Helper functions loaded successfully!\n\n")
