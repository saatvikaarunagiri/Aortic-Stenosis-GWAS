#Association Testing & Analysis
#Load libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(gridExtra)

setwd("/path/to/my/data")  
output_dir <- "statistical_genetics_output"

#QUESTION 5: POWER CALCULATIONS FOR "AS" REPLICATION
#Study parameters
n_cases <- 1100
n_controls <- 5900
n_total <- n_cases + n_controls
prevalence <- n_cases / n_total
alpha <- 0.00333  #Bonferroni-corrected for 15 tests

#SNP 1: rs59030006 (MECOM)
snp1 <- list(
  name = "rs59030006",
  gene = "MECOM",
  maf_sa = 0.498,
  maf_eur = 0.500,
  or_eur = 1.07,
  p_eur = 1.34e-7
)

#SNP 2: rs3901734 (TEX41)
snp2 <- list(
  name = "rs3901734",
  gene = "TEX41",
  maf_sa = 0.250,
  maf_eur = 0.250,
  or_eur = 1.13,
  p_eur = 5.92e-15
)

#Function to calculate genetic power (case-control)
calculate_power_cc <- function(n_cases, n_controls, maf, or, alpha, model = "additive") {
  # Calculate genotype relative risks
  if (model == "additive") {
    grr_het <- or
    grr_hom <- or^2
  }
  
  #Expected genotype frequencies in controls (HWE)
  p <- maf
  q <- 1 - p
  
  freq_aa_control <- q^2
  freq_ab_control <- 2*p*q
  freq_bb_control <- p^2
  
  #Expected genotype frequencies in cases
  prev <- n_cases / (n_cases + n_controls)
  k <- prev  #Population prevalence
  
  #Calculate genotype frequencies in cases using Bayes' theorem
  freq_aa_case <- freq_aa_control * 1 * k / 
    (freq_aa_control * 1 + freq_ab_control * grr_het + freq_bb_control * grr_hom)
  freq_ab_case <- freq_ab_control * grr_het * k / 
    (freq_aa_control * 1 + freq_ab_control * grr_het + freq_bb_control * grr_hom)
  freq_bb_case <- freq_bb_control * grr_hom * k / 
    (freq_aa_control * 1 + freq_ab_control * grr_het + freq_bb_control * grr_hom)
  
  #Normalize
  total <- freq_aa_case + freq_ab_case + freq_bb_case
  freq_aa_case <- freq_aa_case / total
  freq_ab_case <- freq_ab_case / total
  freq_bb_case <- freq_bb_case / total
  
  #Non-centrality parameter (Cochran-Armitage trend test)
  p_case <- freq_ab_case/2 + freq_bb_case
  p_control <- freq_ab_control/2 + freq_bb_control
  
  ncp <- n_cases * n_controls / (n_cases + n_controls) * 
    (p_case - p_control)^2 / (p * q)
  
  #Critical value
  critical_chi <- qchisq(1 - alpha, df = 1)
  
  #Power
  power <- 1 - pchisq(critical_chi, df = 1, ncp = ncp)
  
  #Calculate sample size for 80% power
  n_80_ncp <- (qnorm(0.80) + qnorm(1 - alpha/2))^2
  n_80_cases <- n_80_ncp * (n_cases + n_controls) / ncp * n_cases / (n_cases + n_controls)
  
  return(list(
    ncp = ncp,
    power = power,
    n_cases_80 = ceiling(n_80_cases)
  ))
}

#Calculate power for both SNPs
power_snp1 <- calculate_power_cc(n_cases, n_controls, snp1$maf_sa, snp1$or_eur, alpha)
power_snp2 <- calculate_power_cc(n_cases, n_controls, snp2$maf_sa, snp2$or_eur, alpha)

#Create results table
power_as_results <- data.frame(
  SNP = c(snp1$name, snp2$name),
  Gene = c(snp1$gene, snp2$gene),
  MAF_SA = c(snp1$maf_sa, snp2$maf_sa),
  OR_EUR = c(snp1$or_eur, snp2$or_eur),
  P_EUR = c(snp1$p_eur, snp2$p_eur),
  Sample_Size = paste0(n_cases, " cases / ", n_controls, " controls"),
  Power_Percent = paste0(round(c(power_snp1$power, power_snp2$power) * 100, 2), "%"),
  N_Cases_80pct = c(power_snp1$n_cases_80, power_snp2$n_cases_80),
  Fold_Increase = round(c(power_snp1$n_cases_80 / n_cases, 
                           power_snp2$n_cases_80 / n_cases), 1)
)

print(power_as_results)

#Save results
write.csv(power_as_results, 
          file.path(output_dir, "Q5_power_calculations_AS.csv"), 
          row.names = FALSE)

#Visualization (Power curve)
power_range <- seq(0.01, 0.50, by = 0.01)

power_curves <- lapply(power_range, function(maf) {
  p1 <- calculate_power_cc(n_cases, n_controls, maf, 1.07, alpha)$power
  p2 <- calculate_power_cc(n_cases, n_controls, maf, 1.13, alpha)$power
  data.frame(MAF = maf, Power_OR1.07 = p1, Power_OR1.13 = p2)
})

power_curves_df <- do.call(rbind, power_curves)

png(file.path(output_dir, "Q5_power_curves_AS.png"), 
    width = 12, height = 8, units = "in", res = 300)
ggplot(power_curves_df) +
  geom_line(aes(x = MAF, y = Power_OR1.07, color = "OR = 1.07 (rs59030006)"), 
            size = 1.5) +
  geom_line(aes(x = MAF, y = Power_OR1.13, color = "OR = 1.13 (rs3901734)"), 
            size = 1.5) +
  geom_point(data = data.frame(MAF = snp1$maf_sa, Power = power_snp1$power),
             aes(x = MAF, y = Power), color = "#e74c3c", size = 4) +
  geom_point(data = data.frame(MAF = snp2$maf_sa, Power = power_snp2$power),
             aes(x = MAF, y = Power), color = "#3498db", size = 4) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("#e74c3c", "#3498db")) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Statistical Power for AS Replication",
    subtitle = paste0("N = ", n_cases, " cases, ", n_controls, " controls; α = ", alpha),
    x = "Minor Allele Frequency",
    y = "Statistical Power",
    color = "SNP"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 16)
  )
dev.off()

#QUESTION 6: POWER CALCULATIONS FOR SPP ASSOCIATION
#Study parameters
n_spp <- 7000  #All individuals have SPP measurements
alpha_spp <- 0.00333  #Bonferroni-corrected for 15 tests

#Function to calculate power for quantitative trait
calculate_power_quant <- function(n, maf, beta, alpha) {
  #Variance explained
  var_explained <- 2 * maf * (1 - maf) * beta^2
  
  #Non-centrality parameter
  ncp <- n * var_explained / (1 - var_explained)
  
  #Critical value
  critical_chi <- qchisq(1 - alpha, df = 1)
  
  #Power
  power <- 1 - pchisq(critical_chi, df = 1, ncp = ncp)
  
  #Sample size for 80% power
  n_80 <- ceiling((qnorm(0.80) + qnorm(1 - alpha/2))^2 / var_explained)
  
  return(list(
    var_explained = var_explained,
    ncp = ncp,
    power = power,
    n_80 = n_80
  ))
}

#Test effect sizes
effect_sizes <- c(0.03, 0.05, 0.08)  #SD units
effect_names <- c("Small", "Moderate", "Large")

#Calculate power for both SNPs across effect sizes
power_spp_list <- list()

for (i in seq_along(effect_sizes)) {
  beta <- effect_sizes[i]
  
  #SNP 1
  p1 <- calculate_power_quant(n_spp, snp1$maf_sa, beta, alpha_spp)
  
  #SNP 2
  p2 <- calculate_power_quant(n_spp, snp2$maf_sa, beta, alpha_spp)
  
  power_spp_list[[i]] <- data.frame(
    SNP = c(snp1$name, snp2$name),
    Gene = c(snp1$gene, snp2$gene),
    MAF = c(snp1$maf_sa, snp2$maf_sa),
    Effect_Size = effect_names[i],
    Beta_SD = beta,
    Var_Explained = paste0(round(c(p1$var_explained, p2$var_explained) * 100, 2), "%"),
    Power = paste0(round(c(p1$power, p2$power) * 100, 2), "%"),
    N_for_80pct = c(p1$n_80, p2$n_80)
  )
}

power_spp_results <- do.call(rbind, power_spp_list)
print(power_spp_results)

#Save results
write.csv(power_spp_results, 
          file.path(output_dir, "Q6_power_calculations_SPP.csv"), 
          row.names = FALSE)

#Visualization: Power curves for different effect sizes
beta_range <- seq(0.02, 0.10, by = 0.005)

power_spp_curves <- lapply(beta_range, function(beta) {
  p1 <- calculate_power_quant(n_spp, snp1$maf_sa, beta, alpha_spp)$power
  p2 <- calculate_power_quant(n_spp, snp2$maf_sa, beta, alpha_spp)$power
  data.frame(Beta = beta, 
             Power_rs59030006 = p1, 
             Power_rs3901734 = p2)
})

power_spp_curves_df <- do.call(rbind, power_spp_curves)

png(file.path(output_dir, "Q6_power_curves_SPP.png"), 
    width = 12, height = 8, units = "in", res = 300)
ggplot(power_spp_curves_df) +
  geom_line(aes(x = Beta, y = Power_rs59030006, 
                color = "rs59030006 (MAF=0.50)"), size = 1.5) +
  geom_line(aes(x = Beta, y = Power_rs3901734, 
                color = "rs3901734 (MAF=0.25)"), size = 1.5) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = c(0.03, 0.05, 0.08), linetype = "dotted", 
             color = "gray70", alpha = 0.5) +
  scale_color_manual(values = c("#e74c3c", "#3498db")) +
  scale_y_continuous(labels = scales::percent) +
  annotate("text", x = 0.03, y = 0.05, label = "Small", angle = 90, size = 3) +
  annotate("text", x = 0.05, y = 0.05, label = "Moderate", angle = 90, size = 3) +
  annotate("text", x = 0.08, y = 0.05, label = "Large", angle = 90, size = 3) +
  labs(
    title = "Statistical Power for SPP Associations",
    subtitle = paste0("N = ", n_spp, "; α = ", alpha_spp),
    x = "Effect Size (β, SD units)",
    y = "Statistical Power",
    color = "SNP"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 16)
  )
dev.off()

# QUESTIONS 8-10: ASSOCIATION TESTING
#Process PLINK logistic regression results

process_plink_logistic <- function(results_file, 
                                    sig_threshold = 0.00333,
                                    output_prefix = "Q8_AS_associations") {
  
  #Read PLINK results
  results <- fread(results_file)
  
  #Filter for ADD (additive) model
  results <- results[TEST == "ADD"]
  
  #Add significance flags
  results$Significant <- results$P < sig_threshold
  results$Suggestive <- results$P >= sig_threshold & results$P < 0.05
  
  #Sort by p-value
  results <- results[order(P)]
  
  #Create summary
  summary_stats <- data.frame(
    Total_SNPs = nrow(results),
    Significant = sum(results$Significant, na.rm = TRUE),
    Suggestive = sum(results$Suggestive, na.rm = TRUE),
    Non_significant = sum(results$P >= 0.05, na.rm = TRUE)
  )

  print(summary_stats)
  print(results[1:5, .(SNP, CHR, OR, P, Significant, Suggestive)])

  #Manhattan plot
  results$log10p <- -log10(results$P)
  
  png(paste0(output_dir, "/", output_prefix, "_manhattan.png"), 
      width = 14, height = 8, units = "in", res = 300)
  
  p <- ggplot(results, aes(x = seq_along(SNP), y = log10p)) +
    geom_point(aes(color = as.factor(CHR)), size = 3, alpha = 0.7) +
    geom_hline(yintercept = -log10(sig_threshold), linetype = "dashed", 
               color = "#e74c3c", size = 1) +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", 
               color = "gray50", size = 0.8) +
    scale_x_continuous(breaks = seq_along(results$SNP), 
                       labels = results$SNP) +
    labs(
      title = "SNP-AS Association Results",
      subtitle = paste0("Red line: Bonferroni threshold (α = ", sig_threshold, ")"),
      x = "SNP",
      y = expression(-log[10](P))
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 16)
    )
  
  print(p)
  dev.off()
  
  #Save results
  write.csv(results, 
            paste0(output_dir, "/", output_prefix, "_full_results.csv"), 
            row.names = FALSE)
  
  return(list(results = results, summary = summary_stats))
}

#Process PLINK linear regression results (SPP)
process_plink_linear <- function(results_file, 
                                   sig_threshold = 0.00333,
                                   output_prefix = "Q9_SPP_associations") {
  
  #Read PLINK results
  results <- fread(results_file)
  
  #Filter for ADD (additive) model
  results <- results[TEST == "ADD"]
  
  #Add significance flags
  results$Significant <- results$P < sig_threshold
  results$Suggestive <- results$P >= sig_threshold & results$P < 0.05
  
  #Sort by p-value
  results <- results[order(P)]
  
  #Create summary
  summary_stats <- data.frame(
    Total_SNPs = nrow(results),
    Significant = sum(results$Significant, na.rm = TRUE),
    Suggestive = sum(results$Suggestive, na.rm = TRUE),
    Non_significant = sum(results$P >= 0.05, na.rm = TRUE)
  )

  print(summary_stats)
  print(results[1:5, .(SNP, CHR, BETA, P, Significant, Suggestive)])
 
  #Effect size plot
  results$log10p <- -log10(results$P)
  
  png(paste0(output_dir, "/", output_prefix, "_effect_sizes.png"), 
      width = 12, height = 8, units = "in", res = 300)
  
  p <- ggplot(results, aes(x = reorder(SNP, -log10p), y = BETA)) +
    geom_bar(stat = "identity", aes(fill = Significant), alpha = 0.8) +
    geom_errorbar(aes(ymin = BETA - 1.96*SE, ymax = BETA + 1.96*SE), 
                  width = 0.3) +
    scale_fill_manual(values = c("FALSE" = "gray70", "TRUE" = "#e74c3c")) +
    labs(
      title = "SNP Effect Sizes on Serum Phosphate",
      subtitle = "Error bars: 95% confidence intervals",
      x = "SNP",
      y = "Effect Size (β, mg/dL per allele)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold", size = 16)
    ) +
    coord_flip()
  
  print(p)
  dev.off()
  
  #Save results
  write.csv(results, 
            paste0(output_dir, "/", output_prefix, "_full_results.csv"), 
            row.names = FALSE)
  
  return(list(results = results, summary = summary_stats))
}
