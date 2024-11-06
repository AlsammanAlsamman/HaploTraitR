#' Plotting Haplotype Combination Boxplot for specific SNP
#' @param snp The SNP of interest
#' @param SNPcombTables A list of data frames containing haplotype combinations and phenotype values
#' @param t_test_snpComp A list of t-test results for haplotype combinations
#' @import ggplot2
#' @import ggpubr
#' @export
plotHapCombBoxPlot <- function(snp, SNPcombTables, t_test_snpComp) {
  # Extract the data frame for the specified SNP
  comb_sample <- SNPcombTables[[snp]]

  # Create a boxplot with significance annotations
  p <- ggplot(comb_sample, aes(x = Comb, y = Pheno)) +
    geom_violin(aes(fill = Comb), alpha = 0.5) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste("Phenotype Distribution by Haplotype Combination for SNP", snp),
         x = "Haplotype Combination",
         y = "Phenotype Value") +
    theme(legend.position = "none")

  # Add significance brackets to the plot
  # if the snp has a t test result
  if (snp %in% names(t_test_snpComp)) {
    p <- p + stat_pvalue_manual(t_test_snpComp[[snp]], label = "p.adj.signif", tip.length = 0.01)
  }
  # Display the plot
  p
}


