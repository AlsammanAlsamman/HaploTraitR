#' Plotting Haplotype Combination Distribution for specific SNP
#' @param snp The SNP of interest
#' @param SNPcombTables A list of data frames containing haplotype combinations and phenotype values
#' @importFrom ggplot2 ggplot geom_density theme_minimal labs
#' @export
plotHapCombDistribution <- function(snp, SNPcombTables) {
  # Extract the data frame for the specified SNP
  comb_sample <- SNPcombTables[[snp]]

  # Create a boxplot with significance annotations
  p <- ggplot(comb_sample, aes(x = Pheno, fill = Comb)) +
    geom_density(alpha = 0.5) +
    theme_minimal() +
    labs(title = paste("Phenotype Distribution by Haplotype Combination for SNP", snp),
         x = "Phenotype Value",
         y = "Density")

  # Display the plot
  p
}
