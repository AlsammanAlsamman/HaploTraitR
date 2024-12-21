#' Generate boxplots for genotype-phenotype data
#' @param geno_pheno A data frame containing genotype-phenotype data
#' @param outfolder An optional output folder for saving the plots
#' @return A list of ggplot2 objects or NULL if outfolder is specified
#' @export
#' @import ggplot2
#' @import ggpubr
#' @import dplyr
boxplot_geno_pheno <- function(geno_pheno, outfolder = NULL) {
  # Extract unique SNP values
  snp_var <- unique(geno_pheno$value)
  # Generate consistent colors for SNP values
  snp_colors <- rainbow(length(snp_var))
  names(snp_colors) <- snp_var  # Assign names for consistent mapping

  snps_target <- unique(geno_pheno$variable)
  plots <- list()

  # Loop through the SNP targets
  for (snp in snps_target) {
    # Subset data for the SNP
    sub_data <- geno_pheno[geno_pheno$variable %in% snp, ]
    colnames(sub_data) <- c("Trait", "snp_var", "value")
    sub_data$value <- as.factor(sub_data$value)
    # Remove the NN values
    sub_data <- sub_data[!sub_data$value == "NN", ]

    # Perform the t-test using compare_means()
    stat.test <- tryCatch({
      test_result <- compare_means(
        formula = Trait ~ value,  # Define the formula for comparison
        data = sub_data,          # Use the subset data
        method = "t.test"         # Specify Welch Two Sample t-test
      )
      test_result$y.position <- max(sub_data$Trait) + 1  # Position above the boxplot
      test_result
    }, error = function(e) {
      print(paste("Error in", snp, ":", e))
      NULL
    })

    # Create the boxplot with annotations
    p <- ggplot(sub_data, aes(x = value, y = Trait)) +
      geom_violin(aes(fill = value, alpha = 0.01)) +
      geom_jitter(width = 0.2, alpha = 0.3) +
      scale_fill_manual(values = snp_colors) +
      theme(legend.position = "none") +  # No legend
      labs(title = paste("Phenotype Distribution by SNP", snp),
           x = "Haplotype Combination",
           y = "Phenotype Value")

    if (!is.null(stat.test)) {
      p <- p + stat_pvalue_manual(stat.test[1, ], label = "p.signif", tip.length = 0.01)
    }

    # Save or store the plot
    if (!is.null(outfolder)) {
      ggsave(file.path(outfolder, paste0(snp, "_boxplot.png")), plot = p, width = 6, height = 6)
    } else {
      plots[[snp]] <- p
    }
  }

  if (is.null(outfolder)) {
    return(plots)
  }
}
