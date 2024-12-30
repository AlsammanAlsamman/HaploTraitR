#' Generate boxplots for genotype-phenotype data
#'
#' @param pheno A data frame containing phenotype data
#' @param gwas A data frame containing GWAS data
#' @param hapmap A data frame containing hapmap data
#' @return A list of ggplot2 objects or NULL if outfolder is specified
#' @export
#' @import ggplot2
#' @import ggpubr
#' @import dplyr
boxplot_genotype_phenotype <- function(pheno, gwas, hapmap, plot_width = 5, plot_height = 5) {

  # Get the output folder from the config
  outfolder <- get_config("outfolder")
  # Get the phenotype name from the config
  pheno_name <- get_config("phenotypename")
  # Ensure the output directory for boxplots exists
  boxplot_dir <- file.path(outfolder, paste0("boxplots_", pheno_name))
  if (!dir.exists(boxplot_dir)) {
    dir.create(boxplot_dir, showWarnings = FALSE)
  }

  # Subset hapmap data
  subhapmap <- extract_hapmap(hapmap, gwas)
  # Create a data frame with phenotype and genotype data
  genotype_phenotype_data <- get_pheno_geno(subhapmap, pheno)

  # Access configuration parameters from the environment
  method <- get_config("t_test_method")
  y_label <- paste(get_config("phenotypename"))
  # Append unit to y_label if provided in the config
  if (!is.null(get_config("phenotypeunit"))) {
    y_label <- paste(y_label, "(", get_config("phenotypeunit"), ")")
  }

  # Ensure 'value' is treated as a factor
  genotype_phenotype_data$value <- as.factor(genotype_phenotype_data$value)

  # Generate consistent colors for SNP values
  snp_var <- unique(genotype_phenotype_data$value)
  snp_colors <- setNames(rainbow(length(snp_var)), snp_var)

  # Extract unique SNP targets
  snps_target <- unique(genotype_phenotype_data$variable)
  plots <- list()

  # Loop through the SNP targets
  for (snp in snps_target) {
    # Subset data for the SNP and remove 'NN' values
    sub_data <- subset(genotype_phenotype_data, variable == snp & value != "NN")

    # Check if sub_data is not empty
    if (nrow(sub_data) == 0) {
      message(paste("No data available for SNP:", snp))
      next
    }

    # Rename columns for clarity
    colnames(sub_data) <- c("Trait", "snp_var", "value")

    # Perform the t-test or Wilcoxon test using compare_means()
    stat.test <- tryCatch({
      test_result <- compare_means(
        formula = Trait ~ value,  # Define the formula for comparison
        data = sub_data,          # Use the subset data
        method = method           # Use the specified method
      )
      test_result$y.position <- max(sub_data$Trait) + 1  # Position above the boxplot
      test_result
    }, error = function(e) {
      message(paste("Error in SNP", snp, ":", e$message))
      NULL
    })

    # Create the boxplot with annotations
    p <- ggplot(sub_data, aes(x = value, y = Trait)) +
      geom_violin(aes(fill = value), alpha = 0.3) +
      geom_jitter(width = 0.2, alpha = 0.5) +
      scale_fill_manual(values = snp_colors) +
      theme(legend.position = "none") +
      labs(title = paste("Phenotype Distribution by SNP", snp),
           x = "Genotype",
           y = y_label)

    if (!is.null(stat.test)) {
      p <- p + stat_pvalue_manual(stat.test, label = "p.signif", tip.length = 0.01)
    }

    # Save or store the plot
    if (!is.null(outfolder)) {
      ggsave(filename = file.path(boxplot_dir, paste0(snp, "_boxplot.png")), plot = p, width = plot_width, height = plot_height)
    }
    plots[[snp]] <- p
  }

  return(plots)
}
