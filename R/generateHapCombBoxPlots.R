#' Generate and save haplotype combination boxplots
#'
#' @param SNPcombTables A list containing the SNP combination tables
#' @param t_test_snpComp A list containing the t-test results for SNP comparisons
#' @param pwidth The width of the saved plots (default is 6)
#' @param pheight The height of the saved plots (default is 6)
#' @import ggplot2
#' @export
generateHapCombBoxPlots <- function(SNPcombTables, t_test_snpComp, pwidth = 6, pheight = 6) {
  # Get the output folder from the configuration
  outfolder <- get_config("outfolder")
  # check if the outfolder is not NULL and does exist
  if (!is.null(outfolder) && !dir.exists(outfolder)) {
    stop("The output folder does not exist, please create it first.")
  }
  # Evaluate the SNPs to three categories: significant and positive, significant and negative, not significant
  evaluated_snps <- evaluate_test_group0(t_test_snpComp)
  snps_a <- evaluated_snps %>% dplyr::filter(eval0 == TRUE & dir0 == 1) %>% dplyr::pull(snp)
  snps_b <- evaluated_snps %>% dplyr::filter(eval0 == TRUE & dir0 == -1) %>% dplyr::pull(snp)
  snps_c <- evaluated_snps %>% dplyr::filter(eval0 == FALSE) %>% dplyr::pull(snp)

  # Create directories for the boxplots
  dir_paths <- create_dirs(outfolder)

  # Define a helper function to save plots in the appropriate directory
  save_plot <- function(cls_snp, category) {
    outdir <- dir_paths[[category]]
    HaploTraitR::plotHapCombBoxPlot(cls_snp, SNPcombTables, t_test_snpComp, outfolder = outdir)
    #ggplot2::ggsave(file.path(outdir, paste0(cls_snp, "_LD_boxplot.png")), width = pwidth, height = pheight)
  }

  # Loop over the SNP comparisons and save the plots in the corresponding directories
  lapply(snps_a, save_plot, category = "a")
  lapply(snps_b, save_plot, category = "b")
  lapply(snps_c, save_plot, category = "c")
}

#' Create directories for haplotype combination boxplots
#'
#' @param outfolder The output folder where the directories will be created
#' @return A list containing the paths of the created directories for groups a, b, and c
#' @export
create_dirs <- function(outfolder) {
  # Define the main boxplot directory
  boxplot_dir <- file.path(outfolder, "haplotype_combination_boxplots")

  # Define the subdirectory names and their corresponding labels
  subdirs <- list(
    a = "significant_positive",
    b = "significant_negative",
    c = "not_significant"
  )

  # Initialize a list to store the paths of the created directories
  dir_paths <- list()

  # Create the main boxplot directory if it doesn't exist
  if (!dir.exists(boxplot_dir)) {
    dir.create(boxplot_dir, showWarnings = FALSE)
  }

  # Loop through the subdirectories and create them if they don't exist
  for (key in names(subdirs)) {
    subdir_path <- file.path(boxplot_dir, subdirs[[key]])
    if (!dir.exists(subdir_path)) {
      dir.create(subdir_path, showWarnings = FALSE)
    }
    dir_paths[[key]] <- subdir_path
  }

  return(dir_paths)
}
