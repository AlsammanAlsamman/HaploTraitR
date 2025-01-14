#' Convert t-test results to a data frame and save to a file
#'
#' @param t_test_snpComp A list containing the t-test results for SNP comparisons
#' @param outfolder The output folder for saving the file
#' @export
saveTTestResultsToFile <- function(t_test_snpComp) {
  # Get the phenotype name from the configuration
  phenotype_name <- get_config("phenotypename")
  outfolder <- get_config("outfolder")
  # check if the outfolder is not NULL and does exist
  if (!is.null(outfolder) && !dir.exists(outfolder)) {
    stop("The output folder does not exist, please create it first.")
  }
  # Convert t-test results to a data frame using merge_with_list_name
  t_test_df <- merge_with_list_name(t_test_snpComp, phenotype_name)

  # Save the data frame to a CSV file
  file_path <- file.path(outfolder, paste0(phenotype_name, "_comb_t_test_results.csv"))
  write.csv(t_test_df, file_path, row.names = FALSE)

  message("T-test results saved to ", file_path)
}


#' Merge a list of data frames with a list name
#'
#' @param df_list A list of data frames to merge
#' @param pheno The phenotype name to add to the merged data frame
#' @return A merged data frame with an added phenotype column
#' @importFrom dplyr mutate select
#' @importFrom purrr imap_dfr
#' @export
merge_with_list_name <- function(df_list, pheno) {
  df_list %>%
    purrr::imap_dfr(~ dplyr::mutate(.x, snp = .y)) %>%
    dplyr::select(snp, group1, group2, n1, n2, statistic, df, p, p.adj, p.adj.signif) %>%
    dplyr::mutate(pheno = pheno)
}

#' Evaluate haplotype combinations versus group 0 in each trait
#'
#' @param t_test_snpComp A list of tibbles containing t-test results for different SNP combinations
#' @return A data frame with additional columns 'eval0' and 'dir0' indicating the evaluation and direction of group 0 comparison
#' @importFrom dplyr group_by mutate case_when ungroup bind_rows if_else
#' @importFrom purrr map_dfr
#' @export
evaluate_test_group0 <- function(t_test_snpComp) {
  # Function to evaluate a single tibble
  evaluate_single <- function(data) {
    if (!"p.adj.signif" %in% colnames(data)) {
      data <- data %>% dplyr::mutate(p.adj.signif = dplyr::if_else(p < 0.05, "*", "ns"))
    }

    data %>%
      dplyr::group_by(snp) %>%
      dplyr::mutate(
        eval0 = any(p.adj.signif[group1 == "0"] %in% c("****", "***", "**", "*")),
        dir0 = dplyr::case_when(
          group1 == "0" & p.adj.signif %in% c("****", "***", "**", "*") & statistic > 0 ~ 1,
          group1 == "0" & p.adj.signif %in% c("****", "***", "**", "*") & statistic < 0 ~ -1,
          TRUE ~ 0
        )
      ) %>%
      dplyr::ungroup()
  }

  # Apply the evaluation to each tibble in the list and combine the results
  results <- purrr::map_dfr(t_test_snpComp, evaluate_single)

  return(results)
}
