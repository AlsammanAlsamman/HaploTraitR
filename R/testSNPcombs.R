#' Calculate the p-value for a SNP combinations
#' @param comb_sample_Tables A list of sample tables for each SNP combination
#' @param savecopy A boolean indicating whether to save a copy of the t-test results to the project folder
#' @return A data frame with the p-values for the SNP combinations
#' @importFrom rstatix t_test
#' @importFrom rstatix add_xy_position
#' @importFrom dplyr mutate
#' @export
testSNPcombs <- function(comb_sample_Tables, savecopy = TRUE) {

  snp_clss <- unique(comb_sample_Tables$SNP)
  snp_test_results <- list()
  snp_cls_names <- c()

  for (snp_cls in snp_clss) {
    pheno_data <- comb_sample_Tables[comb_sample_Tables$SNP == snp_cls, ]
    pheno_data$comb <- as.factor(pheno_data$comb)

    tryCatch({
      if (length(unique(pheno_data$comb)) > 1 && nrow(pheno_data) > 2) {
        stat.test <- t_test(pheno_data, Pheno ~ comb)
        print(stat.test)  # Debug statement to print the t_test result

        # Ensure the p.adj column exists before proceeding
        if (!"p.adj" %in% colnames(stat.test)) {
          # skip this SNP combination if p.adj column is missing
          next
        }

        stat.test <- add_xy_position(stat.test, x = "Pheno")
        stat.test <- mutate(stat.test, xmin = as.numeric(factor(group1, levels = levels(pheno_data$comb))),
                            xmax = as.numeric(factor(group2, levels = levels(pheno_data$comb))))
        stat.test$snp <- snp_cls

        if (nrow(stat.test) > 0) {
          snp_test_results <- append(snp_test_results, list(stat.test))
          snp_cls_names <- append(snp_cls_names, snp_cls)
        }
      }
    }, error = function(e) {
      message(paste("Error in SNP", snp_cls, ": ", e$message))
    })
  }
  # if no results, return NULL
  if (length(snp_test_results) == 0) {
    return(NULL)
  }

  names(snp_test_results) <- snp_cls_names
  if (savecopy) {
    saveTTestResultsToFile(snp_test_results)
  }

  return(snp_test_results)
}
