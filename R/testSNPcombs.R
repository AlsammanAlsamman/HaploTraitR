#' Calculate the p-value for a SNP combinations
#' @param comb_sample_Tables A list of sample tables for each SNP combination
#' @return A data frame with the p-values for the SNP combinations
#' @importFrom magrittr %>%
#' @importFrom rstatix t_test
#' @importFrom rstatix add_xy_position
#' @export
testSNPcombs <- function(comb_sample_Tables) {
  snps <- names(comb_sample_Tables)
  hap_comb_tests <- list()
  snpadded<-c()
  for (snp in snps) {
    # Initialize comb_sample as a data frame with specified column types
    comb_sample <- comb_sample_Tables[[snp]]
    # Convert "Comb" to factor
    comb_sample$Comb <- factor(comb_sample$Comb)
    # Perform statistical test if there are enough samples
    tryCatch({
      # Ensure there are at least two levels in "Comb" and sufficient data for testing
      if (length(unique(comb_sample$Comb)) > 1 && nrow(comb_sample) > 2) {
        stat.test <- comb_sample %>% t_test(Pheno ~ Comb)
        stat.test <- stat.test %>% add_xy_position(x = "Pheno")
        # Set xmin and xmax based on the group1 and group2 column values
        stat.test <- stat.test %>%
          mutate(xmin = as.numeric(factor(group1, levels = levels(comb_sample$Comb))),
                 xmax = as.numeric(factor(group2, levels = levels(comb_sample$Comb))))
        stat.test$snp <- snp
        if (nrow(stat.test) > 0) {
          hap_comb_tests <- append(hap_comb_tests, list(stat.test))
          snpadded<-append(snpadded, snp)
        }
      }
    }, error = function(e) {
      message(paste("Error in SNP", snp, ": ", e$message))

    })
  }

  names(hap_comb_tests) <- snpadded
  return(hap_comb_tests)
}
