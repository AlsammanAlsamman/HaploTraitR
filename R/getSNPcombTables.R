#' Get SNP combinations sample tables
#' @param haplotypes A list of haplotype combinations
#' @param hapmap A hapmap object
#' @return A list of sample tables for each SNP combination
#' @export
getSNPcombTables <- function(haplotypes, hapmap) {
  snps <- unique(snpcombsample$snp)
  hap_comb_samples <- list()
  snpadded<-c()
  comb_sample_Tables<-list()

  for (snp in snps) {
    # Get the haplotype combination
    hap_comb <- snpcombsample[snpcombsample$snp == snp, ]

    # Initialize comb_sample as a data frame with specified column types
    comb_sample <- data.frame(Sample = character(),
                              Pheno = numeric(),
                              Comb = character(),
                              SNP = character(),
                              stringsAsFactors = FALSE)

    # Loop to fill in comb_sample
    for (i in 1:nrow(hap_comb)) {
      comb <- hap_comb[i, ]
      samples <- unlist(strsplit(comb$samples, "\\|"))

      # Subset pheno data for matching samples
      pheno_samples <- pheno[pheno[, 1] %in% samples, ]

      # Only proceed if pheno_samples has matching rows
      if (nrow(pheno_samples) > 0) {
        # Add the "Comb" and "SNP" columns
        pheno_samples$Comb <- paste0("comb_", i)
        pheno_samples$SNP <- comb$snp

        # Rename columns to match exactly with comb_sample
        colnames(pheno_samples) <- c("Sample", "Pheno", "Comb", "SNP")

        # Append to comb_sample
        comb_sample <- rbind(comb_sample, pheno_samples)
      }
    }
    comb_sample_Tables[[snp]]<-comb_sample

  }
  return(comb_sample_Tables)
}
