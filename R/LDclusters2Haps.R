#' Convert LD clusters to haplotypes and calculate the frequency of each haplotype combination
#'
#' @param hapmap a hapmap object containing the genotypic data
#' @param clusterLDs a list of clusters
#' @param savecopy a boolean indicating whether to save a copy of the haplotype combinations to the project folder
#' @return a data frame containing the haplotype combinations and their frequencies
#' @export
convertLDclusters2Haps <- function(hapmap, clusterLDs, savecopy=TRUE) {
  # Get the configuration
  comb_freq_threshold <- get_config("comb_freq_threshold")

  # Create a list to store the haplotype combinations
  hap_comb <- list()

  # Loop over the clusters
  for (i in seq_along(clusterLDs)) {
    cluster <- clusterLDs[[i]]

    # If there is no cluster, skip
    if (length(cluster) == 0) {
      next
    }

    # Get the chromosome
    chr <- strsplit(cluster[1], ":")[[1]][1]

    # Select the data of these SNPs
    cluster.data <- hapmap[[chr]][cluster, ]

    # Remove all non-variant columns (first 11 columns)
    cluster.data <- cluster.data[, -c(1:11)]

    # Transpose the data
    cluster.data <- t(cluster.data)

    # Concatenate the data in one column
    clusterComb <- apply(cluster.data, 1, paste, collapse = "|")

    # Calculate the frequency of each combination
    snpsComb <- as.data.frame(table(clusterComb))

    totalsample <- sum(snpsComb$Freq)
    snpsComb$Freq <- snpsComb$Freq / totalsample
    snpsComb$snps <- paste(cluster, collapse = "|")
    snpsComb$chr <- chr
    snpsComb$snp <- names(clusterLDs)[i]
    # Remove rows with frequency less than comb_freq_threshold
    snpsComb <- snpsComb[snpsComb$Freq > comb_freq_threshold, ]
    # if there is no combination, skip
    if (nrow(snpsComb) == 0) {
      next
    }
    snpsComb$comb <- 1:nrow(snpsComb)

    # Append to hap_comb
    hap_comb[[length(hap_comb) + 1]] <- snpsComb
  }

  # Combine all data frames in the list into one data frame
  if (length(hap_comb) > 0) {
    hap_comb <- do.call(rbind, hap_comb)
  } else {
    hap_comb <- data.frame()
  }
  # save a copy to the project folder
  if (savecopy) {
    outfolder <- get_config("outfolder")
    if (!is.null(outfolder) && !dir.exists(outfolder)) {
      print("The output folder does not exist, please create it first, meanwhile the haplotype combinations will not be saved.")
      return(hap_comb)
    }
    write.csv(hap_comb, file.path(outfolder, "hap_comb.csv"), row.names = FALSE)
    print("Haplotype combinations saved to hap_comb.csv")
  }
  return(hap_comb)
}
