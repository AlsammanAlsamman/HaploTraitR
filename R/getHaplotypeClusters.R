#' Get Distance clusters from GWAS and hapmap data By chromosome
#' @param gwas A list of GWAS data frames split by chromosome
#' @param hapmap A list of hapmap data frames split by chromosome
#' @param chromosome The chromosome of interest
#' @return A list of Dist clusters
getDistClustersByChr <- function(gwas, hapmap, chromosome) {
  # Access configuration parameters from the environment
  dist_threshold <- get_config("dist_threshold")
  dist_cluster_count <- get_config("dist_cluster_count")

  # Select chromosome-specific data from GWAS and hapmap
  gwas_chr <- gwas[[chromosome]]
  hapmap_chr <- hapmap[[chromosome]]

  # Get the positions of the SNPs
  gwas_positions <- gwas_chr$pos
  hapmap_positions <- hapmap_chr$pos

  # Get nearby SNPs within the specified distance threshold
  nearby_snps <- find_nearby_snps(gwas_positions, hapmap_positions, dist_threshold)

  # Check if there are any nearby SNPs
  if (length(nearby_snps) == 0) {
    message(paste("No nearby SNPs found for chromosome", chromosome))
    return(NULL)
  }

  # Name the list elements by SNP rsid
  names(nearby_snps) <- gwas_chr$rs

  # Filter clusters by the specified minimum number of SNPs
  filtered_clusters <- nearby_snps[sapply(nearby_snps, length) > dist_cluster_count]

  # Check if there are any filtered clusters
  if (length(filtered_clusters) == 0) {
    message(paste("No filtered clusters found for chromosome", chromosome))
    return(NULL)
  }

  # Return a list of clusters
  return(filtered_clusters)
}

#' Get haplotype clusters from GWAS and hapmap data
#' @param gwas A list of GWAS data frames split by chromosome
#' @param hapmap A list of hapmap data frames split by chromosome
#' @return A list of haplotype clusters
#' @export
getDistClusters <- function(gwas, hapmap) {
  # Initialize an empty list to store haplotype clusters
  dist_clusters <- list()

  # Iterate over chromosomes
  for (chromosome in names(gwas)) {
    # Get haplotype clusters for the current chromosome
    clusters <- getDistClustersByChr(gwas, hapmap, chromosome)

    # Add the clusters to the list
    if (!is.null(clusters)) {
      dist_clusters[[chromosome]] <- clusters
    }
  }

  # Return the list of haplotype clusters
  return(dist_clusters)
}
