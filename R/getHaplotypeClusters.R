#' Get haplotype clusters from GWAS and hapmap data By chromosome
#' @param gwas A list of GWAS data frames split by chromosome
#' @param hapmap A list of hapmap data frames split by chromosome
#' @param chromosome The chromosome of interest
#' @param dist_threshold The maximum distance between SNPs to consider for clustering
#' @param dist_cluster_count The minimum number of SNPs required to form a cluster
#' @return A list of haplotype clusters

getHaplotypeClustersByChr <- function(gwas, hapmap, chromosome, dist_threshold, dist_cluster_count) {
  # Select chromosome-specific data from GWAS and hapmap
  gwas_chr <- gwas[[chromosome]]
  hapmap_chr <- hapmap[[chromosome]]

  # Get the positions of the SNPs
  gwas_positions <- gwas_chr$pos
  hapmap_positions <- hapmap_chr$pos

  # Get nearby SNPs within the specified distance threshold
  nearby_snps <- find_nearby_snps(gwas_positions, hapmap_positions, dist_threshold)
  # check if there are any nearby SNPs
  if (length(nearby_snps) == 0) {
    message("No nearby SNPs found within the specified distance threshold.")
    return(NULL)
  }
  # name the list elements by snp rsid
  names(nearby_snps) <- gwas_chr$rs
  # Filter clusters by the specified minimum number of SNPs
  filtered_clusters <- nearby_snps[sapply(nearby_snps, length) > dist_cluster_count]
  # check if there are any filtered clusters
  if (length(filtered_clusters) == 0) {
    message("No clusters found with the specified minimum number of SNPs.")
    return(NULL)
  }
  # Return a list of clusters
  return(filtered_clusters)
}

#' Get haplotype clusters from GWAS and hapmap data
#' @param gwas A list of GWAS data frames split by chromosome
#' @param hapmap A list of hapmap data frames split by chromosome
#' @param dist_threshold The maximum distance between SNPs to consider for clustering
#' @param dist_cluster_count The minimum number of SNPs required to form a cluster
#' @return A list of haplotype clusters
#' @export

getHaplotypeClusters <- function(gwas, hapmap, dist_threshold, dist_cluster_count) {
  # Initialize an empty list to store haplotype clusters
  haplotype_clusters <- list()
  # Iterate over chromosomes
  for (chromosome in names(gwas)) {
    # Get haplotype clusters for the current chromosome
    clusters <- getHaplotypeClustersByChr(gwas, hapmap, chromosome, dist_threshold, dist_cluster_count)
    # Add the clusters to the list
    if (!is.null(clusters)) {
      haplotype_clusters[[chromosome]] <- clusters
    }
  }
  # Return the list of haplotype clusters
  return(haplotype_clusters)
}
