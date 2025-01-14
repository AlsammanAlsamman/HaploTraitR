#' Compute LD matrix for a cluster
#'
#' @param hapmap a hapmap object containing the genotypic data
#' @param haplotype_clusters a list of haplotype clusters
#' @param chr the chromosome
#' @param cls the cluster
#' @param outfolder the output folder
computeLDclusterByChr <- function(hapmap, haplotype_clusters, chr, cls, outfolder) {
  geno_cluster <- haplotype_clusters[[chr]][[cls]]

  # Extract the genotypic data of this cluster
  geno_cluster_data <- hapmap[[chr]][paste(chr, geno_cluster, sep = ":"), ]

  # Get the SNPs
  cluster_snps <- rownames(geno_cluster_data)

  # Convert the data to numeric
  geno_cluster_data_numeric <- convertGenoBi2Numeric(geno_cluster_data)

  # Calculate the LD matrix
  ld_matrix <- calculate_ld_matrix(geno_cluster_data_numeric)
  rownames(ld_matrix) <- cluster_snps
  colnames(ld_matrix) <- cluster_snps

  # Save the LD matrix
  ld_matrix_file <- file.path(outfolder, paste0(cls, "_ld_matrix.csv"))
  write.csv(ld_matrix, ld_matrix_file)
}

#' Compute LD matrix for clusters
#'
#' @param hapmap a hapmap object containing the genotypic data
#' @param haplotype_clusters a list of haplotype clusters
#' @return a list containing the path to the folder with LD matrices and information about the clusters
#' @export
computeLDclusters <- function(hapmap, haplotype_clusters) {
  # Get the output folder from the configuration
  outfolder <- get_config("outfolder")
  # check if the outfolder is not NULL and does exist
  if (!is.null(outfolder) && !dir.exists(outfolder)) {
    stop("The output folder does not exist, please create it first.")
  }
  # Create the "LD_matrices" folder inside the specified output folder
  ld_folder <- file.path(outfolder, "LD_matrices")
  if (!dir.exists(ld_folder)) {
    dir.create(ld_folder, showWarnings = FALSE, recursive = TRUE)
  }

  out_info <- list()

  for (chr in names(haplotype_clusters)) {
    for (cls in names(haplotype_clusters[[chr]])) {
      out_info[[length(out_info) + 1]] <- c(chr, cls)
      computeLDclusterByChr(hapmap, haplotype_clusters, chr, cls, ld_folder)
    }
  }

  # Convert the list to a data frame
  out_info <- do.call(rbind, out_info)

  # Return list
  out_list <- list(ld_folder = ld_folder, out_info = out_info)

  return(out_list)
}



#' Get LD matrix information from a folder
#'
#' @param folder The folder where the LD matrices are stored
#' @return A list containing the path to the folder with LD matrices and information about the clusters
#' @export
retrieveLDMatricesFromFolder <- function(folder) {
  # Check if the folder exists
  if (!dir.exists(folder)) {
    stop("The specified folder does not exist.")
  }

  # List all files in the folder
  files <- list.files(folder, full.names = TRUE)

  # Extract chromosome and cluster information from file names
  chr_cls <- do.call(rbind, strsplit(basename(files), ":"))
  chr_cls[, 2] <- gsub("_ld_matrix.csv", "", chr_cls[, 2])
  chr_cls[, 2] <- paste(chr_cls[, 1], chr_cls[, 2], sep = ":")

  # Create LD info list
  LD_info <- list(
    ld_folder = folder,
    out_info = chr_cls
  )
  return(LD_info)
}
