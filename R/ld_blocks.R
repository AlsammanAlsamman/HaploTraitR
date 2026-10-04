#' Get distance clusters around the significant GWAS SNPs
#'
#' For every significant GWAS SNP, collects the genotyped SNPs within
#' `dist_threshold` bp. Windows with fewer than `dist_cluster_count` SNPs are dropped.
#' @param gwas A list of GWAS data frames split by chromosome (from [filter_gwas_data()])
#' @param hapmap A list of HapMap data frames split by chromosome (from [readHapmap()])
#' @return A list (per chromosome) of named lists (per GWAS SNP) of SNP positions.
#' @export
getDistClusters <- function(gwas, hapmap) {
  dist_threshold <- get_config("dist_threshold")
  min_count <- get_config("dist_cluster_count")
  dist_clusters <- list()
  for (chromosome in names(gwas)) {
    if (is.null(hapmap[[chromosome]])) {
      cli::cli_alert_warning("Chromosome {.val {chromosome}} is not in the genotype file; skipped.")
      next
    }
    nearby <- find_nearby_snps(gwas[[chromosome]]$pos, hapmap[[chromosome]]$pos, dist_threshold)
    names(nearby) <- gwas[[chromosome]]$rs
    nearby <- nearby[lengths(nearby) > min_count]
    if (length(nearby) == 0) {
      cli::cli_alert_warning("No SNP window with more than {min_count} SNPs on chromosome {.val {chromosome}}.")
      next
    }
    dist_clusters[[chromosome]] <- nearby
  }
  dist_clusters
}

#' Compute LD matrices for the distance clusters
#'
#' Computes the r2 (squared correlation of allele dosages) between all SNPs of
#' each window and saves the matrices in the `LD_matrices` folder.
#' @param hapmap A list of HapMap data frames (from [readHapmap()])
#' @param haplotype_clusters A list of distance clusters (from [getDistClusters()])
#' @return A list with the LD folder (`ld_folder`), cluster information
#'   (`out_info`) and the LD matrices (`matrices`).
#' @export
computeLDclusters <- function(hapmap, haplotype_clusters) {
  outfolder <- get_outfolder(required = FALSE)
  ld_folder <- if (!is.null(outfolder)) out_subdir("LD_matrices") else NULL

  out_info <- list()
  matrices <- list()
  for (chr in names(haplotype_clusters)) {
    for (cls in names(haplotype_clusters[[chr]])) {
      ids <- paste(chr, haplotype_clusters[[chr]][[cls]], sep = ":")
      dosage <- convertGenoBi2Numeric(hapmap[[chr]][ids, , drop = FALSE])
      ld <- ld_r2(dosage)
      dimnames(ld) <- list(ids, ids)
      matrices[[cls]] <- ld
      file <- paste0(safe_name(cls), "_ld_matrix.csv")
      out_info[[length(out_info) + 1]] <- c(chr, cls, file)
      if (!is.null(ld_folder)) utils::write.csv(ld, file.path(ld_folder, file))
    }
  }
  out_info <- do.call(rbind, out_info)
  if (!is.null(out_info)) colnames(out_info) <- c("chr", "snp", "file")
  if (!is.null(ld_folder) && !is.null(out_info)) {
    utils::write.csv(out_info, file.path(ld_folder, "LD_index.csv"), row.names = FALSE)
  }
  list(ld_folder = ld_folder, out_info = out_info, matrices = matrices)
}

#' Get LD matrix information from a folder
#'
#' @param folder The folder where the LD matrices are stored (written by [computeLDclusters()])
#' @return A list as returned by [computeLDclusters()].
#' @export
retrieveLDMatricesFromFolder <- function(folder) {
  if (!dir.exists(folder)) stop("The specified folder does not exist: ", folder)
  index <- file.path(folder, "LD_index.csv")
  if (file.exists(index)) {
    out_info <- as.matrix(utils::read.csv(index, stringsAsFactors = FALSE))
  } else {
    # files written by HaploTraitR < 0.2 are named "<chr>:<pos>_ld_matrix.csv"
    files <- list.files(folder, pattern = "_ld_matrix\\.csv$")
    snp <- sub("_ld_matrix\\.csv$", "", files)
    out_info <- cbind(chr = sub(":.*$", "", snp), snp = snp, file = files)
  }
  matrices <- lapply(out_info[, "file"], function(f) {
    m <- as.matrix(utils::read.csv(file.path(folder, f), row.names = 1, check.names = FALSE))
    colnames(m) <- rownames(m)
    m
  })
  names(matrices) <- out_info[, "snp"]
  list(ld_folder = folder, out_info = out_info, matrices = matrices)
}

# Connected components of a logical adjacency matrix
ld_components <- function(adj) {
  n <- nrow(adj)
  membership <- rep(NA_integer_, n)
  k <- 0L
  for (i in seq_len(n)) {
    if (!is.na(membership[i])) next
    k <- k + 1L
    membership[i] <- k
    queue <- i
    while (length(queue)) {
      v <- queue[1]
      queue <- queue[-1]
      nb <- which(adj[v, ] & is.na(membership))
      membership[nb] <- k
      queue <- c(queue, nb)
    }
  }
  membership
}

#' Cluster SNPs based on LD matrices
#'
#' SNPs are linked when their r2 is at least `ld_threshold`; the connected
#' group of SNPs that contains the GWAS SNP forms its LD block. Blocks smaller
#' than `dist_cluster_count` SNPs are discarded.
#' @param LDsInfo A list returned by [computeLDclusters()] or [retrieveLDMatricesFromFolder()]
#' @return A named list (per GWAS SNP) of the SNP identifiers in its LD block.
#' @export
getLDclusters <- function(LDsInfo) {
  ld_threshold <- get_config("ld_threshold")
  min_count <- get_config("dist_cluster_count")
  matrices <- LDsInfo$matrices
  if (is.null(matrices)) matrices <- retrieveLDMatricesFromFolder(LDsInfo$ld_folder)$matrices

  cluster_info <- list()
  for (lead in names(matrices)) {
    ld <- matrices[[lead]]
    if (!lead %in% rownames(ld)) {
      cli::cli_alert_warning("GWAS SNP {.val {lead}} is not genotyped; its window is skipped.")
      cluster_info[[lead]] <- character(0)
      next
    }
    adj <- !is.na(ld) & ld >= ld_threshold
    membership <- ld_components(adj)
    block <- rownames(ld)[membership == membership[rownames(ld) == lead]]
    if (length(block) < min_count) block <- character(0)
    pos <- as.numeric(sub("^.*:", "", block))
    cluster_info[[lead]] <- block[order(pos)]
  }
  cluster_info
}
