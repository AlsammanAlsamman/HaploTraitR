#' Convert HapMap genotypes to allele dosages
#'
#' @param habmap A HapMap data frame (11 annotation columns followed by
#'   genotype columns), or a character matrix of genotype calls.
#' @return An integer matrix (SNPs x accessions) counting copies of the minor
#'   allele (0, 1, 2); missing calls are `NA`.
#' @export
#' @examples
#' g <- data.frame(matrix(NA, 2, 11), G1 = c("AA", "CC"), G2 = c("AG", "CT"), G3 = c("GG", "NN"))
#' convertGenoBi2Numeric(g)
convertGenoBi2Numeric <- function(habmap) {
  geno <- if (is.matrix(habmap)) habmap else as.matrix(habmap[, -(1:11), drop = FALSE])
  a1 <- substr(geno, 1, 1)
  a2 <- substr(geno, 2, 2)
  missing <- a1 == "N" | a2 == "N" | is.na(geno)
  minor <- vapply(seq_len(nrow(geno)), function(i) {
    alleles <- c(a1[i, !missing[i, ]], a2[i, !missing[i, ]])
    if (length(alleles) == 0) return(NA_character_)
    tab <- sort(table(alleles), decreasing = TRUE)
    if (length(tab) > 1) names(tab)[2] else names(tab)[1]
  }, character(1))
  dosage <- (a1 == minor) + (a2 == minor)
  dosage[missing] <- NA
  storage.mode(dosage) <- "integer"
  dimnames(dosage) <- dimnames(geno)
  dosage
}

# Squared correlation between genotype dosages (r2, unphased LD)
ld_r2 <- function(dosage) {
  r2 <- suppressWarnings(stats::cor(t(dosage), use = "pairwise.complete.obs")^2)
  diag(r2) <- 1
  r2
}

# Major and minor allele of every SNP (rows of a genotype call matrix)
allele_summary <- function(geno) {
  out <- t(vapply(seq_len(nrow(geno)), function(i) {
    calls <- geno[i, geno[i, ] != "NN"]
    alleles <- c(substr(calls, 1, 1), substr(calls, 2, 2))
    if (length(alleles) == 0) return(c(NA_character_, NA_character_, NA_character_))
    tab <- sort(table(alleles), decreasing = TRUE)
    minor <- if (length(tab) > 1) names(tab)[2] else NA_character_
    maf <- if (length(tab) > 1) tab[[2]] / sum(tab) else 0
    c(names(tab)[1], minor, format(round(maf, 4)))
  }, character(3)))
  data.frame(major = out[, 1], minor = out[, 2], maf = as.numeric(out[, 3]),
             stringsAsFactors = FALSE)
}

#' Find nearby SNPs to a given SNP list assuming they are in the same chromosome
#' @param snp_list1 A vector of SNP positions (reference)
#' @param snp_list2 A vector of SNP positions (query)
#' @param threshold The maximum distance between the SNPs
#' @return A list with, for each position in `snp_list1`, the positions of
#'   `snp_list2` within `threshold`.
#' @export
#' @examples
#' find_nearby_snps(c(1000, 2000, 3000), c(1500, 2500, 3500, 5000), 1000)
find_nearby_snps <- function(snp_list1, snp_list2, threshold) {
  lapply(snp_list1, function(p) snp_list2[abs(snp_list2 - p) <= threshold])
}

#' Extract HapMap data for significant SNPs
#' @param hapmap A list of HapMap data frames (from [readHapmap()])
#' @param gwas A list of GWAS data frames (from [filter_gwas_data()])
#' @return A list of HapMap data frames containing only the GWAS SNPs.
#' @export
extract_hapmap <- function(hapmap, gwas) {
  sub_hapmap <- list()
  for (chromosome in names(gwas)) {
    if (is.null(hapmap[[chromosome]])) next
    ids <- paste(gwas[[chromosome]]$chr, gwas[[chromosome]]$pos, sep = ":")
    ids <- ids[ids %in% rownames(hapmap[[chromosome]])]
    if (length(ids)) sub_hapmap[[chromosome]] <- hapmap[[chromosome]][ids, , drop = FALSE]
  }
  sub_hapmap
}

#' Get the phenotype and genotype data
#' @param subhapmap A list of HapMap data frames (from [extract_hapmap()])
#' @param pheno A phenotype data frame (accession, trait)
#' @return A long data frame with columns `Phenotype`, `variable` (SNP),
#'   `value` (genotype call) and `Sample`.
#' @export
get_pheno_geno <- function(subhapmap, pheno) {
  pheno <- prepare_pheno(pheno)
  names(subhapmap) <- NULL
  geno <- do.call(rbind, subhapmap)
  geno <- as.matrix(geno[, -(1:11), drop = FALSE])
  samples <- intersect(colnames(geno), pheno$Sample)
  geno <- geno[, samples, drop = FALSE]
  data.frame(
    Phenotype = rep(pheno$Pheno[match(samples, pheno$Sample)], times = nrow(geno)),
    variable = rep(rownames(geno), each = length(samples)),
    value = as.vector(t(geno)),
    Sample = rep(samples, times = nrow(geno)),
    stringsAsFactors = FALSE
  )
}
