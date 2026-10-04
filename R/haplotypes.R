# Readable block identifier, e.g. "2H:579.59-581.37Mb"
make_block_id <- function(chr, pos) {
  paste0(chr, ":", sprintf("%.2f", min(pos) / 1e6), "-", sprintf("%.2f", max(pos) / 1e6), "Mb")
}

# Per-block SNP information stored as an attribute of the haplotype table
block_info <- function(haplotypes, block) {
  info <- attr(haplotypes, "blocks")
  if (is.null(info) || is.null(info[[block]])) {
    stop("No SNP information for block ", block,
         ". Use the haplotype table returned by convertLDclusters2Haps().")
  }
  info[[block]]
}

# Accept a block id or a GWAS lead SNP and return the block id
resolve_block <- function(haplotypes, id) {
  if (id %in% haplotypes$block) return(id)
  hit <- unique(haplotypes$block[vapply(strsplit(haplotypes$lead_snps, ";"), function(l) id %in% l, logical(1))])
  if (length(hit) == 0) stop("'", id, "' is neither a block id nor a GWAS SNP of a block.")
  hit[1]
}

# Keep the 'blocks' attribute when a haplotype table is modified
with_blocks <- function(df, from) {
  attr(df, "blocks") <- attr(from, "blocks")
  df
}

#' Convert LD clusters to haplotypes
#'
#' GWAS SNPs whose LD blocks contain exactly the same SNPs are merged into one
#' haplotype block. Within each block, every distinct combination of alleles
#' carried by at least `comb_freq_threshold` of the accessions becomes a
#' haplotype (H1 = most frequent).
#' @param hapmap A list of HapMap data frames (from [readHapmap()])
#' @param clusterLDs A list of LD blocks (from [getLDclusters()])
#' @param savecopy Save a copy (`tables/haplotypes_alleles.csv`) in the output folder
#' @return A data frame with one row per haplotype: `block`, `chr`, `start`,
#'   `end`, `n_snps`, `lead_snps`, `snps`, `hap`, `alleles` and `freq`.
#'   SNP details and LD of each block are kept in the `"blocks"` attribute.
#' @export
convertLDclusters2Haps <- function(hapmap, clusterLDs, savecopy = TRUE) {
  freq_thr <- get_config("comb_freq_threshold")
  clusterLDs <- clusterLDs[lengths(clusterLDs) > 0]
  if (length(clusterLDs) == 0) {
    warning("No LD block found. Try a lower ld_threshold or dist_cluster_count.")
    return(data.frame())
  }
  keys <- vapply(clusterLDs, paste, character(1), collapse = "|")

  rows <- list()
  info <- list()
  for (key in unique(keys)) {
    leads <- names(clusterLDs)[keys == key]
    snps <- clusterLDs[[which(keys == key)[1]]]
    chr <- sub(":.*$", "", snps[1])
    data <- hapmap[[chr]][snps, , drop = FALSE]
    pos <- as.numeric(data$pos)
    ord <- order(pos)
    data <- data[ord, , drop = FALSE]
    snps <- snps[ord]
    pos <- pos[ord]
    geno <- as.matrix(data[, -(1:11), drop = FALSE])

    combos <- apply(geno, 2, paste, collapse = "|")
    complete <- !grepl("NN", combos, fixed = TRUE)
    tab <- sort(table(combos[complete]), decreasing = TRUE)
    freq <- as.numeric(tab) / ncol(geno)
    keep <- freq > freq_thr
    if (!any(keep)) {
      cli::cli_alert_warning("Block of GWAS SNP{?s} {.val {leads}}: no haplotype above the frequency threshold; skipped.")
      next
    }
    block <- make_block_id(chr, pos)
    if (block %in% names(info)) block <- paste0(block, "_", sum(startsWith(names(info), block)) + 1)

    alleles <- allele_summary(geno)
    ld <- ld_r2(convertGenoBi2Numeric(geno))
    dimnames(ld) <- list(snps, snps)
    info[[block]] <- data.frame(snp = snps, marker = data[[1]], pos = pos,
                                alleles, stringsAsFactors = FALSE)
    attr(info[[block]], "ld") <- ld
    attr(info[[block]], "leads") <- leads

    rows[[block]] <- data.frame(
      block = block, chr = chr, start = min(pos), end = max(pos), n_snps = length(snps),
      lead_snps = paste(leads, collapse = ";"), snps = paste(snps, collapse = "|"),
      hap = paste0("H", seq_len(sum(keep))), alleles = names(tab)[keep], freq = freq[keep],
      stringsAsFactors = FALSE
    )
  }
  haplotypes <- do.call(rbind, rows)
  if (is.null(haplotypes)) return(data.frame())
  rownames(haplotypes) <- NULL
  attr(haplotypes, "blocks") <- info
  n_blocks <- length(info)
  n_leads <- length(clusterLDs)
  cli::cli_alert_info("{.strong {n_blocks}} haplotype block{?s} from {n_leads} GWAS SNP{?s}.")
  if (savecopy) save_table(haplotypes, "haplotypes_alleles.csv")
  haplotypes
}

#' Assign accessions to haplotypes
#'
#' Accessions whose alleles match a haplotype exactly are assigned to it.
#' Accessions with some missing calls (up to `hap_max_missing`) are assigned
#' to the single haplotype compatible with their observed alleles. All other
#' accessions are pooled as `"Other"`.
#' @param haplotypes A haplotype table (from [convertLDclusters2Haps()])
#' @param hapmap A list of HapMap data frames (from [readHapmap()])
#' @param savecopy Save a copy (`tables/haplotype_accessions.csv`) in the output folder
#' @return The haplotype table with an extra `"Other"` row per block and the
#'   columns `n` (number of accessions), `freq` and `samples` (`|`-separated).
#' @export
getHapCombSamples <- function(haplotypes, hapmap, savecopy = TRUE) {
  max_missing <- get_config("hap_max_missing")
  out <- list()
  for (block in unique(haplotypes$block)) {
    h <- haplotypes[haplotypes$block == block & haplotypes$hap != "Other", ]
    info <- block_info(haplotypes, block)
    geno <- as.matrix(hapmap[[h$chr[1]]][info$snp, -(1:11), drop = FALSE])
    hap_alleles <- do.call(cbind, strsplit(h$alleles, "|", fixed = TRUE))

    assigned <- apply(geno, 2, function(calls) {
      observed <- calls != "NN"
      if (mean(!observed) > max_missing) return("Other")
      match_hap <- colSums(hap_alleles[observed, , drop = FALSE] == calls[observed]) == sum(observed)
      if (sum(match_hap) == 1) h$hap[match_hap] else "Other"
    })

    groups <- c(h$hap, "Other")
    samples <- vapply(groups, function(g) paste(names(assigned)[assigned == g], collapse = "|"), character(1))
    n <- vapply(groups, function(g) sum(assigned == g), integer(1))
    other <- h[1, ]
    other$hap <- "Other"
    other$alleles <- ""
    h <- rbind(h, other)
    h$n <- n
    h$freq <- n / length(assigned)
    h$samples <- samples
    out[[block]] <- h[h$n > 0, ]
  }
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out <- with_blocks(out, haplotypes)
  if (savecopy) save_table(out[, setdiff(names(out), "snps")], "haplotype_accessions.csv")
  out
}

#' Build the accession-level haplotype and phenotype table
#'
#' @param haplotypes A haplotype table with samples (from [getHapCombSamples()])
#' @param pheno A phenotype data frame (accession, trait); see [prepare_pheno()]
#' @param savecopy Save a copy (`tables/accession_haplotype_phenotype.csv`) in the output folder
#' @return A long data frame with columns `Sample`, `Pheno`, `block`, `hap`
#'   and `lead_snps` (phenotyped accessions only).
#' @export
getSNPcombTables <- function(haplotypes, pheno, savecopy = TRUE) {
  pheno <- prepare_pheno(pheno)
  long <- do.call(rbind, lapply(seq_len(nrow(haplotypes)), function(i) {
    s <- strsplit(haplotypes$samples[i], "|", fixed = TRUE)[[1]]
    if (length(s) == 0) return(NULL)
    data.frame(Sample = s, block = haplotypes$block[i], hap = haplotypes$hap[i],
               lead_snps = haplotypes$lead_snps[i], stringsAsFactors = FALSE)
  }))
  n_geno <- length(unique(long$Sample))
  long$Pheno <- pheno$Pheno[match(long$Sample, pheno$Sample)]
  long <- long[!is.na(long$Pheno), c("Sample", "Pheno", "block", "hap", "lead_snps")]
  n_both <- length(unique(long$Sample))
  if (n_both == 0) {
    stop("No accession name is shared by the genotype and phenotype files. ",
         "Genotype examples: ", paste(head(unique(unlist(strsplit(haplotypes$samples, "|", fixed = TRUE))), 3), collapse = ", "),
         "; phenotype examples: ", paste(head(pheno$Sample, 3), collapse = ", "))
  }
  cli::cli_alert_info("{n_both} of {n_geno} genotyped accessions have phenotype data.")
  rownames(long) <- NULL
  if (savecopy) save_table(long, "accession_haplotype_phenotype.csv")
  long
}
