# Guess the column separator from the first line of a text file
detect_sep <- function(file) {
  first <- readLines(file, n = 1, warn = FALSE)
  if (grepl("\t", first)) "\t" else if (grepl(",", first)) "," else if (grepl(";", first)) ";" else ""
}

read_table_auto <- function(file, sep = NULL, header = TRUE, ...) {
  if (!file.exists(file)) stop("File not found: ", file)
  if (is.null(sep)) sep <- detect_sep(file)
  utils::read.table(file, sep = sep, header = header, check.names = FALSE,
                    stringsAsFactors = FALSE, quote = "\"", comment.char = "",
                    na.strings = c("NA", ""), ...)
}

#' Read GWAS results file with specified columns
#'
#' The column names are taken from the configuration (`rsid_col`, `chr_col`,
#' `pos_col`, `pval_col`); see [set_config()].
#' @param gwasfile Path to the GWAS file (tab, comma or semicolon separated; may be gzipped)
#' @param header Logical, whether the file contains a header row
#' @param sep Separator used in the file. `NULL` (default) detects it automatically.
#' @return A data frame with the standardised columns `rsid`, `chr`, `pos`, `p`
#'   and `rs` (`chr:pos` identifier).
#' @export
#' @examples
#' gwas <- read_gwas_file(system.file("extdata", "gwas_area_ann19.csv.gz", package = "HaploTraitR"))
#' head(gwas)
read_gwas_file <- function(gwasfile, header = TRUE, sep = NULL) {
  cols <- c(rsid = get_config("rsid_col"), pos = get_config("pos_col"),
            chr = get_config("chr_col"), p = get_config("pval_col"))
  gwas <- read_table_auto(gwasfile, sep = sep, header = header)

  missing_cols <- cols[!cols %in% colnames(gwas)]
  if (length(missing_cols) > 0) {
    stop("The following required columns are missing from the GWAS file: ",
         paste(missing_cols, collapse = ", "), "\nColumns found: ",
         paste(colnames(gwas), collapse = ", "),
         "\nUse set_config(list(rsid_col = ..., chr_col = ..., pos_col = ..., pval_col = ...)).")
  }
  for (std in names(cols)) colnames(gwas)[colnames(gwas) == cols[[std]]] <- std

  gwas$chr <- as.character(gwas$chr)
  gwas$pos <- as.numeric(gwas$pos)
  gwas$p <- as.numeric(gwas$p)
  gwas <- gwas[!is.na(gwas$pos) & !is.na(gwas$p), ]
  if (!"rs" %in% colnames(gwas)) gwas$rs <- paste(gwas$chr, gwas$pos, sep = ":")
  gwas
}

#' Filter GWAS data by FDR threshold
#'
#' @param gwas Data frame with the GWAS data (from [read_gwas_file()])
#' @return A list of data frames, one per chromosome, with the significant SNPs.
#' @export
#' @examples
#' gwas <- read_gwas_file(system.file("extdata", "gwas_area_ann19.csv.gz", package = "HaploTraitR"))
#' set_config(list(fdr_threshold = 0.1))
#' sig <- filter_gwas_data(gwas)
#' reset_config()
filter_gwas_data <- function(gwas) {
  fdr_col <- get_config("fdr_col")
  fdr_threshold <- get_config("fdr_threshold")

  if (!fdr_col %in% colnames(gwas)) {
    gwas[[fdr_col]] <- p.adjust(gwas$p, method = get_config("fdr_method"))
    cli::cli_alert_info("FDR column not found; FDR computed with method {.val {get_config('fdr_method')}}.")
  }
  gwas <- gwas[!is.na(gwas[[fdr_col]]) & gwas[[fdr_col]] < fdr_threshold, ]
  n_sig <- nrow(gwas)
  cli::cli_alert_info("{.strong {n_sig}} significant GWAS SNP{?s} (FDR < {fdr_threshold}).")
  if (nrow(gwas) == 0) {
    warning("No GWAS SNP passed the FDR threshold. Consider set_config(list(fdr_threshold = 0.1)).")
  }
  split(gwas, gwas$chr)
}

# Standardise genotype calls to two-letter codes (AA, AG, ..., NN for missing)
normalize_calls <- function(x) {
  iupac <- c(A = "AA", C = "CC", G = "GG", T = "TT", R = "AG", Y = "CT", S = "CG",
             W = "AT", K = "GT", M = "AC", N = "NN", "-" = "NN", "0" = "NN")
  x <- toupper(x)
  x <- gsub("[/|: ]", "", x)
  one <- !is.na(x) & nchar(x) == 1
  x[one] <- iupac[x[one]]
  bad <- is.na(x) | nchar(x) != 2 | grepl("[^ACGT]", x)
  x[bad] <- "NN"
  x
}

#' Read a HapMap genotype file
#'
#' Reads a HapMap file (11 annotation columns followed by one column per
#' accession). Genotypes can be two-letter (`AA`, `AG`) or single-letter IUPAC
#' codes (`A`, `R`); missing calls (`NN`, `N`, `-`) are standardised to `NN`.
#' @param path Path to the HapMap file (may be gzipped)
#' @param sep Separator (default tab)
#' @param header Logical indicating if the file has a header
#' @param ... Additional arguments passed to [utils::read.table()]
#' @return A list of data frames, one per chromosome, with row names `chr:pos`.
#' @export
#' @examples
#' hapmap <- readHapmap(system.file("extdata", "Barley_50K.tsv.gz", package = "HaploTraitR"))
#' names(hapmap)
readHapmap <- function(path, sep = "\t", header = TRUE, ...) {
  hap <- read_table_auto(path, sep = sep, header = header, colClasses = "character", ...)
  if (ncol(hap) < 12) stop("A HapMap file needs 11 annotation columns plus at least one accession.")
  colnames(hap)[3:4] <- c("chrom", "pos")
  hap$pos <- as.numeric(hap$pos)

  geno <- as.matrix(hap[, -(1:11), drop = FALSE])
  if (any(nchar(geno) != 2, na.rm = TRUE) || any(is.na(geno)) || any(grepl("[^ACGTN]", geno))) {
    geno[] <- normalize_calls(geno)
    hap[, -(1:11)] <- as.data.frame(geno, stringsAsFactors = FALSE)
  }

  rsid <- paste(hap$chrom, hap$pos, sep = ":")
  if (any(duplicated(rsid))) {
    n_dup <- sum(duplicated(rsid))
    cli::cli_alert_info("{n_dup} SNP{?s} share{?s/} a position with another SNP; only the first is kept.")
    hap <- hap[!duplicated(rsid), ]
    rsid <- rsid[!duplicated(rsid)]
  }
  rownames(hap) <- rsid
  hap <- hap[order(hap$chrom, hap$pos), ]
  split(hap, hap$chrom)
}

#' Read a VCF file into HapMap format
#'
#' Converts the `GT` field of a (possibly gzipped) VCF file into the HapMap
#' structure used by HaploTraitR, so VCF data can be used directly in the
#' pipeline. Only bi- and multi-allelic SNPs are kept (indels are skipped);
#' heterozygous and missing calls are preserved (`AG`, `NN`).
#' @param vcfpath Path to the VCF file
#' @return A list of data frames, one per chromosome, as returned by [readHapmap()].
#' @export
readVCF <- function(vcfpath) {
  lines <- readLines(vcfpath)
  lines <- lines[!startsWith(lines, "##")]
  header <- strsplit(sub("^#", "", lines[1]), "\t")[[1]]
  body <- strsplit(lines[-1], "\t")
  vcf <- do.call(rbind, body)
  colnames(vcf) <- header
  ref <- vcf[, "REF"]
  alt <- vcf[, "ALT"]
  snp <- nchar(ref) == 1 & grepl("^[ACGT](,[ACGT])*$", alt)
  vcf <- vcf[snp, , drop = FALSE]
  if (nrow(vcf) == 0) stop("No SNPs found in the VCF file.")

  samples <- header[10:length(header)]
  alleles <- strsplit(paste(vcf[, "REF"], vcf[, "ALT"], sep = ","), ",")
  gt_index <- vapply(strsplit(vcf[, "FORMAT"], ":"), function(f) match("GT", f), integer(1))
  geno <- matrix("NN", nrow(vcf), length(samples), dimnames = list(NULL, samples))
  for (i in seq_len(nrow(vcf))) {
    gt <- vapply(strsplit(vcf[i, samples], ":"), `[`, character(1), gt_index[i])
    idx <- strsplit(gt, "[/|]")
    geno[i, ] <- vapply(idx, function(a) {
      a <- suppressWarnings(as.integer(a))
      if (length(a) == 1) a <- c(a, a)
      if (length(a) != 2 || anyNA(a)) "NN" else paste0(alleles[[i]][a + 1], collapse = "")
    }, character(1))
  }
  ids <- ifelse(vcf[, "ID"] %in% c(".", ""), paste(vcf[, "CHROM"], vcf[, "POS"], sep = ":"), vcf[, "ID"])
  hap <- data.frame(`rs#` = ids, alleles = paste(vcf[, "REF"], vcf[, "ALT"], sep = "/"),
                    chrom = vcf[, "CHROM"], pos = as.numeric(vcf[, "POS"]), strand = "+",
                    `assembly#` = NA, center = NA, protLSID = NA, assayLSID = NA,
                    panelLSID = NA, QCcode = NA, check.names = FALSE, stringsAsFactors = FALSE)
  hap <- cbind(hap, as.data.frame(geno, stringsAsFactors = FALSE))
  rsid <- paste(hap$chrom, hap$pos, sep = ":")
  hap <- hap[!duplicated(rsid), ]
  rownames(hap) <- rsid[!duplicated(rsid)]
  hap <- hap[order(hap$chrom, hap$pos), ]
  split(hap, hap$chrom)
}

#' Prepare phenotype data
#'
#' Selects the trait column, removes missing values and averages replicated
#' measurements of the same accession.
#' @param pheno A data frame whose first column holds accession names and
#'   other columns hold traits. If it has more than two columns the trait is
#'   taken from the `pheno_col` configuration parameter.
#' @return A data frame with columns `Sample` and `Pheno`.
#' @export
#' @examples
#' pheno <- data.frame(Taxa = c("G1", "G1", "G2"), Yield = c(1, 3, 5))
#' prepare_pheno(pheno)
prepare_pheno <- function(pheno) {
  pheno <- as.data.frame(pheno)
  if (ncol(pheno) > 2) {
    col <- get_config("pheno_col")
    if (!col %in% colnames(pheno)) {
      stop("Trait column '", col, "' not found in the phenotype data. Available columns: ",
           paste(colnames(pheno)[-1], collapse = ", "),
           "\nUse set_config(list(pheno_col = \"<column>\")).")
    }
    values <- pheno[[col]]
  } else {
    values <- pheno[[2]]
  }
  out <- data.frame(Sample = as.character(pheno[[1]]),
                    Pheno = suppressWarnings(as.numeric(values)),
                    stringsAsFactors = FALSE)
  out <- out[!is.na(out$Pheno) & !is.na(out$Sample), ]
  if (anyDuplicated(out$Sample)) {
    cli::cli_alert_info("Replicated measurements found; the mean per accession is used.")
    out <- stats::aggregate(Pheno ~ Sample, data = out, FUN = mean)
  }
  out
}

#' Create a unique result folder for HaploTraitR
#' @param base_folder_name The base name for the result folder
#' @param location The location where the folder should be created
#'   (default: the current working directory)
#' @return The path to the created folder
#' @export
#' @examples
#' result_folder <- create_unique_result_folder(location = tempdir())
create_unique_result_folder <- function(base_folder_name = "haplotraitR_run", location = NULL) {
  if (is.null(location)) {
    location <- getwd()
    cli::cli_alert_info("No location provided. Using the working directory {.path {location}}.")
  }
  if (!dir.exists(location)) dir.create(location, recursive = TRUE)
  result_folder <- file.path(location, base_folder_name)
  counter <- 1
  while (dir.exists(result_folder)) {
    counter <- counter + 1
    result_folder <- file.path(location, paste(base_folder_name, counter, sep = "_"))
  }
  dir.create(result_folder, showWarnings = FALSE)
  result_folder
}
