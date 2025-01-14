#' Read GWAS results file with specified columns
#'
#' @param gwasfile Path to the GWAS file
#' @param header Logical, whether the file contains a header row
#' @param sep Separator used in the file
#' @return A data frame with the GWAS data
#' @export
#' @examples
#' gwasfile <- "SignificantSNP_GWAS.csv"
#' gwas <- read_gwas_file(gwasfile)
#' print(gwas)
read_gwas_file <- function(gwasfile, header = TRUE, sep = "\t") {
  # Access configuration parameters from the environment
  rsid_col <- get_config("rsid_col")
  pos_col <- get_config("pos_col")
  chr_col <- get_config("chr_col")
  pval_col <- get_config("pval_col")

  # Read the GWAS file
  gwas <- read.csv(gwasfile, header = header, sep = sep)

  # Check if the required columns are present
  missing_cols <- c()
  if (!(rsid_col %in% colnames(gwas))) missing_cols <- c(missing_cols, rsid_col)
  if (!(pos_col %in% colnames(gwas))) missing_cols <- c(missing_cols, pos_col)
  if (!(chr_col %in% colnames(gwas))) missing_cols <- c(missing_cols, chr_col)
  if (!(pval_col %in% colnames(gwas))) missing_cols <- c(missing_cols, pval_col)
  if (length(missing_cols) > 0) {
    stop(paste("The following required columns are missing from the GWAS file:", paste(missing_cols, collapse = ", ")))
  }

  # Rename columns according to the specified configuration variables
  colnames(gwas)[colnames(gwas) == rsid_col] <- "rsid"
  colnames(gwas)[colnames(gwas) == pos_col] <- "pos"
  colnames(gwas)[colnames(gwas) == chr_col] <- "chr"
  colnames(gwas)[colnames(gwas) == pval_col] <- "p"

  # Add rs id if not already present
  if (!"rs" %in% colnames(gwas)) {
    gwas$rs <- paste(gwas$chr, gwas$pos, sep = ":")
  }

  return(gwas)
}

#' Filter GWAS data by FDR threshold
#'
#' @param gwas Data frame with the GWAS data
#' @return A list of data frames split by chromosome with filtered GWAS data
#' @importFrom stats p.adjust
#' @export
#' @examples
#' gwasfile <- "SignificantSNP_GWAS.csv"
#' gwas <- read_gwas_file(gwasfile)
#' filtered_gwas <- filter_gwas_data(gwas)
#' print(filtered_gwas)
filter_gwas_data <- function(gwas) {
  # Access configuration parameters from the environment
  fdr_col <- get_config("fdr_col")
  fdr_threshold <- get_config("fdr_threshold")
  fdr_method <- get_config("fdr_method")

  # Calculate FDR if the column is not present
  if (!fdr_col %in% colnames(gwas)) {
    gwas[[fdr_col]] <- p.adjust(gwas$p, method = fdr_method)
    message("FDR column not found, calculated FDR")
    message("If you want to use a different method, please specify it in the configuration parameters")
  }

  # Filter by FDR threshold
  gwas <- gwas[gwas[[fdr_col]] < fdr_threshold, ]
  message(paste("SNPs with FDR <", fdr_threshold, ":", nrow(gwas)))
  message(paste("Number of significant SNPs:", nrow(gwas)))

  # Split by chromosome
  gwas <- split(gwas, gwas$chr)

  return(gwas)
}
