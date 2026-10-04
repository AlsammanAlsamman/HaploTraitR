# Add a formatted data sheet to an openxlsx workbook
add_sheet <- function(wb, name, df, widths = NULL, number_cols = NULL, digits = "0.00") {
  openxlsx::addWorksheet(wb, name, gridLines = FALSE)
  if (is.null(df) || nrow(df) == 0) {
    openxlsx::writeData(wb, name, data.frame(Note = "No results for this sheet."))
    return(invisible(wb))
  }
  header <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#1F4E78", fontColour = "white",
                                  halign = "center", valign = "center", wrapText = TRUE,
                                  border = "Bottom")
  openxlsx::writeData(wb, name, df, headerStyle = header, withFilter = TRUE)
  openxlsx::freezePane(wb, name, firstRow = TRUE)
  openxlsx::setRowHeights(wb, name, rows = 1, heights = 32)
  openxlsx::setColWidths(wb, name, cols = seq_len(ncol(df)), widths = widths %||% "auto")
  if (!is.null(number_cols)) {
    cols <- which(names(df) %in% number_cols)
    openxlsx::addStyle(wb, name, openxlsx::createStyle(numFmt = digits), rows = 2:(nrow(df) + 1),
                       cols = cols, gridExpand = TRUE, stack = TRUE)
  }
  zebra <- openxlsx::createStyle(fgFill = "#F3F6FA")
  even <- seq(3, nrow(df) + 1, by = 2)
  if (length(even)) {
    openxlsx::addStyle(wb, name, zebra, rows = even, cols = seq_len(ncol(df)), gridExpand = TRUE, stack = TRUE)
  }
  invisible(wb)
}

pretty_names <- function(df, map) {
  keep <- intersect(names(map), names(df))
  df <- df[, keep, drop = FALSE]
  names(df) <- unname(map[keep])
  df
}

#' Write the HaploTraitR results to a formatted Excel workbook
#'
#' The workbook contains the sheets *Read me*, *Block summary*,
#' *Haplotype effects*, *Diagnostic markers*, *Candidate accessions*,
#' *Pairwise tests*, *GWAS SNPs* and *Parameters*. Superior haplotypes are
#' highlighted in green.
#' @param results The list returned by [run_haplotraitr()] (or a list with the
#'   elements `summary`, `haplotypes`, `pairwise`, `markers`, `ranking`, `gwas_significant`)
#' @param file Path of the `.xlsx` file. Default: `<outfolder>/HaploTraitR_<trait>.xlsx`
#' @return Invisibly, the file path.
#' @export
write_excel_report <- function(results, file = NULL) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' is needed for the Excel report: install.packages(\"openxlsx\")")
  }
  trait <- get_config("phenotypename")
  unit <- get_config("phenotypeunit") %||% ""
  if (is.null(file)) file <- file.path(get_outfolder(), paste0("HaploTraitR_", safe_name(trait), ".xlsx"))
  wb <- openxlsx::createWorkbook()
  green <- openxlsx::createStyle(fgFill = "#C7E9C0", textDecoration = "bold")
  grey <- openxlsx::createStyle(fontColour = "#7F7F7F")

  # Read me ------------------------------------------------------------------
  s <- results$summary
  n_sig <- sum(s$blocks$strength != "None", na.rm = TRUE)
  readme <- c(
    paste("HaploTraitR results -", trait, if (nzchar(unit)) paste0("(", unit, ")") else ""),
    paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M"), "| HaploTraitR", as.character(utils::packageVersion("HaploTraitR"))),
    "",
    sprintf("%d haplotype block(s) analysed; %d with a significant haplotype effect.", nrow(s$blocks), n_sig),
    paste("Favourable direction:", if (isTRUE(get_config("higher_is_better"))) "higher trait values" else "lower trait values"),
    "",
    "SHEETS",
    "Block summary - one row per haplotype block (group of SNPs in LD around GWAS hits), with the superior haplotype and a selection recommendation.",
    "Haplotype effects - trait statistics of every haplotype: accessions, frequency, mean, 95% CI, difference from the population mean and letter groups.",
    "Haplotype alleles - the allele of every haplotype at every SNP of the block (major/minor allele as in the population).",
    "Diagnostic markers - SNPs whose allele distinguishes the superior haplotype; 'fully diagnostic' SNPs are candidates for KASP / marker-assisted selection.",
    "Candidate accessions - the haplotype carried by each accession in each block; ranked by the number of superior haplotypes carried (green). Good parents stack many superior haplotypes.",
    "Pairwise tests - all pairwise comparisons between haplotypes, with adjusted p-values.",
    "GWAS SNPs - the significant GWAS SNPs used to build the blocks.",
    "Parameters - the settings used for this analysis.",
    "",
    "HOW TO READ",
    "Haplotype H1 is the most frequent allele combination in a block; 'Other' pools rare combinations and accessions with too many missing calls.",
    "Letter groups: haplotypes sharing a letter do not differ significantly; 'a' is the best group.",
    "Evidence: 'Strong' = the superior haplotype is significantly better than all other haplotypes; 'Moderate' = better than some; 'None' = no significant effect.",
    "Effect vs population (%) = (haplotype mean - population mean) / population mean x 100."
  )
  openxlsx::addWorksheet(wb, "Read me", gridLines = FALSE)
  openxlsx::writeData(wb, "Read me", readme, colNames = FALSE)
  openxlsx::addStyle(wb, "Read me", openxlsx::createStyle(fontSize = 15, textDecoration = "bold", fontColour = "#1F4E78"), rows = 1, cols = 1)
  openxlsx::addStyle(wb, "Read me", openxlsx::createStyle(textDecoration = "bold"), rows = c(7, 17), cols = 1)
  openxlsx::setColWidths(wb, "Read me", cols = 1, widths = 140)

  # Block summary ------------------------------------------------------------
  blocks <- s$blocks
  blocks$start_mb <- blocks$start / 1e6
  blocks$end_mb <- blocks$end / 1e6
  blocks$length_kb <- (blocks$end - blocks$start) / 1e3
  bs <- pretty_names(blocks, c(
    block = "Block", chr = "Chromosome", start_mb = "Start (Mb)", end_mb = "End (Mb)",
    length_kb = "Length (kb)", n_snps = "SNPs", lead_snps = "GWAS SNPs", n_haplotypes = "Haplotypes",
    n_accessions = "Phenotyped accessions", test = "Test", p_value = "P value", fdr = "FDR (across blocks)",
    strength = "Evidence", superior = "Superior haplotype", superior_mean = paste0("Superior mean ", unit),
    population_mean = paste0("Population mean ", unit), superior_vs_pop_pct = "Superior vs population (%)",
    best_worst_diff = "Best - worst haplotype", superior_n = "Accessions with superior haplotype",
    recommendation = "Recommendation"))
  add_sheet(wb, "Block summary", bs, widths = c(22, 11, 10, 10, 10, 7, 30, 11, 12, 13, 10, 11, 10, 11, 12, 12, 13, 12, 13, 90),
            number_cols = c("Start (Mb)", "End (Mb)", "Length (kb)", names(bs)[15:18]), digits = "0.00")
  openxlsx::addStyle(wb, "Block summary", openxlsx::createStyle(numFmt = "0.00E+00"), rows = 2:(nrow(bs) + 1),
                     cols = which(names(bs) %in% c("P value", "FDR (across blocks)")), gridExpand = TRUE, stack = TRUE)
  sig_rows <- which(blocks$strength != "None") + 1
  if (length(sig_rows)) {
    openxlsx::addStyle(wb, "Block summary", green, rows = sig_rows, cols = which(names(bs) %in% c("Evidence", "Superior haplotype")),
                       gridExpand = TRUE, stack = TRUE)
  }

  # Haplotype effects ----------------------------------------------------------
  e <- s$effects
  e$freq_pct <- 100 * e$freq
  e$superior <- ifelse(e$superior, "yes", "")
  e$alleles <- gsub("|", " ", e$alleles, fixed = TRUE)
  he <- pretty_names(e, c(
    block = "Block", hap = "Haplotype", superior = "Superior", letters = "Group", rank = "Rank",
    n = "Accessions (genotyped)", freq_pct = "Frequency (%)", n_pheno = "Accessions (phenotyped)",
    mean = "Mean", ci_low = "95% CI low", ci_high = "95% CI high", sd = "SD", se = "SE", median = "Median",
    min = "Min", max = "Max", diff_vs_pop = "Difference vs population", pct_vs_pop = "Effect vs population (%)",
    alleles = "Alleles (SNPs in position order)"))
  add_sheet(wb, "Haplotype effects", he, widths = c(22, 10, 9, 7, 6, 11, 10, 11, rep(10, 10), 60),
            number_cols = names(he)[7:18])
  rows <- which(e$superior == "yes") + 1
  if (length(rows)) openxlsx::addStyle(wb, "Haplotype effects", green, rows = rows, cols = 1:5, gridExpand = TRUE, stack = TRUE)
  rows <- which(e$hap == "Other") + 1
  if (length(rows)) openxlsx::addStyle(wb, "Haplotype effects", grey, rows = rows, cols = seq_len(ncol(he)), gridExpand = TRUE, stack = TRUE)
  pct_col <- which(names(he) == "Effect vs population (%)")
  openxlsx::conditionalFormatting(wb, "Haplotype effects", cols = pct_col, rows = 2:(nrow(he) + 1),
                                  style = c("#F8696B", "#FFFFFF", "#63BE7B"), type = "colourScale")

  # Haplotype alleles ---------------------------------------------------------
  h <- results$haplotypes
  hap_cols <- hap_levels(setdiff(h$hap, "Other"))
  alleles <- do.call(rbind, lapply(unique(h$block), function(bl) {
    info <- block_info(h, bl)
    hb <- h[h$block == bl & h$hap != "Other", ]
    al <- do.call(rbind, strsplit(hb$alleles, "|", fixed = TRUE))
    calls <- matrix("", nrow(info), length(hap_cols), dimnames = list(NULL, hap_cols))
    calls[, hb$hap] <- t(al)
    data.frame(Block = bl, Marker = info$marker, SNP = info$snp, `Position (bp)` = info$pos,
               `GWAS SNP` = ifelse(info$snp %in% attr(info, "leads"), "yes", ""),
               `Major allele` = info$major, `Minor allele` = info$minor, MAF = info$maf,
               calls, check.names = FALSE, stringsAsFactors = FALSE)
  }))
  add_sheet(wb, "Haplotype alleles", alleles, widths = c(22, 24, 16, 13, 9, 8, 8, 7, rep(6, length(hap_cols))))

  # Diagnostic markers ----------------------------------------------------------
  dm <- results$markers
  if (!is.null(dm) && nrow(dm)) {
    dm$fully_diagnostic <- ifelse(dm$fully_diagnostic, "yes", "")
    dm$is_gwas_snp <- ifelse(dm$is_gwas_snp, "yes", "")
  }
  dm <- pretty_names(dm, c(
    block = "Block", superior = "Superior haplotype", marker = "Marker", snp = "SNP", pos = "Position (bp)",
    favourable_allele = "Favourable allele", other_alleles = "Other allele(s)",
    distinguishes = "Distinguishes from N haplotypes", fully_diagnostic = "Fully diagnostic",
    is_gwas_snp = "GWAS SNP", r2_with_gwas_snp = "LD (r2) with GWAS SNP"))
  add_sheet(wb, "Diagnostic markers", dm, widths = c(22, 11, 24, 16, 13, 11, 11, 13, 10, 9, 12))
  if (nrow(dm)) {
    rows <- which(dm[["Fully diagnostic"]] == "yes") + 1
    if (length(rows)) openxlsx::addStyle(wb, "Diagnostic markers", green, rows = rows, cols = 3:6, gridExpand = TRUE, stack = TRUE)
  }

  # Candidate accessions ----------------------------------------------------------
  r <- results$ranking
  rank_df <- r
  names(rank_df)[1:5] <- c("Rank", "Accession", paste0(trait, if (nzchar(unit)) paste0(" (", unit, ")") else ""),
                           "Superior haplotypes carried", "Significant blocks")
  add_sheet(wb, "Candidate accessions", rank_df, widths = c(6, 16, 12, 12, 11, rep(14, ncol(r) - 5)),
            number_cols = names(rank_df)[3])
  for (bl in s$blocks$block[s$blocks$strength != "None"]) {
    col <- which(names(r) == bl)
    rows <- which(!is.na(r[[bl]]) & r[[bl]] == s$blocks$superior[s$blocks$block == bl]) + 1
    if (length(rows)) openxlsx::addStyle(wb, "Candidate accessions", green, rows = rows, cols = col, stack = TRUE)
  }

  # Pairwise tests, GWAS SNPs, parameters -----------------------------------------
  pw <- do.call(rbind, results$pairwise)
  if (!is.null(pw)) {
    pw <- pretty_names(pw, c(block = "Block", group1 = "Haplotype 1", group2 = "Haplotype 2", n1 = "n1", n2 = "n2",
                             mean1 = "Mean 1", mean2 = "Mean 2", diff = "Difference", statistic = "Statistic",
                             df = "df", p = "P value", p.adj = "Adjusted P", p.adj.signif = "Significance"))
  }
  add_sheet(wb, "Pairwise tests", pw, widths = c(22, 11, 11, 6, 6, 10, 10, 10, 10, 8, 11, 11, 11),
            number_cols = c("Mean 1", "Mean 2", "Difference", "Statistic", "df"))
  if (!is.null(pw)) {
    openxlsx::addStyle(wb, "Pairwise tests", openxlsx::createStyle(numFmt = "0.00E+00"), rows = 2:(nrow(pw) + 1),
                       cols = 11:12, gridExpand = TRUE, stack = TRUE)
  }
  gw <- results$gwas_significant
  if (!is.null(gw) && !is.data.frame(gw)) gw <- do.call(rbind, gw)
  add_sheet(wb, "GWAS SNPs", gw)
  add_sheet(wb, "Parameters", list_config_df())

  openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
  cli::cli_alert_success("Excel report written to {.file {file}}")
  invisible(file)
}
