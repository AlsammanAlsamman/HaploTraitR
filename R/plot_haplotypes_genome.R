#' Plot Haplotypes SNPs Across the Genome
#'
#' @param gwas_data A data frame or list of data frames containing the GWAS results.
#' @param haplotype_data A data frame containing the haplotype information.
#' @param use_facet Logical, whether to use facet wrapping for chromosomes.
#' @return A list of ggplot2 objects representing the plots.
#' @import ggplot2
#' @export
plot_haplotypes_genome <- function(gwas_data, haplotype_data, use_facet = TRUE, pwidth = 10, pheight = 5, pdpi = 300) {
  # Ensure output folder exists
  outfolder <- get_config("outfolder")
  if (!dir.exists(outfolder)) {
    stop("The output folder does not exist, please create it first.")
  }

  # Prepare haplotype data
  haplotype_data <- haplotype_data[, c("snps", "snp")]
  hap_df <- data.frame()
  for (i in 1:nrow(haplotype_data)) {
    snps <- strsplit(as.character(haplotype_data$snps[i]), "\\|")
    snp <- haplotype_data$snp[i]
    for (j in 1:length(snps[[1]])) {
      hap_df <- rbind(hap_df, data.frame(snps = snps[[1]][j], haplotype = snp))
    }
  }
  colnames(hap_df) <- c("snps", "haplotype")
  # remove snps with NA in haplotype
  hap_df <- hap_df[!is.na(hap_df$haplotype), ]

  hap_df$chr <- sapply(strsplit(as.character(hap_df$snps), ":"), function(x) x[1])


  # Combine GWAS data if it's a list
  if (!is.data.frame(gwas_data)) {
    gwas_data <- do.call(rbind, gwas_data)
  }

  combined_gwas <- gwas_data
  combined_gwas$chr <- as.factor(combined_gwas$chr)
  combined_gwas$pos <- as.numeric(combined_gwas$pos) / 1e6
  combined_gwas$p <- -log10(combined_gwas$p)
  combined_gwas$haplotype <- hap_df$haplotype[match(combined_gwas$rs, hap_df$snps)]

  # Filter combined_gwas to include only chromosomes with haplotypes
  combined_gwas <- combined_gwas[combined_gwas$chr %in% unique(hap_df$chr), ]

  # Plotting
  plots <- list()
  if (use_facet) {
    p <- ggplot() +
      geom_point(data = combined_gwas, aes(x = pos, y = p, color = haplotype), size = 3, alpha = 0.7) +
      labs(title = "Haplotypes SNPs Across the Genome",
           x = "Position (Mb)",
           y = "-log10(P-value)",
           color = "Haplotypes") +
      facet_wrap(~chr, scales = "free_x")
    plots[["All_Chromosomes"]] <- p
    ggsave(filename = file.path(outfolder, "Haplotypes_SNPs_Across_Genome.png"), plot = p, width = pwidth, height = pheight, dpi = pdpi)
  } else {
    chromosomes <- unique(combined_gwas$chr)
    for (chr in chromosomes) {
      chr_data <- combined_gwas[combined_gwas$chr == chr, ]
      p <- ggplot() +
        geom_point(data = chr_data, aes(x = pos, y = p, color = haplotype), size = 3, alpha = 0.7) +
        labs(title = paste("Haplotypes SNPs Across Chromosome", chr),
             x = "Position (Mb)",
             y = "-log10(P-value)",
             color = "Haplotypes")
      plots[[as.character(chr)]] <- p
      ggsave(filename = file.path(outfolder, paste0("Haplotypes_SNPs_Chromosome_", chr, ".png")), plot = p, width = pwidth, height = pheight, dpi = pdpi)
    }
  }

  return(plots)
}
