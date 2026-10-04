# Colours shared by all HaploTraitR plots
.hc <- list(
  superior = "#1A9850", hap = "#4575B4", other = "#B0B0B0",
  major = "#C6DBEF", minor = "#FD8D3C", het = "#9E9AC8", missing = "#FFFFFF",
  lead = "#D7301F", ld = c("#FFFFFF", "#FEE0D2", "#FC9272", "#EF3B2C", "#99000D")
)
.role_colours <- c(Superior = .hc$superior, Haplotype = .hc$hap, Other = .hc$other)

clamp <- function(x, lo, hi) max(lo, min(hi, x))

theme_haplo <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(colour = "grey30"),
          plot.caption = element_text(colour = "grey40", size = rel(0.75)),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "white", colour = NA))
}

# Group statistics used by the haplotype trait plot
group_stats <- function(d, pw) {
  alpha <- get_config("t_test_threshold")
  higher <- isTRUE(get_config("higher_is_better"))
  lv <- hap_levels(d$hap)
  st <- do.call(rbind, lapply(lv, function(g) {
    v <- d$Pheno[d$hap == g]
    data.frame(hap = g, n = length(v), mean = mean(v), max = max(v), stringsAsFactors = FALSE)
  }))
  testable <- st$hap[st$n >= 2]
  ranked <- testable[order(st$mean[match(testable, st$hap)] * if (higher) -1 else 1)]
  sig <- if (is.null(pw)) NULL else as.matrix(pw[!is.na(pw$p.adj) & pw$p.adj < alpha, c("group1", "group2")])
  st$letters <- ""
  if (length(ranked)) st$letters[match(ranked, st$hap)] <- cld_letters(ranked, sig)
  superior <- setdiff(ranked, "Other")[1]
  # only highlight the best haplotype when it differs from another haplotype
  if (!is.na(superior)) {
    sup_letters <- strsplit(st$letters[st$hap == superior], "")[[1]]
    others <- st$letters[st$hap %in% setdiff(ranked, c(superior, "Other"))]
    differs <- vapply(strsplit(others, ""), function(l) !any(l %in% sup_letters), logical(1))
    if (!any(differs)) superior <- NA
  }
  st$role <- ifelse(st$hap == "Other", "Other", ifelse(st$hap %in% superior, "Superior", "Haplotype"))
  st
}

#' Plot the trait distribution of each haplotype
#'
#' Violin + box plot of the trait for every haplotype of a block, with the
#' number of accessions, the mean (diamond), the population mean (dashed line)
#' and letter groups. The superior haplotype is shown in green.
#' @param cls_snp A block id or GWAS SNP. If `NULL`, the first block is plotted.
#' @param SNPcombTables The accession table (from [getSNPcombTables()])
#' @param t_test_snpComp Pairwise tests (from [testSNPcombs()])
#' @param with_dots Show the individual accessions
#' @param outfolder Folder in which to save the plot (`NULL`: do not save)
#' @param show_pairs Draw brackets for significant pairs (`NULL`: only when
#'   there are at most four haplotypes)
#' @return A ggplot object.
#' @export
plotHapCombBoxPlot <- function(cls_snp = NULL, SNPcombTables, t_test_snpComp, with_dots = TRUE,
                               outfolder = NULL, show_pairs = NULL) {
  blocks <- unique(SNPcombTables$block)
  if (is.null(cls_snp)) cls_snp <- blocks[1]
  block <- if (cls_snp %in% blocks) cls_snp else
    unique(SNPcombTables$block[vapply(strsplit(SNPcombTables$lead_snps, ";"), function(l) cls_snp %in% l, logical(1))])[1]
  if (is.na(block)) stop("'", cls_snp, "' is not a block or GWAS SNP in SNPcombTables.")

  d <- SNPcombTables[SNPcombTables$block == block, ]
  pw <- t_test_snpComp[[block]]
  st <- group_stats(d, pw)
  lv <- st$hap
  d$hap <- factor(d$hap, lv)
  d$role <- st$role[match(as.character(d$hap), st$hap)]
  pop <- mean(d$Pheno)
  rng <- diff(range(d$Pheno))
  if (rng == 0) rng <- 1

  p_block <- NA
  testable <- st$hap[st$n >= 2]
  if (length(testable) >= 2) {
    dt <- d[d$hap %in% testable, ]
    p_block <- tryCatch(
      if (get_config("t_test_method") == "wilcox.test") stats::kruskal.test(Pheno ~ droplevels(hap), dt)$p.value
      else stats::oneway.test(Pheno ~ droplevels(hap), dt)$p.value, error = function(e) NA)
  }
  test_name <- if (get_config("t_test_method") == "wilcox.test") "Kruskal-Wallis" else "Welch ANOVA"
  leads <- gsub(";", ", ", d$lead_snps[1])

  p <- ggplot(d, aes(x = hap, y = Pheno)) +
    geom_hline(yintercept = pop, linetype = "dashed", colour = "grey45") +
    geom_violin(aes(fill = role), alpha = 0.3, colour = NA, scale = "width", width = 0.85)
  if (with_dots) {
    p <- p + geom_jitter(aes(colour = role), width = 0.12, height = 0, size = 1.3, alpha = 0.55)
  }
  p <- p +
    geom_boxplot(aes(colour = role), fill = NA, width = 0.22, outlier.shape = NA, linewidth = 0.6) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white", colour = "black") +
    geom_text(data = st, aes(x = hap, y = max + 0.07 * rng, label = letters),
              fontface = "bold", size = 4.5, inherit.aes = FALSE)

  if (is.null(show_pairs)) show_pairs <- length(testable) <= 4
  if (show_pairs && !is.null(pw)) {
    sig <- pw[!is.na(pw$p.adj) & pw$p.adj < get_config("t_test_threshold"), ]
    if (nrow(sig)) {
      base <- max(d$Pheno) + 0.18 * rng
      br <- data.frame(xmin = match(sig$group1, lv), xmax = match(sig$group2, lv),
                       y = base + (seq_len(nrow(sig)) - 1) * 0.1 * rng, label = sig$p.adj.signif)
      p <- p +
        geom_segment(data = br, aes(x = xmin, xend = xmax, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.4) +
        geom_segment(data = br, aes(x = xmin, xend = xmin, y = y, yend = y - 0.025 * rng), inherit.aes = FALSE, linewidth = 0.4) +
        geom_segment(data = br, aes(x = xmax, xend = xmax, y = y, yend = y - 0.025 * rng), inherit.aes = FALSE, linewidth = 0.4) +
        geom_text(data = br, aes(x = (xmin + xmax) / 2, y = y + 0.025 * rng, label = label),
                  inherit.aes = FALSE, size = 3.5)
    }
  }

  p <- p +
    scale_fill_manual(values = .role_colours) +
    scale_colour_manual(values = .role_colours) +
    scale_x_discrete(labels = sprintf("%s\nn = %d", st$hap, st$n)) +
    labs(title = paste0(get_config("phenotypename"), " by haplotype - block ", block),
         subtitle = sprintf("%s p = %s  |  GWAS SNP(s): %s", test_name, format_p(p_block), leads),
         x = NULL, y = trait_label(),
         caption = paste("Green = superior haplotype; diamond = mean; dashed line = population mean.",
                         "Haplotypes sharing a letter do not differ significantly.", sep = "\n")) +
    theme_haplo() +
    theme(legend.position = "none", panel.grid.major.x = element_blank())

  if (!is.null(outfolder) && nzchar(outfolder)) {
    save_plot(p, file.path(outfolder, paste0(safe_name(block), "_trait")),
              width = clamp(2.5 + 1.1 * length(lv), 5, 12), height = 6)
  }
  p
}

#' Generate and save the haplotype trait plots of all blocks
#'
#' Plots are saved in `plots/haplotype_trait_plots/significant` or
#' `.../not_significant` according to the overall test of the block.
#' @param SNPcombTables The accession table (from [getSNPcombTables()])
#' @param t_test_snpComp Pairwise tests (from [testSNPcombs()])
#' @param pwidth,pheight Not used anymore; the size adapts to the number of haplotypes.
#' @return Invisibly, a named list of ggplot objects.
#' @export
generateHapCombBoxPlots <- function(SNPcombTables, t_test_snpComp, pwidth = NULL, pheight = NULL) {
  alpha <- get_config("t_test_threshold")
  plots <- list()
  for (block in unique(SNPcombTables$block)) {
    d <- SNPcombTables[SNPcombTables$block == block, ]
    groups <- names(which(table(d$hap) >= 2))
    p_block <- NA
    if (length(groups) >= 2) {
      p_block <- tryCatch(stats::oneway.test(Pheno ~ factor(hap), d[d$hap %in% groups, ])$p.value,
                          error = function(e) NA)
    }
    folder <- out_subdir("plots", "haplotype_trait_plots",
                         if (!is.na(p_block) && p_block < alpha) "significant" else "not_significant")
    plots[[block]] <- plotHapCombBoxPlot(block, SNPcombTables, t_test_snpComp, outfolder = folder)
  }
  invisible(plots)
}

#' Generate boxplots of the trait for each significant GWAS SNP
#'
#' @param pheno A data frame containing phenotype data (accession, trait)
#' @param gwas A list of significant GWAS SNPs (from [filter_gwas_data()])
#' @param hapmap A list of HapMap data frames (from [readHapmap()])
#' @param plot_width The width of the saved plots (inches)
#' @param plot_height The height of the saved plots (inches)
#' @param plot_all If `TRUE` plot all SNPs, otherwise only SNPs with a
#'   significant genotype effect
#' @return A named list of ggplot objects.
#' @export
boxplot_genotype_phenotype <- function(pheno, gwas, hapmap, plot_width = 5, plot_height = 5.5, plot_all = TRUE) {
  alpha <- get_config("t_test_threshold")
  outfolder <- get_outfolder(required = FALSE)
  if (!is.null(outfolder)) folder <- out_subdir("plots", "gwas_snp_plots")
  data <- get_pheno_geno(extract_hapmap(hapmap, gwas), pheno)
  gw <- do.call(rbind, gwas)
  plots <- list()
  for (snp in unique(data$variable)) {
    d <- data[data$variable == snp & data$value != "NN", ]
    if (nrow(d) == 0) next
    counts <- table(d$value)
    groups <- names(counts)[counts >= 2]
    p_snp <- NA
    if (length(groups) >= 2) {
      p_snp <- tryCatch(stats::oneway.test(Phenotype ~ factor(value), d[d$value %in% groups, ])$p.value,
                        error = function(e) NA)
    }
    if (!plot_all && (is.na(p_snp) || p_snp >= alpha)) next

    alleles <- allele_summary(matrix(d$value, nrow = 1))
    type <- ifelse(substr(d$value, 1, 1) != substr(d$value, 2, 2), "Heterozygous",
                   ifelse(substr(d$value, 1, 1) == alleles$major, "Major allele", "Minor allele"))
    d$type <- type
    lv <- names(sort(tapply(d$Phenotype, d$value, mean)))
    d$value <- factor(d$value, lv)
    marker <- gw$rsid[match(snp, gw$rs)]
    title_snp <- if (!is.na(marker) && marker != snp) paste0(marker, " (", snp, ")") else snp

    p <- ggplot(d, aes(x = value, y = Phenotype)) +
      geom_violin(aes(fill = type), alpha = 0.3, colour = NA, scale = "width", width = 0.8) +
      geom_jitter(aes(colour = type), width = 0.12, height = 0, size = 1.2, alpha = 0.55) +
      geom_boxplot(aes(colour = type), fill = NA, width = 0.22, outlier.shape = NA, linewidth = 0.6) +
      stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
      scale_fill_manual(values = c("Major allele" = .hc$hap, "Minor allele" = .hc$minor, Heterozygous = .hc$het)) +
      scale_colour_manual(values = c("Major allele" = .hc$hap, "Minor allele" = .hc$minor, Heterozygous = .hc$het)) +
      scale_x_discrete(labels = sprintf("%s\nn = %d", lv, as.integer(counts[lv]))) +
      labs(title = title_snp, subtitle = paste0("Welch test p = ", format_p(p_snp)),
           x = "Genotype", y = trait_label(), fill = NULL, colour = NULL) +
      theme_haplo() +
      theme(legend.position = "bottom", panel.grid.major.x = element_blank())
    if (!is.null(outfolder)) save_plot(p, file.path(folder, paste0(safe_name(snp), "_genotypes")), plot_width, plot_height)
    plots[[snp]] <- p
  }
  plots
}
