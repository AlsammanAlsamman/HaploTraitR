#' Plot the haplotype blocks across the genome
#'
#' Manhattan plot of the GWAS results in which the haplotype blocks are
#' highlighted (green: significant haplotype effect, orange: not significant)
#' and labelled.
#' @param gwas_data The full GWAS results (data frame from [read_gwas_file()],
#'   or a list of data frames)
#' @param haplotype_data A haplotype table (from [getHapCombSamples()])
#' @param use_facet If `TRUE` a genome-wide plot is drawn; if `FALSE` one
#'   zoomed plot per chromosome carrying a block
#' @param pwidth,pheight,pdpi Size (inches) and resolution of the saved plot(s)
#' @param summary Optional output of [summarize_haplotypes()] used to colour the blocks
#' @return A named list of ggplot objects (also saved in `plots/`).
#' @export
plot_haplotypes_genome <- function(gwas_data, haplotype_data, use_facet = TRUE, pwidth = 12,
                                   pheight = 5, pdpi = 300, summary = NULL) {
  old <- set_config(list(plot_dpi = pdpi))
  on.exit(set_config(old))
  if (!is.data.frame(gwas_data)) gwas_data <- do.call(rbind, gwas_data)
  g <- gwas_data[!is.na(gwas_data$p) & gwas_data$p > 0, ]
  g$chr <- as.character(g$chr)
  g$logp <- -log10(g$p)
  fdr <- p.adjust(g$p, method = get_config("fdr_method"))
  sig <- fdr < get_config("fdr_threshold")
  threshold <- if (any(sig)) -log10(max(g$p[sig])) else NA

  blocks <- unique(haplotype_data[, c("block", "chr", "start", "end")])
  blocks$status <- "Haplotype block"
  if (!is.null(summary)) {
    s <- summary$blocks[match(blocks$block, summary$blocks$block), ]
    blocks$status <- ifelse(!is.na(s$strength) & s$strength != "None",
                            "Significant haplotype effect", "No significant effect")
  }
  block_cols <- c("Significant haplotype effect" = .hc$superior, "No significant effect" = "#F28E2B",
                  "Haplotype block" = .hc$hap)
  leads <- unique(unlist(strsplit(haplotype_data$lead_snps, ";")))
  ymax <- max(g$logp) * 1.12
  folder <- out_subdir("plots")
  plots <- list()

  if (use_facet) {
    chrs <- chr_order(g$chr)
    chr_len <- tapply(g$pos, g$chr, max)[chrs]
    gap <- 0.01 * sum(chr_len)
    offset <- stats::setNames(c(0, cumsum(chr_len + gap)[-length(chrs)]), chrs)
    g$pos_cum <- g$pos + offset[g$chr]
    g$chr_col <- factor(match(g$chr, chrs) %% 2)
    centers <- offset + chr_len / 2
    blocks$xmin <- blocks$start + offset[blocks$chr]
    blocks$xmax <- blocks$end + offset[blocks$chr]
    min_w <- 0.003 * sum(chr_len)
    widen <- pmax(0, min_w - (blocks$xmax - blocks$xmin)) / 2
    blocks$xmin <- blocks$xmin - widen
    blocks$xmax <- blocks$xmax + widen
    lead_pts <- g[g$rs %in% leads, ]

    p <- ggplot(g, aes(x = pos_cum, y = logp)) +
      geom_rect(data = blocks, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = ymax, fill = status),
                inherit.aes = FALSE, alpha = 0.25) +
      geom_point(aes(colour = chr_col), size = 0.7, alpha = 0.7) +
      scale_colour_manual(values = c("0" = "grey35", "1" = "grey70"), guide = "none") +
      geom_point(data = lead_pts, colour = .hc$lead, size = 2) +
      geom_text(data = blocks, aes(x = (xmin + xmax) / 2, y = ymax, label = block),
                inherit.aes = FALSE, angle = 90, hjust = 1, vjust = -0.3, size = 2.8, check_overlap = TRUE) +
      scale_fill_manual(values = block_cols, name = NULL) +
      scale_x_continuous(breaks = centers, labels = chrs, expand = c(0.01, 0)) +
      scale_y_continuous(limits = c(0, ymax), expand = c(0, 0)) +
      labs(title = paste("GWAS signals and haplotype blocks -", get_config("phenotypename")),
           subtitle = "Red points: GWAS SNPs used to build the haplotype blocks",
           x = "Chromosome", y = expression(-log[10](p))) +
      theme_haplo() +
      theme(legend.position = "bottom", panel.grid.major.x = element_blank())
    if (!is.na(threshold)) {
      p <- p + geom_hline(yintercept = threshold, linetype = "dashed", colour = .hc$lead, linewidth = 0.4)
    }
    save_plot(p, file.path(folder, "genome_overview"), pwidth, pheight)
    plots[["All_Chromosomes"]] <- p
  } else {
    for (chr in chr_order(blocks$chr)) {
      gc <- g[g$chr == chr, ]
      bc <- blocks[blocks$chr == chr, ]
      p <- ggplot(gc, aes(x = pos / 1e6, y = logp)) +
        geom_rect(data = bc, aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 0, ymax = ymax, fill = status),
                  inherit.aes = FALSE, alpha = 0.3) +
        geom_point(size = 0.8, colour = "grey45", alpha = 0.7) +
        geom_point(data = gc[gc$rs %in% leads, ], colour = .hc$lead, size = 2) +
        geom_text(data = bc, aes(x = (start + end) / 2e6, y = ymax, label = block), inherit.aes = FALSE,
                  angle = 90, hjust = 1, vjust = -0.3, size = 3, check_overlap = TRUE) +
        scale_fill_manual(values = block_cols, name = NULL) +
        scale_y_continuous(limits = c(0, ymax), expand = c(0, 0)) +
        labs(title = paste("Chromosome", chr, "-", get_config("phenotypename")),
             x = "Position (Mb)", y = expression(-log[10](p))) +
        theme_haplo() + theme(legend.position = "bottom")
      if (!is.na(threshold)) p <- p + geom_hline(yintercept = threshold, linetype = "dashed", colour = .hc$lead)
      save_plot(p, file.path(folder, paste0("chromosome_", safe_name(chr))), pwidth, pheight)
      plots[[chr]] <- p
    }
  }
  plots
}

#' Plot the haplotype effects of all blocks
#'
#' One row per block; each point is a haplotype placed at its difference from
#' the population mean (in %), with its 95% confidence interval and a size
#' proportional to its frequency. The superior haplotype is green.
#' @param summary The list returned by [summarize_haplotypes()]
#' @return A ggplot object (also saved as `plots/haplotype_effects_summary`).
#' @export
plot_effect_summary <- function(summary) {
  e <- summary$effects[!is.na(summary$effects$mean), ]
  b <- summary$blocks
  significant <- b$block[b$strength != "None"]
  pop <- e$mean - e$diff_vs_pop
  e$lo <- 100 * (e$ci_low - pop) / pop
  e$hi <- 100 * (e$ci_high - pop) / pop
  e$role <- ifelse(e$hap == "Other", "Other",
                   ifelse(e$superior & e$block %in% significant, "Superior", "Haplotype"))
  e$freq_pct <- 100 * e$freq
  star <- ifelse(b$strength == "None", "", "  *")
  labels <- stats::setNames(paste0(b$block, star), b$block)
  order <- b$block[order(b$strength == "None", -abs(b$superior_vs_pop_pct))]
  e$block_label <- factor(labels[e$block], levels = rev(labels[order]))
  higher <- isTRUE(get_config("higher_is_better"))

  p <- ggplot(e, aes(x = pct_vs_pop, y = block_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45") +
    geom_errorbar(aes(xmin = lo, xmax = hi, colour = role), orientation = "y", width = 0.25, alpha = 0.6,
                   position = position_dodge(width = 0.6)) +
    geom_point(aes(size = freq_pct, fill = role, group = hap), shape = 21, colour = "black",
               position = position_dodge(width = 0.6)) +
    geom_text(aes(label = hap, group = hap), position = position_dodge(width = 0.6), vjust = -1.1, size = 3) +
    scale_fill_manual(values = .role_colours, name = NULL) +
    scale_colour_manual(values = .role_colours, guide = "none") +
    scale_size_area(max_size = 7, name = "Frequency (%)") +
    labs(title = paste("Haplotype effects on", get_config("phenotypename")),
         subtitle = paste0("Difference from the population mean (%) with 95% CI. ",
                           if (higher) "Higher" else "Lower", " values are favourable. * = significant block.\nGreen = superior haplotype; grey = rare haplotypes pooled as 'Other'."),
         x = "Difference from population mean (%)", y = NULL) +
    theme_haplo() +
    theme(legend.position = "bottom", panel.grid.major.y = element_line(colour = "grey92"))
  if (!is.null(get_outfolder(required = FALSE))) {
    save_plot(p, file.path(out_subdir("plots"), "haplotype_effects_summary"), 10, 2 + 0.7 * nrow(b))
  }
  p
}

#' Plot the best candidate accessions
#'
#' Heat map of the haplotypes carried by the top-ranked accessions (see
#' [rank_accessions()]) in the significant blocks; green cells mark the
#' superior haplotypes. Useful to choose parents that stack favourable haplotypes.
#' @param ranking The data frame returned by [rank_accessions()]
#' @param summary The list returned by [summarize_haplotypes()]
#' @param top Number of accessions to show
#' @return A ggplot object (also saved as `plots/candidate_accessions`).
#' @export
plot_candidate_accessions <- function(ranking, summary, top = 30) {
  b <- summary$blocks
  blocks <- b$block[b$strength != "None"]
  subtitle <- "Green: superior haplotype of a block with a significant effect"
  if (length(blocks) == 0) {
    blocks <- b$block
    subtitle <- "No block has a significant effect; superior haplotypes are shown for information only"
  }
  r <- head(ranking, top)
  long <- do.call(rbind, lapply(blocks, function(bl) {
    sup <- b$superior[b$block == bl]
    data.frame(accession = r$accession, block = bl, hap = r[[bl]],
               status = ifelse(is.na(r[[bl]]), "Not genotyped",
                               ifelse(r[[bl]] == sup, "Superior haplotype",
                                      ifelse(r[[bl]] == "Other", "Rare haplotype", "Other haplotype"))),
               stringsAsFactors = FALSE)
  }))
  long$accession <- factor(long$accession, levels = rev(r$accession))
  long$block <- factor(long$block, levels = blocks)
  trait <- data.frame(accession = factor(r$accession, levels = rev(r$accession)),
                      label = ifelse(is.na(r$trait), "NA", sprintf("%.2f", r$trait)),
                      score = sprintf("%d/%d", r$n_superior, r$n_significant_blocks))
  nb <- length(blocks)
  p <- ggplot(long, aes(x = block, y = accession)) +
    geom_tile(aes(fill = status), colour = "white", linewidth = 0.6) +
    geom_text(aes(label = hap), size = 3, colour = "grey15") +
    geom_text(data = trait, aes(x = nb + 0.9, y = accession, label = label), inherit.aes = FALSE,
              hjust = 0, size = 3.2) +
    geom_text(data = trait, aes(x = nb + 1.9, y = accession, label = score), inherit.aes = FALSE,
              hjust = 0, size = 3.2) +
    annotate("text", x = c(nb + 0.9, nb + 1.9), y = nrow(r) + 0.9, label = c(get_config("phenotypename"), "Superior"),
             hjust = 0, fontface = "bold", size = 3.3) +
    scale_fill_manual(values = c("Superior haplotype" = .hc$superior, "Other haplotype" = "#C6DBEF",
                                 "Rare haplotype" = "#EEEEEE", "Not genotyped" = "white"), name = NULL) +
    scale_x_discrete(expand = expansion(add = c(0.6, 2.6))) +
    coord_cartesian(clip = "off", ylim = c(0.5, nrow(r) + 0.5)) +
    labs(title = paste("Top", nrow(r), "candidate accessions -", get_config("phenotypename")),
         subtitle = subtitle, x = NULL, y = NULL) +
    theme_haplo() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid = element_blank(),
          legend.position = "bottom", plot.margin = margin(25, 10, 10, 10))
  if (!is.null(get_outfolder(required = FALSE))) {
    save_plot(p, file.path(out_subdir("plots"), "candidate_accessions"),
              clamp(3.5 + 1.3 * nb, 5.5, 20), clamp(2.5 + 0.24 * nrow(r), 5, 30))
  }
  p
}
