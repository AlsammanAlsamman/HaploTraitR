# A figure made of vertically stacked ggplots that share the same x axis
new_figure <- function(plots, heights, width) {
  structure(list(plots = plots, heights = heights, width = width, height = sum(heights)),
            class = "haplotraitr_figure")
}

#' @export
print.haplotraitr_figure <- function(x, ...) {
  grid::grid.newpage()
  layout <- grid::grid.layout(length(x$plots), 1, heights = grid::unit(x$heights, "null"))
  grid::pushViewport(grid::viewport(layout = layout))
  for (i in seq_along(x$plots)) {
    grid::pushViewport(grid::viewport(layout.pos.row = i, layout.pos.col = 1))
    grid::grid.draw(ggplotGrob(x$plots[[i]]))
    grid::popViewport()
  }
  grid::popViewport()
  invisible(x)
}

# Allele class of each call relative to the major/minor allele of its SNP
allele_class <- function(call, major) {
  a1 <- substr(call, 1, 1)
  a2 <- substr(call, 2, 2)
  ifelse(call == "" | is.na(call), "Mixed",
         ifelse(a1 == "N", "Missing",
                ifelse(a1 != a2, "Heterozygous",
                       ifelse(a1 == major, "Major allele", "Minor allele"))))
}

# Rows of the haplotype matrix (best haplotype first, 'Other' last)
block_rows <- function(haplotypes, block, summary) {
  h <- haplotypes[haplotypes$block == block, ]
  e <- NULL
  if (!is.null(summary)) {
    e <- summary$effects[summary$effects$block == block, ]
    named <- e[e$hap != "Other", ]
    order <- c(named$hap[order(named$rank)], intersect("Other", e$hap))
    h <- h[match(intersect(order, h$hap), h$hap), ]
    e <- e[match(h$hap, e$hap), ]
  } else {
    h <- h[match(hap_levels(h$hap), h$hap), ]
  }
  list(h = h, e = e)
}

# Top panel: haplotype allele matrix (+ trait means when available)
block_matrix_panel <- function(haplotypes, block, summary, geom) {
  info <- block_info(haplotypes, block)
  rows <- block_rows(haplotypes, block, summary)
  h <- rows$h
  e <- rows$e
  n <- nrow(info)
  k <- nrow(h)
  y <- rev(seq_len(k))

  al <- lapply(h$alleles, function(a) if (nzchar(a)) strsplit(a, "|", fixed = TRUE)[[1]] else rep("", n))
  al <- do.call(rbind, al)
  tiles <- data.frame(x = rep(seq_len(n), each = k), y = rep(y, times = n),
                      call = as.vector(al), major = rep(info$major, each = k), stringsAsFactors = FALSE)
  tiles$fill <- factor(allele_class(tiles$call, tiles$major),
                       levels = c("Major allele", "Minor allele", "Heterozygous", "Missing", "Mixed"))
  tiles$label <- ifelse(tiles$call == "", "",
                        ifelse(substr(tiles$call, 1, 1) == substr(tiles$call, 2, 2),
                               substr(tiles$call, 1, 1), paste0(substr(tiles$call, 1, 1), "/", substr(tiles$call, 2, 2))))
  tiles$label[tiles$fill == "Missing"] <- "-"

  n_txt <- if (!is.null(h$n)) sprintf("n=%d, %.0f%%", h$n, 100 * h$freq) else sprintf("%.0f%%", 100 * h$freq)
  hap_names <- ifelse(h$hap == "Other", "Other (rare)", h$hap)
  left <- data.frame(x = 0.3, y = y, label = paste0(hap_names, "  ", n_txt))
  lead_x <- which(info$snp %in% attr(info, "leads"))

  p <- ggplot() +
    geom_tile(data = tiles, aes(x = x, y = y, fill = fill), colour = "white", linewidth = 0.4, height = 0.9) +
    scale_fill_manual(values = c("Major allele" = .hc$major, "Minor allele" = .hc$minor,
                                 Heterozygous = .hc$het, Missing = .hc$missing, Mixed = "#F2F2F2"),
                      drop = TRUE, name = NULL) +
    geom_text(data = left, aes(x = x, y = y, label = label), hjust = 1, size = 3.4) +
    annotate("text", x = 0.3, y = k + 0.95, label = "Haplotype (accessions)", hjust = 1,
             fontface = "bold", size = 3.3) +
    annotate("text", x = (n + 1) / 2, y = k + 0.95, label = "Alleles at the SNPs of the block",
             fontface = "bold", size = 3.3)
  if (n <= 90) {
    p <- p + geom_text(data = tiles, aes(x = x, y = y, label = label), size = geom$tile_text, colour = "grey15")
  }
  if (length(lead_x)) {
    p <- p + annotate("point", x = lead_x, y = k + 0.6, shape = 25, size = 2.8,
                      fill = .hc$lead, colour = .hc$lead)
  }

  if (!is.null(e)) {
    pop <- (e$mean - e$diff_vs_pop)[1]
    fx0 <- n + 1.5
    fx1 <- fx0 + geom$forest_width
    lo <- ifelse(is.na(e$ci_low), e$mean, e$ci_low)
    hi <- ifelse(is.na(e$ci_high), e$mean, e$ci_high)
    rng <- range(c(lo, hi, pop), na.rm = TRUE)
    if (diff(rng) == 0) rng <- rng + c(-1, 1)
    rng <- rng + c(-0.05, 0.05) * diff(rng)
    sx <- function(v) fx0 + (v - rng[1]) / diff(rng) * (fx1 - fx0)
    strong <- summary$blocks$strength[summary$blocks$block == block] != "None"
    role <- ifelse(e$hap == "Other", "Other", ifelse(e$superior & strong, "Superior", "Haplotype"))
    fd <- data.frame(y = y, x = sx(e$mean), lo = sx(lo), hi = sx(hi), role = role,
                     letters = e$letters, mean = sprintf("%.2f", e$mean), stringsAsFactors = FALSE)
    fd <- fd[!is.na(fd$x), ]
    ticks <- pretty(rng, 4)
    ticks <- ticks[ticks >= rng[1] & ticks <= rng[2]]
    p <- p +
      annotate("segment", x = sx(pop), xend = sx(pop), y = 0.5, yend = k + 0.5,
               linetype = "dashed", colour = "grey45") +
      geom_segment(data = fd, aes(x = lo, xend = hi, y = y, yend = y, colour = role), linewidth = 1) +
      geom_point(data = fd, aes(x = x, y = y), shape = 21, size = 3.4, colour = "black",
                 fill = .role_colours[fd$role]) +
      geom_text(data = fd, aes(x = fx1 + 0.6, y = y, label = letters), hjust = 0, fontface = "bold", size = 3.6) +
      geom_text(data = fd, aes(x = fx1 + 2.4, y = y, label = mean), hjust = 0, size = 3.3) +
      annotate("segment", x = fx0, xend = fx1, y = 0.12, yend = 0.12, colour = "grey30") +
      annotate("segment", x = sx(ticks), xend = sx(ticks), y = 0.12, yend = 0.0, colour = "grey30") +
      annotate("text", x = sx(ticks), y = -0.22, label = format(ticks), size = 2.8, colour = "grey20") +
      annotate("text", x = (fx0 + fx1) / 2, y = k + 0.95, label = paste0(trait_label(), ": mean \u00b1 95% CI"),
               fontface = "bold", size = 3.3) +
      annotate("text", x = fx1 + 0.6, y = k + 0.95, label = "Group", hjust = 0, fontface = "bold", size = 3.3) +
      annotate("text", x = fx1 + 2.4, y = k + 0.95, label = "Mean", hjust = 0, fontface = "bold", size = 3.3) +
      scale_colour_manual(values = .role_colours, guide = "none")
  }

  marker <- ifelse(is.na(info$marker) | info$marker %in% c("", "."), info$snp, info$marker)
  p + scale_x_continuous(breaks = seq_len(n), labels = marker) +
    coord_cartesian(xlim = c(0.5, geom$xmax), ylim = c(-0.45, k + 1.3), expand = FALSE, clip = "off") +
    theme_void(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = geom$label_pt, colour = "grey20"),
          legend.position = "top", legend.justification = "left",
          plot.title = element_text(face = "bold", size = 14, margin = margin(b = 4)),
          plot.subtitle = element_text(colour = "grey30", size = 10.5, margin = margin(b = 6)),
          plot.title.position = "plot",
          plot.background = element_rect(fill = "white", colour = NA),
          plot.margin = margin(0.15, geom$right_in, 0.05, geom$left_in, unit = "in"))
}

# LD triangle (diamond heat map) with a physical position map on top
ld_panel <- function(ld, pos, lead_x, geom, show_values) {
  n <- nrow(ld)
  idx <- which(upper.tri(ld), arr.ind = TRUE)
  i <- idx[, 1]
  j <- idx[, 2]
  cx <- (i + j) / 2
  cy <- -(j - i - 1) / 2
  r2 <- ld[idx]
  poly <- data.frame(id = rep(seq_along(cx), each = 4),
                     x = rep(cx, each = 4) + c(-0.5, 0, 0.5, 0),
                     y = rep(cy, each = 4) + c(0, 0.5, 0, -0.5),
                     R2 = rep(r2, each = 4))
  px <- if (diff(range(pos)) > 0) 1 + (pos - min(pos)) / diff(range(pos)) * (n - 1) else seq_len(n)
  is_lead <- seq_len(n) %in% lead_x
  map <- data.frame(x = px, xend = seq_len(n), y = 2.05, yend = 0.55, is_lead = is_lead)

  p <- ggplot() +
    geom_polygon(data = poly, aes(x = x, y = y, group = id, fill = R2),
                 colour = if (n <= 70) "white" else NA, linewidth = 0.2) +
    scale_fill_gradientn(colours = .hc$ld, limits = c(0, 1), na.value = "grey92",
                         name = expression(LD~r^2)) +
    annotate("rect", xmin = 1, xmax = max(n, 1.01), ymin = 2.05, ymax = 2.3, fill = "grey88", colour = "grey40") +
    geom_segment(data = map[!map$is_lead, ], aes(x = x, xend = xend, y = y, yend = yend),
                 colour = "grey55", linewidth = 0.3) +
    geom_segment(data = map[map$is_lead, ], aes(x = x, xend = xend, y = y, yend = yend),
                 colour = .hc$lead, linewidth = 0.7) +
    annotate("segment", x = px, xend = px, y = 2.05, yend = 2.3, colour = ifelse(is_lead, .hc$lead, "grey30"),
             linewidth = 0.3) +
    annotate("text", x = 1, y = 2.55, label = sprintf("%.3f Mb", min(pos) / 1e6), hjust = 0, size = 3.2) +
    annotate("text", x = n, y = 2.55, label = sprintf("%.3f Mb", max(pos) / 1e6), hjust = 1, size = 3.2) +
    annotate("text", x = 0.3, y = 2.17, label = "Physical position", hjust = 1, size = 3.2) +
    annotate("text", x = 0.3, y = -n / 8, label = "Linkage disequilibrium", hjust = 1, size = 3.2)
  if (show_values) {
    vals <- data.frame(x = cx, y = cy, label = ifelse(is.na(r2), "", round(100 * r2)))
    p <- p + geom_text(data = vals, aes(x = x, y = y, label = label), size = geom$ld_text, colour = "grey15")
  }
  p + coord_cartesian(xlim = c(0.5, geom$xmax), ylim = c(-(n - 2) / 2 - 0.6, 2.8),
                      expand = FALSE, clip = "off") +
    theme_void(base_size = 11) +
    theme(legend.position = "inside", legend.position.inside = geom$ld_legend,
          legend.justification = c(1, 0.5),
          plot.background = element_rect(fill = "white", colour = NA),
          plot.margin = margin(0.05, geom$right_in, 0.15, geom$left_in, unit = "in"))
}

# Sizes of a block figure, so the matrix columns and the LD diamonds line up
block_geometry <- function(n, k, has_effects, label_chars) {
  forest_width <- if (has_effects) max(6, 0.25 * n) else 0
  xmax <- if (has_effects) n + 1.5 + forest_width + 4.5 else n + 0.5
  range_x <- xmax - 0.5
  unit_in <- clamp(14 / range_x, 0.1, 0.35)
  left_in <- 2.1
  right_in <- 0.3
  label_pt <- clamp(unit_in * 72 * 0.75, 4, 8.5)
  row_in <- clamp(unit_in, 0.3, 0.42)
  top_h <- 1.25 + (k + 1.75) * row_in + label_chars * 0.6 * label_pt / 72 + 0.15
  ld_h <- (n / 2 + 2.2) * unit_in + 0.2
  list(forest_width = forest_width, xmax = xmax, unit_in = unit_in, left_in = left_in,
       right_in = right_in, label_pt = label_pt, top_h = top_h, ld_h = ld_h,
       width = left_in + right_in + range_x * unit_in,
       tile_text = clamp(unit_in * 9, 1.4, 3.4), ld_text = clamp(unit_in * 7, 1.2, 2.6),
       ld_legend = c(if (has_effects) 1 - 4.5 / range_x else 0.98, 0.45))
}

#' Plot a haplotype block: alleles, trait effects and LD
#'
#' Produces a three-part figure for one haplotype block:
#' * top: the alleles of every haplotype at every SNP (major allele in blue,
#'   minor allele in orange; red triangles mark the GWAS SNPs);
#' * right: the trait mean and 95% confidence interval of each haplotype with
#'   letter groups (green = superior haplotype, dashed line = population mean);
#' * bottom: the physical position of the SNPs and their pairwise LD (r2).
#' @param block A block id or GWAS SNP of the block
#' @param haplotypes A haplotype table (from [getHapCombSamples()] or [convertLDclusters2Haps()])
#' @param summary Optional output of [summarize_haplotypes()]; adds the trait panel
#' @return A `haplotraitr_figure` object (print it to draw it). Its `width`
#'   and `height` elements give a suitable size in inches.
#' @export
plot_haplotype_block <- function(block, haplotypes, summary = NULL) {
  block <- resolve_block(haplotypes, block)
  info <- block_info(haplotypes, block)
  n <- nrow(info)
  k <- sum(haplotypes$block == block)
  marker <- ifelse(is.na(info$marker), info$snp, info$marker)
  geom <- block_geometry(n, k, !is.null(summary), max(nchar(marker)))

  top <- block_matrix_panel(haplotypes, block, summary, geom)
  title <- sprintf("Haplotype block %s  |  %d SNPs  |  %s", block, n, get_config("phenotypename"))
  leads <- gsub(";", ", ", haplotypes$lead_snps[haplotypes$block == block][1])
  subtitle <- paste0("GWAS SNP(s) (red triangles): ", leads)
  if (!is.null(summary)) {
    b <- summary$blocks[summary$blocks$block == block, ]
    if (nrow(b)) {
      subtitle <- paste0(subtitle, sprintf("\n%s p = %s  |  superior haplotype: %s (%+.1f%% vs population mean)  |  evidence: %s",
                                           b$test, format_p(b$p_value), b$superior, b$superior_vs_pop_pct, b$strength))
    }
  }
  top <- top + labs(title = title, subtitle = subtitle)
  ld <- attr(info, "ld")
  bottom <- ld_panel(ld, info$pos, which(info$snp %in% attr(info, "leads")), geom, show_values = n <= 40)
  fig <- new_figure(list(top, bottom), c(geom$top_h, geom$ld_h), geom$width)
  fig
}

#' Plot LD and haplotype matrices of all blocks
#'
#' Saves one [plot_haplotype_block()] figure per block in
#' `plots/haplotype_blocks`.
#' @param clusterLDs Not used; kept for compatibility with HaploTraitR < 0.2
#' @param haplotypes A haplotype table (from [getHapCombSamples()])
#' @param gwas Not used; kept for compatibility with HaploTraitR < 0.2
#' @param summary Optional output of [summarize_haplotypes()]
#' @return Invisibly, a named list of `haplotraitr_figure` objects.
#' @export
plotLDCombMatrix <- function(clusterLDs = NULL, haplotypes, gwas = NULL, summary = NULL) {
  folder <- out_subdir("plots", "haplotype_blocks")
  figs <- list()
  for (block in unique(haplotypes$block)) {
    fig <- plot_haplotype_block(block, haplotypes, summary)
    save_plot(fig, file.path(folder, paste0(safe_name(block), "_block")), fig$width, fig$height)
    figs[[block]] <- fig
  }
  invisible(figs)
}

#' Plot the haplotype allele matrix of a block
#'
#' @param cluster_id A block id or GWAS SNP
#' @param haplotypes A haplotype table (from [convertLDclusters2Haps()] or [getHapCombSamples()])
#' @param gwas Not used; kept for compatibility
#' @param snps Not used; kept for compatibility
#' @param outfolder Folder in which to save the plot (`NULL`: do not save)
#' @return A ggplot object.
#' @export
plotCombMatrix <- function(cluster_id, haplotypes, gwas = NULL, snps = NULL, outfolder = NULL) {
  block <- resolve_block(haplotypes, cluster_id)
  info <- block_info(haplotypes, block)
  k <- sum(haplotypes$block == block)
  geom <- block_geometry(nrow(info), k, FALSE, max(nchar(info$marker)))
  p <- block_matrix_panel(haplotypes, block, NULL, geom) + labs(title = paste("Haplotypes of block", block))
  if (!is.null(outfolder)) {
    save_plot(p, file.path(outfolder, paste0(safe_name(block), "_haplotypes")), geom$width, geom$top_h)
  }
  p
}

#' Plot an LD heat map
#'
#' @param LD_matrix A square LD (r2) matrix whose row names are `chr:pos`
#' @return A ggplot object.
#' @export
plotLDheatmap <- function(LD_matrix) {
  LD_matrix <- as.matrix(LD_matrix)
  pos <- as.numeric(sub("^.*:", "", rownames(LD_matrix)))
  ord <- order(pos)
  LD_matrix <- LD_matrix[ord, ord]
  n <- nrow(LD_matrix)
  geom <- block_geometry(n, 1, FALSE, 1)
  geom$left_in <- 1.4
  ld_panel(LD_matrix, pos[ord], integer(0), geom, show_values = n <= 40)
}

#' Plot LD heat maps for all LD blocks
#'
#' @param clusterLDs A list of LD blocks (from [getLDclusters()])
#' @param LDsInfo Optional output of [computeLDclusters()]; otherwise the
#'   matrices are read from the `LD_matrices` folder of the output folder
#' @return Invisibly, a list of ggplot objects (also saved in `plots/LD_heatmaps`).
#' @export
plotLDForClusters <- function(clusterLDs, LDsInfo = NULL) {
  if (is.null(LDsInfo)) LDsInfo <- retrieveLDMatricesFromFolder(file.path(get_outfolder(), "LD_matrices"))
  folder <- out_subdir("plots", "LD_heatmaps")
  plots <- list()
  for (cls in names(clusterLDs)) {
    snps <- clusterLDs[[cls]]
    if (length(snps) == 0 || is.null(LDsInfo$matrices[[cls]])) next
    p <- plotLDheatmap(LDsInfo$matrices[[cls]][snps, snps])
    n <- length(snps)
    geom <- block_geometry(n, 1, FALSE, 1)
    save_plot(p, file.path(folder, paste0(safe_name(cls), "_LD")), geom$width, geom$ld_h + 0.3)
    plots[[cls]] <- p
  }
  invisible(plots)
}
