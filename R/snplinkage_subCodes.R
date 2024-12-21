# This file is part of the snplinkage package.
# The author of this software is Thomas Charlon.
# This software is under the MIT license.
# Copyright 2021 Thomas Charlon.
# See the LICENSE file at the root of the project "snplinkage" for details.
# Link: https://github.com/ThomasChln

#' Ggplot linkage disequilibrium
#'
#' Display SNP r2 correlations using points or diamonds with text.
#'
#' @param df_ld       Data frame with columns SNP_A, SNP_B, and R2.
#'                    As returned by the snprelate_ld function.
#' @param diamonds    Should the values be displayed as diamonds or points ?
#'                    Default is TRUE for less than 40 SNPs.
#' @param reverse     Reverse the display (horizontal symmetry)
#' @param reindex     If FALSE, SNPs are positionned following their IDs
#' @param point_size  Size for geom_point. Ignored if diamonds is TRUE.
#' @return ggplot
#' @export
ggplot_ld <- function(df_ld, diamonds = length(unique(df_ld$SNP_A)) < 40,
  point_size = 120 / sqrt(nrow(df_ld)),
  reverse = FALSE, reindex = TRUE) {

#  check_fx_args(df_ld = '!d+', diamonds = '!B1')
  required_names <- c('SNP_A', 'SNP_B', 'R2')
  names_idxs <- match(required_names, names(df_ld))
  missing <- is.na(names_idxs)
  # die_if(any(missing), paste('Df_ld must have columns:',
  #     paste0(required_names[missing], collapse = ', ')))
  df_ld <- as.data.frame(lapply(df_ld[required_names], as.numeric))

  if (reindex) {
    uniq_ids <- unique(unlist(df_ld[1:2]))
    df_ld[1:2] <- lapply(df_ld[1:2],
      function(ids) as.numeric(factor(ids, uniq_ids)))
  }
  args <- list(df_ld, reverse)
  fun <- paste0('.ggplot_ld_', if (diamonds) {
      'diamonds'
    } else {
      args[[3]] <- point_size
      'points'
    })
  ggplt <- do.call(fun, args)

  ggplt + theme(panel.background = element_blank(),
    panel.grid = element_blank(), axis.text = element_blank(),
    axis.title = element_blank(), axis.ticks = element_blank())
}

diamond_grid <- function(df_grid_snp, r2, window) {
  grid_snp = df_grid_snp
  grid_snp <- grid_snp[apply(grid_snp, 1, diff) <= window, ]
  grid_snp[, 1] <- grid_snp[, 1] - grid_snp[, 2]
  grid_snp[, 2] <- grid_snp[, 2] + .5 * grid_snp[, 1]
  grid_snp[, 1] <- grid_snp[, 1] - .5 * grid_snp[, 1]

  df_ld <- cbind.data.frame(grid_snp, r2)[c(2, 1, 3)]

  stats::setNames(df_ld, c('x', 'y', 'R2'))
}

#' @import ggplot2
.ggplot_ld_diamonds <- function(df_ld, reverse) {
  all_snps <- unlist(df_ld[1:2])
  seq_snp <- unique(all_snps)
  n_snp <- max(seq_snp)
  window <- max(apply(df_ld[1:2], 1, diff))
  df_ld = diamond_grid(df_ld[1:2], df_ld$R2, window)

  # Plot LD values
  inv_r2 <- pmin(1, pmax(0, 1 - df_ld$R2))
  df_ld$color <- grDevices::rgb(1, inv_r2, inv_r2)
  ggplt <- ggplot(df_ld, aes(.data$x, .data$y)) +
    diamond_annots(df_ld) +
    ggplot2::geom_text(label = round(df_ld$R2 * 100), size = 2.5)

  # Plot SNPs indexes
  df_snp <- data.frame(x = seq_snp, y = 0, color = 'gray')
  ggplt <- ggplt + diamond_annots(df_snp) +
    ggplot2::annotate('text', df_snp$x, df_snp$y,
      label = df_snp$x, size = 2.5, fontface = 'bold')

  ywindow <- c(- window / 2 - 0.5, 0.5)
  if (reverse) ywindow %<>% rev

  ggplt +
    ggplot2::coord_cartesian(xlim = c(0.5, n_snp + 0.5), ylim = ywindow, expand = FALSE)
}
#' @import ggplot2
.ggplot_ld_points <- function(df_ld, reverse, size) {

  # add snp indexes
  uniq_snps <- unique(unlist(df_ld[1:2]))
  df_snps <- data.frame(SNP_A = uniq_snps, SNP_B = uniq_snps)
  df_snps <- rbind(df_ld[if (reverse) 1:2 else 2:1], df_snps)

  # rotate points
  m_ld <- matrix(unlist(df_snps), 2, byrow = TRUE)
  theta <- -pi / 4
  rot_mat <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2)
  m_ld <- rot_mat %*% m_ld

  # range values
  ranges <- apply(m_ld, 1, range)
  df_snps <- lapply(list(x = 1, y = 2),
    function(idx) ((m_ld[idx, ] - ranges[1, idx]) / diff(ranges[, idx]) * length(uniq_snps)) + 0.5)

  # separate indexes and correlations
  seq_rows <- seq_len(nrow(df_ld))
  df_ld[1:2] <- lapply(df_snps, '[', seq_rows)
  df_snps <- lapply(df_snps, '[', -seq_rows)

  map <- aes(.data$SNP_A, .data$SNP_B, color = .data$R2)
  xframe <- c(0, length(uniq_snps) + 1)
  yframe <- if (reverse) rev(xframe) else xframe

  ggplot(df_ld, map) + geom_point(shape = 18, size = size) +
    ggplot2::annotate('point', x = df_snps$x, y = df_snps$y, shape = 18,
      color = grDevices::grey(.6), size = size * 0.9) +
    ggplot2::annotate('segment', x = xframe + c(-0.3, 0.3), xend = xframe[2] / 2,
      y = yframe[2], yend = yframe[1] + 0.6, lineend = 'round') +
    scale_colour_gradient(limits = 0:1, breaks = c(0, 0.5, 1), low = 'white',
      high = '#972D15', na.value = 'grey',
      guide = guide_colorbar(barheight = 5)) +
    ggplot2::coord_cartesian(xlim = xframe, ylim = xframe + c(0, 0.6), expand = FALSE) +
    labs(color = 'R2') +
    theme(legend.justification = c(1, reverse), legend.position = c(1, reverse))
}


#' Get diamond ggplot layer.
#'
#' Diamond ggplot layer for ggplot_ld
#'
#' @param data Data frame of 3 columns defining the diamonds
#' @param x Name of the column for horizontal positions
#' @param y Name of the column for vertical positions
#' @param color Name of the column for color values
#' @param size Radius of the diamonds
#' @return gglayers
#' @export
diamond_annots <- function(data, x = 'x', y = 'y', color = 'color', size = .5) {

#  check_fx_args(x = '!C1', y = '!C1', color = '!C1')

  # die_unless(all(c(x, y, color) %in% names(data)),
  #   paste('Data does not have columns',
  #     paste(c(x, y, color), collapse = ',')))

  data[c(x, y)] <- lapply(data[c(x, y)], as.numeric)
  bounds <- apply(data[c(x, y)], 1, diamond_bounds, x, y, size)

  color <- data[[color]]
  diamonds <- lapply(seq_along(bounds), function(bounds_idx) {
      bounds <- bounds[[bounds_idx]]
      ggplot2::annotate('polygon', bounds$X, bounds$Y, color = 'black',
        fill = color[bounds_idx])
    })

  diamonds
}


diamond_bounds <- function(diamond_coord, x, y, size) {
  bounds <- sapply(diamond_coord[c(x, y)],
    function(dim) sapply(c('-', '+'), Reduce, c(dim, size)))
  bounds <- data.frame('X' = c(bounds[, 1], rep(diamond_coord[[x]], 2)),
    'Y' = c(rep(diamond_coord[[y]], 2), bounds[, 2]))
  bounds[2:3, ] <- bounds[3:2, ]
  bounds
}


snp_position_colors <- function(n_snps,
  colors = c('#D8B70A', '#02401B', '#A2A475', '#972D15', '#81A88D')) {
  colors[seq_len(n_snps) %% length(colors) + 1]
}

.format_eng_range <- function(nums) {
  bases <- floor(log(nums, 1e3))
  base_strs <- bases
  base_strs[base_strs > 2] <- 2
  base_strs <- list(NULL, 'k', 'M')[base_strs + 1]
  nums <- nums / 1e3 ^ bases

  # keep 3 digits
  digit_base <- 10 ^ (3 - nchar(floor(nums)))
  nums <- nums * digit_base
  nums <- c(floor(nums[1]), ceiling(nums[2])) / digit_base
  nums <- paste(nums, base_strs)
}


#' Function to plot SNP positions
#' @param df_snp Data frame with columns position, labels_colname
#' @param upper_subset Subset of SNPs to highlight
#' @param labels_colname Column name for SNP labels
#' @param colors Colors for SNP positions
#' @importFrom ggplot2 geom_segment
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 coord_cartesian
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 annotate
#' @importFrom ggplot2 element_blank
#' @return ggplot
#' @export
plot_snp_positions <- function(df_snp, upper_subset = NULL, labels_colname = NULL,
                               colors = snp_position_colors(nrow(df_snp))) {
    # check_fx_args(df_snp = "!d+")
    required_names <- c("position", labels_colname)

    n_snp <- nrow(df_snp)
    range_pos <- range(df_snp$position)
    pos_size <- 0.15 * n_snp
    xrange <- c(1 + pos_size, n_snp - pos_size)
    df_snp$position <- (df_snp$position - range_pos[1]) * diff(xrange)
    df_snp$position <- df_snp$position / diff(range_pos) + xrange[1]
    df_snp$index <- seq(xrange[1], xrange[2], length = n_snp)
    df_snp$seq <- seq_len(n_snp)
    map <- aes(.data$position, xend = .data$index)
    yrange <- c(-2, 1.5)
    xrange[1:2] <- c(xrange[1] - 0.25, xrange[2] + 0.25)
    plot <- ggplot(df_snp, map) +
        geom_segment(aes(xend = .data$position), y = 0, yend = 1, color = colors) +
        geom_segment(aes(xend = .data$seq), y = -0.025, yend = yrange[1], color = colors, lineend = "round")
    annots_params <- list("text", y = 0.5)
    xpos <- c(xrange[1] / 2, xrange[2] + (n_snp - xrange[2]) / 2)
    range_pos <- paste0(.format_eng_range(range_pos), "bp")
    annots <- lapply(1:2, function(idx) {
        annots_params[c("x", "label")] <- list(xpos[idx], range_pos[idx])
        do.call(ggplot2::annotate, annots_params)
    })
    plot <- plot + annots
    if (!is.null(labels_colname)) {
        yrange[1] <- -7
        plot <- plot + list(
            geom_text(x = df_snp$seq - 0.25, y = -5, color = colors, label = df_snp[[labels_colname]], size = 2.5, angle = 90),
            geom_segment(x = df_snp$seq, xend = df_snp$seq, y = yrange[1], yend = -2, color = colors, lineend = "round")
        )
    }
    if (!is.null(upper_subset)) {
        colors <- colors[upper_subset]
        df_snp_subset <- df_snp[upper_subset, ]
        yrange[2] <- 3
        plot <- plot + geom_segment(aes(xend = .data$seq), df_snp_subset, y = 1.025, yend = yrange[2] - 0.05, color = colors, lineend = "round")
    }
    plot + ggplot2::annotate("segment", x = rep(xrange[1], 2), xend = rep(xrange[2], 2), y = 0:1, yend = 0:1, lineend = "square") +
        ggplot2::annotate("segment", x = xrange[1:2], xend = xrange[1:2], y = rep(0, 2), yend = rep(1, 2)) +
      ggplot2::coord_cartesian(xlim = c(0.5, n_snp + 0.5), ylim = c(yrange[1], yrange[2]), expand = FALSE) +
        theme(panel.background = element_blank(), text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())
}

#' Function to create LD plot and integrate with SNP positions
#' @param df_ld Data frame with columns SNP_A, SNP_B, and R2
#' @param df_snp Data frame with columns snpID, position, and labels_colname
#' @param biplot_subset Subset of SNPs to highlight in the biplot
#' @param labels_colname Column name for SNP labels
#' @param diamonds Should the values be displayed as diamonds or points ?
#' @param point_size Size for geom_point
#' @param title Title for the plot
#' @param title_biplot Title for the biplot
#' @import ggplot2
#' @import ggpubr
#' @import dplyr
#' @importFrom gtable gtable_add_rows
#' @importFrom gtable gtable_add_grob
#' @return ggplot
#' @export
create_ld_plot <- function(df_ld, df_snp, biplot_subset = NULL, labels_colname = NULL,
                           diamonds = length(unique(df_ld$SNP_A)) < 40, point_size = ifelse(is.null(biplot_subset), 120, 80) / sqrt(nrow(df_ld)),
                           title = "", title_biplot = "", ...) {
    unique_ids <- unique(unlist(df_ld[1:2]))
    df_snp <- df_snp[match(unique_ids, df_snp$snpID), ]
    df_ld[1:2] <- lapply(df_ld[1:2], function(ids) as.numeric(factor(ids, unique_ids)))
    plots <- list(
        snp_pos = plot_snp_positions(df_snp, biplot_subset, labels_colname),
        ld = ggplot_ld(df_ld, ..., point_size = point_size) + labs(title = title)
    )
    if (!is.null(biplot_subset)) {
        snp_subset <- with(df_ld, SNP_A %in% biplot_subset & SNP_B %in% biplot_subset)
        if (diamonds)
            df_ld <- df_ld[snp_subset, ]
        else df_ld$R2[!snp_subset] = NA
        plots$biplot <- ggplot_ld(df_ld, diamonds, ..., reverse = TRUE, reindex = FALSE, point_size = point_size) + labs(title = title_biplot)
    }
    combine_plots(plots, labels_colname, title)
}

#' Function to combine multiple ggplot objects into a single gtable object
#' @param plots List of ggplot objects
#' @param labels_colname Column name for SNP labels
#' @param title Title for the plot
#' @importFrom gtable gtable_add_rows
#' @importFrom gtable gtable_add_grob
#' @importFrom gtable gtable_filter
#' @importFrom ggplot2 ggplotGrob
#' @importFrom grid unit
#' @return gtable
#' @export
combine_plots <- function(plots, labels_colname, title) {
    plot_pos = if (getRversion() >= "4.3.3") {
        list(height = 9, width = 7, bottom = 17)
    }
    else {
        list(height = 7, width = 5, bottom = 13)
    }
    plots <- lapply(plots, ggplot2::ggplotGrob)
    is_biplot <- "biplot" %in% names(plots)
    ld_relative_size <- if (!is.null(labels_colname) || is_biplot)
        3
    else 7
    plots$snp_pos <- gtable::gtable_add_rows(plots$snp_pos, grid::unit(ld_relative_size, "null"), plot_pos$height)
    plots$snp_pos <- gtable::gtable_add_grob(plots$snp_pos, gtable::gtable_filter(plots$ld, "panel|guide-box-inside"), plot_pos$height + 1, plot_pos$width)
    if (title != "") {
        title_pos <- if (is_biplot)
            plot_pos$bottom
        else 2
        plots$snp_pos <- gtable::gtable_add_rows(plots$snp_pos, grid::unit(0.2, "null"), title_pos)
        plots$snp_pos <- gtable::gtable_add_grob(plots$snp_pos, gtable::gtable_filter(plots$ld, "title"), title_pos + 1, plot_pos$width)
    }
    if (is_biplot) {
        plots$snp_pos <- gtable::gtable_add_rows(plots$snp_pos, grid::unit(3, "null"), 2)
        plots$snp_pos <- gtable::gtable_add_grob(plots$snp_pos, gtable::gtable_filter(plots$biplot, "panel|guide-box|title"), 3, plot_pos$width)
    }
    plots$snp_pos
}
