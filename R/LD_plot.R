#' Plotting LD Heatmap for a Specific Haplotype
#'
#' This function plots an LD heatmap for a specific haplotype using an LD matrix and hapmap object.
#' @param LD_matrix The LD matrix.
#' @return A ggplot object showing the LD heatmap.
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom gtable gtable
#' @importFrom gtable gtable_add_grob
#' @importFrom gtable gtable_add_cols
#' @importFrom gtable gtable_add_rows
#' @importFrom snplinkage gtable_ld
#' @export
plotLDheatmap <- function(LD_matrix) {
  if (is.null(LD_matrix) ) {
    stop("LD_matrix must be provided.")
  }

  snps <- rownames(LD_matrix)

  # Process SNP names to extract chromosome and position
  snp_info <- do.call(rbind, strsplit(snps, ":"))

  snp_info <- data.frame(
    chromosome = snp_info[, 1],
    position = as.integer(snp_info[, 2]),
    probe_id = snps,
    stringsAsFactors = FALSE,
    alleleA = NA,
    alleleB = NA,
    snpID = 1:nrow(snp_info)
  )
  # Sort SNP info by chromosome and position
  snp_info <- snp_info[order(snp_info$chromosome, snp_info$position), ]
  colnames(LD_matrix) <- rownames(LD_matrix)

  # Reorder LD matrix to match sorted SNPs
  LD_matrix <- LD_matrix[snp_info$probe_id, snp_info$probe_id]
  rownames(LD_matrix) <- 1:nrow(LD_matrix)
  colnames(LD_matrix) <- 1:ncol(LD_matrix)

  # Set upper triangle of LD matrix to NA to avoid duplicate entries
  LD_matrix[lower.tri(LD_matrix)] <- NA
  diag(LD_matrix) <- NA

  # Melt LD matrix for ggplot
  LD_matrix_melted <- melt(as.matrix(LD_matrix), na.rm = TRUE)
  colnames(LD_matrix_melted) <- c("SNP_A", "SNP_B", "R2")

  # Filter out duplicate SNP pairs
  LD_matrix_melted <- LD_matrix_melted[LD_matrix_melted$SNP_A != LD_matrix_melted$SNP_B, ]
  #
  # Generate LD heatmap plot using custom plotting function
  plt <- create_ld_plot(LD_matrix_melted, df_snp = snp_info)
  return(plt)
}


#' Plot LD heatmap from LDs Info table
#' @param LDsInfo a list containing the LD matrices and the info of the clusters
#' @param clusterLDs a list of LD clusters
#' @param snps a vector of SNP names to plot
#' @param outfolder a folder to save the plots
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 ggsave
#' @return A list of ggplot objects showing the LD heatmaps
#' @export
plotLDForClusters <- function(clusterLDs)
{
  outfolder <- get_config("outfolder")
  # create a subfolder to save the LD matrices plots
  outfolder_sub <- file.path(outfolder, "LD_matrices_plots")
  if (!dir.exists(outfolder_sub)) {
    dir.create(outfolder_sub, showWarnings = FALSE, recursive = TRUE)
  }

  # ld_matrix_folder
  ld_matrix_folder<-file.path(outfolder, "LD_matrices")

  # set the outfolder to save the plots in a subfolder
  outfolder <- outfolder_sub

  ldPlots <- list()
  # # if cluster_ids is not provided, plot all clusters
  # if (is.null(cluster_ids)) {
  #   cluster_ids <- names(clusterLDs)
  # }
  for (cls in names(clusterLDs)) {
    # if there is no cluster, skip
    if (length(clusterLDs[[cls]]) == 0) {
      next
    }
    chr <- strsplit(cls, ":")[[1]][1]
    # Read LD matrix
    LD_matrix_file <- file.path(ld_matrix_folder, paste(cls, "ld_matrix.csv", sep = "_"))
    LD_matrix <- read.csv(LD_matrix_file, header = TRUE, row.names = 1)
    colnames(LD_matrix) <- rownames(LD_matrix)
    # get snps in cluster
    LD_matrix <- LD_matrix[clusterLDs[[cls]], clusterLDs[[cls]]]
    happlot <- plotLDheatmap(LD_matrix)
    ldPlots[[length(ldPlots) + 1]] <- happlot
    if (!is.null(outfolder)) {
      ggsave(file.path(outfolder, paste(cls, "ld_heatmap.png", sep = "_")), happlot, width = 10, height = 10, units = "in", dpi = 300)
    }
  }
  return(ldPlots)
}

#' Plot Haplotype Combinations Matrix
#' @param cluster_id cluster id
#' @param haplotypes a data frame containing haplotype information
#' @param gwas a list of GWAS SNPs
#' @param snps a vector of SNP names to plot
#' @param outfolder a folder to save the plots
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 ggsave
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_minimal
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_blank
#' @import ggtext
#' @return a ggplot objects showing the haplotype combinations matrix
#' @export
plotCombMatrix<-function(cluster_id, haplotypes, gwas, snps=NULL, outfolder = NULL)
{
  chr <- strsplit(cluster_id, ":")[[1]][1]
  # Extract haplotypes
  cluster_haplotypes <- haplotypes[haplotypes$snp == cluster_id, ]$clusterComb

  # Extract haplotypes
  cluster_haplotypes <- haplotypes[haplotypes$snp == cluster_id, ]$clusterComb
  snp_names <- haplotypes[haplotypes$snp == cluster_id, ]$snps[1]
  snp_names <- unlist(strsplit(as.character(snp_names), "\\|"))
  # Split the cluster combinations
  cluster_combs <- lapply(cluster_haplotypes, function(x) unlist(strsplit(as.character(x), "\\|")))
  cluster_combs <- do.call(rbind, cluster_combs)
  cluster_combs <- as.data.frame(cluster_combs)
  # Matrix
  rownames(cluster_combs) <- paste0("Comb", 1:nrow(cluster_combs))
  colnames(cluster_combs) <- snp_names

  # Melt
  cluster_combs_melt <- melt(as.matrix(cluster_combs))
  colnames(cluster_combs_melt) <- c("Comb", "SNP", "Genotype")

  # Define the color for axis labels based on SNP presence
  axis_label_colors <- ifelse(unique(cluster_combs_melt$SNP) %in% gwas[[chr]]$rs, "red", "black")
  snp_color_mapping <- setNames(axis_label_colors, unique(cluster_combs_melt$SNP))

  # Plot
  p <- ggplot(cluster_combs_melt, aes(x = SNP, y = Comb, fill = Genotype)) +
    geom_tile() +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
      legend.position = "top",
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.title.x = element_blank()
    ) +
      scale_x_discrete(
        labels = function(x) {
          sapply(x, function(label) {
            if (label %in% names(snp_color_mapping)) {
              color <- snp_color_mapping[[label]]
              paste0("<span style='color:", color, ";'>", label, "</span>")
            } else {
              label
            }
          })
        }
      ) +
      labs(title = "Cluster Haplotypes", x = "SNP", y = "Haplotype Combination") +
      theme(axis.text.x = ggtext::element_markdown(angle = 90, hjust = 1)) # Requires ggtext for markdown rendering
  if (!is.null(outfolder)) {
    ggsave(file.path(outfolder, paste(cluster_id, "comb_matrix.png", sep = "_")), p, width = 10, height = 10, units = "in", dpi = 300)
  }
  return(p)
}

#' Plot LD and Haplotype with Combination Matrix
#' @param ld_matrix_folder a folder containing LD matrices
#' @param clusterLDs a list of LD clusters
#' @param haplotypes a data frame containing haplotype information
#' @param gwas a list of GWAS SNPs
#' @param snps a vector of SNP names to plot
#' @param outfolder a folder to save the plots
#' @importFrom gridExtra arrangeGrob
#' @importFrom ggplot2 ggsave
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom ggplot2 theme
#' @importFrom gridExtra arrangeGrob
#' @importFrom ggplot2 theme_minimal
#' @export

plotLDCombMatrix<-function(clusterLDs, haplotypes, gwas)
{
  outfolder <- get_config("outfolder")
  ld_matrix_folder<-file.path(outfolder, "LD_matrices")

  # create a subfolder to save the LD matrices plots
  outfolder_sub <- file.path(outfolder, "LD_hap_comb_matrices_plots")
  if (!dir.exists(outfolder_sub)) {
    dir.create(outfolder_sub, showWarnings = FALSE, recursive = TRUE)
  }
  outfolder <- outfolder_sub
  # # if cluster_ids is not provided, plot all clusters
  # if (is.null(cluster_ids)) {
  #   cluster_ids <- names(clusterLDs)
  # }

  ld_comb_plots<-list()
  for (cls in names(clusterLDs)) {
    # if there is no cluster, skip
    if (length(clusterLDs[[cls]]) == 0) {
      next
    }
    chr <- strsplit(cls, ":")[[1]][1]
    # Read LD matrix
    LD_matrix_file <- file.path(ld_matrix_folder, paste(cls, "ld_matrix.csv", sep = "_"))

    LD_matrix <- read.csv(LD_matrix_file, header = TRUE, row.names = 1)
    colnames(LD_matrix) <- rownames(LD_matrix)
    # get snps in cluster
    LD_matrix <- LD_matrix[clusterLDs[[cls]], clusterLDs[[cls]]]
    happlot <- plotLDheatmap(LD_matrix)
    # plot the combination matrix
    comb_plot<-plotCombMatrix(cls, haplotypes, gwas, snps = NULL, outfolder = NULL)
    # Arrange and save the combined plot
    arranged_plot <- arrangeGrob(
      grobs = list(comb_plot,grid::grid.grabExpr(grid::grid.draw(happlot))),
      ncol = 1,
      heights = c(0.4, 0.6) # Allocate 40% for upper plot, 5% for spacer, 55% for lower plot
    )

    ld_comb_plots[[cls]]<-arranged_plot
    # save the plot
    if (!is.null(outfolder)) {
      ggsave(file.path(outfolder, paste(cls, "ld_comb_matrix.png", sep = "_")), arranged_plot, width = 10, height = 10, units = "in", dpi = 300)
    }
  }
}
