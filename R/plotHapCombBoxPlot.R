#' Plotting Haplotype Combination Boxplot for specific SNP
#' @param cls_snp The SNP of interest, if NULL then plots for all SNPs in t_test_snpComp are generated
#' @param SNPcombTables A list of data frames containing haplotype combinations and phenotype values
#' @param t_test_snpComp A list of t-test results for haplotype combinations
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_violin
#' @importFrom ggplot2 theme_minimal
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 position_jitter
#' @importFrom ggplot2 geom_point
#' @import ggpubr
#' @export
plotHapCombBoxPlot <- function(cls_snp = NULL, SNPcombTables, t_test_snpComp, with_dots = TRUE, outfolder) {
  if(is.null(cls_snp)){
    for (i in names(t_test_snpComp)) {
      comb_sample <- SNPcombTables[SNPcombTables$SNP == i, ]
      comb_sample$comb <- paste0("comb", comb_sample$comb)
      p <- ggplot(comb_sample, aes(x = comb, y = Pheno)) +
        geom_violin(aes(fill = comb), alpha = 0.5) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = paste("Phenotype Distribution by Haplotype Combination for SNP", i),
             x = "Haplotype Combination",
             y = "Phenotype Value") +
        theme(legend.position = "none")
    
      if (with_dots) {
        p <- p + geom_point(aes(x = comb, y = Pheno, color=Pheno), position = position_jitter(width = 0.1), alpha = 0.5)
      }
    
      if ("p.adj.signif" %in% colnames(t_test_snpComp[[i]])) {
        p <- p + stat_pvalue_manual(t_test_snpComp[[i]], label = "p.adj.signif", tip.length = 0.01)
      }
    
      p
      ggplot2::ggsave(file.path(outfolder, paste0(i, "_boxplot.png")), width = 8, height = 6)
    }
  }else{
    comb_sample <- SNPcombTables[SNPcombTables$SNP == cls_snp, ]
    comb_sample$comb <- paste0("comb", comb_sample$comb)
    p <- ggplot(comb_sample, aes(x = comb, y = Pheno)) +
      geom_violin(aes(fill = comb), alpha = 0.5) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste("Phenotype Distribution by Haplotype Combination for SNP", cls_snp),
           x = "Haplotype Combination",
           y = "Phenotype Value") +
      theme(legend.position = "none")
  
    if (with_dots) {
      p <- p + geom_point(aes(x = comb, y = Pheno, color=Pheno), position = position_jitter(width = 0.1), alpha = 0.5)
    }
  
    if ("p.adj.signif" %in% colnames(t_test_snpComp[[cls_snp]])) {
      p <- p + stat_pvalue_manual(t_test_snpComp[[cls_snp]], label = "p.adj.signif", tip.length = 0.01)
    }
  
  p
  ggplot2::ggsave(file.path(outfolder, paste0(cls_snp, "_boxplot.png")), width = 8, height = 6)
  }
}


# # Extract the data frame for the specified SNP
# cls<-strsplit(cls_snp, "@")[[1]][1]
# snp<-strsplit(cls_snp, "@")[[1]][2]
# comb_sample <- SNPcombTables[SNPcombTables$ldcls==cls & SNPcombTables$SNP==snp,]
# comb_sample$comb<-paste0("comb", comb_sample$comb)
# # Create a boxplot with significance annotations
# p <- ggplot(comb_sample, aes(x = comb, y = Pheno)) +
#   geom_violin(aes(fill = comb), alpha = 0.5) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   labs(title = paste("Phenotype Distribution by Haplotype Combination for SNP", snp),
#        x = "Haplotype Combination",
#        y = "Phenotype Value") +
#   theme(legend.position = "none")
#
# if (with_dots) {
#   p <- p + geom_point(aes(x = comb, y = Pheno), position = position_jitter(width = 0.1), alpha = 0.5)
# }
# # Add significance brackets to the plot
# # if the snp has a t test result
# # if there p.adj.signif column in the t test result
# if ("p.adj.signif" %in% colnames(t_test_snpComp[[cls_snp]])) {
#   p <- p + stat_pvalue_manual(t_test_snpComp[[cls_snp]], label = "p.adj.signif", tip.length = 0.01)
# }
# # Display the plot
# p

