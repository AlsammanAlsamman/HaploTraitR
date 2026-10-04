#' HaploTraitR: Haplotype-Based Trait Analysis for Plant Breeding
#'
#' HaploTraitR turns GWAS hits into haplotype blocks (SNPs linked by LD),
#' tests how each haplotype affects a trait, and reports the superior
#' haplotypes, diagnostic markers and candidate parents in plots and a
#' formatted Excel workbook.
#'
#' The quickest way to use the package is [run_haplotraitr()]. The individual
#' steps ([read_gwas_file()], [readHapmap()], [getDistClusters()], ...) remain
#' available for custom workflows.
#'
#' @keywords internal
"_PACKAGE"

#' @import ggplot2
#' @importFrom stats aggregate cor kruskal.test median oneway.test p.adjust qt sd setNames t.test wilcox.test
#' @importFrom utils combn head read.table write.csv
#' @importFrom grDevices dev.off pdf png
NULL

utils::globalVariables(c(
  "x", "y", "xend", "yend", "fill", "label", "group", "id", "value", "R2",
  "Pheno", "hap", "hap_label", "mean", "ci_low", "ci_high", "role", "freq_pct",
  "pct_vs_pop", "block_label", "pos_cum", "logp", "chr_col", "is_lead",
  "Genotype", "status", "trait", "accession", "xmin", "xmax", "ymin", "ymax",
  "y.position", "p.adj.signif", "size", "Phenotype", "block", "score",
  "lo", "hi", "pos", "start", "end"
))
