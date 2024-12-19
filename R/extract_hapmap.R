#' Extract hapmap data for significant SNPs
#' @param hapmap A list of hapmap data frames
#' @param gwas A list of GWAS data frames
#' @return A list of hapmap data frames for significant SNPs
#' @export
extract_hapmap<-function(hapmap, gwas){
  sub_hapmap <- list()
  for (chromosome in names(gwas)) {
    gwas_chr <- gwas[[chromosome]]
    hapmap_chr <- hapmap[[chromosome]]
    gwas_rs_id<-paste(gwas_chr$chr, gwas_chr$pos, sep=":")
    # Extract hapmap data for significant SNPs
    sub_hapmap_chr <- hapmap_chr[gwas_rs_id, ]
    sub_hapmap[[chromosome]] <- sub_hapmap_chr
  }
  return(sub_hapmap)
}
