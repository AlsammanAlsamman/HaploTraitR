#' Get the phenotype and genotype data
#' @param subhapmap A subset of the hapmap data containing the significant SNPs
#' @param pheno A phenotype data frame
#' @return A data frame with the phenotype and genotype data
#' @export
#' @importFrom reshape2 melt
get_pheno_geno<-function(subhapmap, pheno){
  colnames(pheno)<-c("Taxa", "Phenotype")
  # rbinding the hapmap
  # remove names
  names(subhapmap)<-NULL
  subhapmap <- do.call(rbind, subhapmap)
  # remove the first 11 columns
  subhapmap[,1:11]<-NULL
  # transpose the data
  subhapmap<-t(subhapmap)
  snps_target<-colnames(subhapmap)
  # merge the pheno data with the hapmap
  geno_pheno<-merge(subhapmap, pheno, by.x="row.names", by.y="Taxa")
  # transfere row.names to the first column
  rownames(geno_pheno)<-geno_pheno[,1]
  geno_pheno<-geno_pheno[,-1]
  geno_pheno_melt<-melt(geno_pheno, id.vars=c("Phenotype"))
  return(geno_pheno_melt)
}
