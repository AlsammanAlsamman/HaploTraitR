#' Get SNP combinations sample tables
#' @param haplotypes A list of haplotype combinations
#' @param pheno A phenotype data frame
#' @param savecopy A boolean indicating whether to save a copy of the sample tables to the project folder+
#' @return A list of sample tables for each SNP combination
#' @export
getSNPcombTables <- function(haplotypes, pheno, savecopy=TRUE){
  comb_sample_Tables<-list()
  colnames(pheno)<-c("Sample", "Pheno")
  for(snp_cls_i in 1:nrow(haplotypes))
  {
    snp_cls<-haplotypes[snp_cls_i,]
    snp_cls_name<-snp_cls$cluster
    snp_cls_chr<-snp_cls$chr
    snp_cls_snp<-snp_cls$snp
    snp_cls_comb<-snp_cls$comb
    snp_cls_samples<-snp_cls$samples
    snp_cls_samples<-unlist(strsplit(as.character(snp_cls_samples), "\\|"))
    snp_cls_samples_pheno<-pheno[pheno$Sample %in% snp_cls_samples,]$Pheno
    # create a data frame
    snp_cls_samples<-data.frame(Sample=snp_cls_samples,
                                Pheno=snp_cls_samples_pheno,
                                ldcls=snp_cls_name,
                                comb=snp_cls_comb,
                                SNP=snp_cls_snp)
    # append to comb_sample_Tables
    comb_sample_Tables[[length(comb_sample_Tables)+1]]<-snp_cls_samples
  }
  # bind
  comb_sample_Tables<-do.call(rbind, comb_sample_Tables)
  if(savecopy){
    write.csv(comb_sample_Tables, file.path(get_config("outfolder"), "SNP_comb_sample_Tables.csv"))
    print("Saved SNP_comb_sample_Tables.csv")
  }
  return(comb_sample_Tables)
}
