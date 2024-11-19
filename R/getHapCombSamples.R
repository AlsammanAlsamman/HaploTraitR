#' Get samples for a haplotype combination
#' @param haplotypes A list of haplotype combinations
#' @param hapmap A hapmap object
#' @return A list of samples for each haplotype combination
#' @export
########## Need to be more optimized
getHapCombSamples<-function(haplotypes, hapmap)
{
  # add column for haplotypes
  haplotypes$samples<-NA
  for(snp_cls_i in 1:nrow(haplotypes))
  {
    snp_cls<-haplotypes[snp_cls_i,]
    snp_cls_chr<-snp_cls$chr
    snp_cls_combstr<-snp_cls$clusterComb
    snp_cls_comb<-unlist(strsplit(as.character(snp_cls_combstr), "\\|"))
    snp_cls_snps<-unlist(strsplit(as.character(snp_cls$snps), "\\|"))
    snp_cls_data<-hapmap[[snp_cls_chr]][snp_cls_snps,]
    for(i in 1:length(snp_cls_comb))
    {
      snpvalue<-snp_cls_comb[i]
      snp_cls_data[i,]<-ifelse(snp_cls_data[i,]!=snpvalue, NA, snp_cls_data[i,])
    }
    snp_cls_data<-snp_cls_data[,colSums(is.na(snp_cls_data))==0]
    samples<-colnames(snp_cls_data)
    samples<-paste(samples, collapse="|")
    # add samples to haplotypes
    haplotypes[snp_cls_i,]$samples<-samples
  }
  # add combination 0, which is the samples that do not have any of the haplotypes for the snps clusters
  snp_clss<-unique(haplotypes$snp)
  hap_samples <- colnames(hapmap[[1]])[12:ncol(hapmap[[1]])]

  for (snp_cls in snp_clss)
  {
    chr<-strsplit(snp_cls, ":")[[1]][1]
    allsamples<-unique(unlist(strsplit(haplotypes[haplotypes$snp==snp_cls,]$samples, "\\|")))
    total_comp_freq<-sum(haplotypes[haplotypes$snp==snp_cls,]$freq)
    no_comb_samples<-setdiff(hap_samples, allsamples)
    no_comb_samples<-paste(no_comb_samples, collapse="|")
    # add samples to haplotypes
    haplotypes<-rbind(haplotypes, data.frame(clusterComb="None", Freq=1-total_comp_freq, chr=chr, snps="",snp=snp_cls, comb="0", samples=no_comb_samples))
  }
  return(haplotypes)
}
