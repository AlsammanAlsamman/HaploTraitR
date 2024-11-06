#' Get samples for a haplotype combination
#' @param haplotypes A list of haplotype combinations
#' @param hapmap A hapmap object
#' @return A list of samples for each haplotype combination
#' @export
########## Need to be more optimized
getHapCombSamples<-function(haplotypes, hapmap)
{
  snpsnames<-names(haplotypes)
  hap_comb_samples<-matrix(ncol= 4)
  for(snpname in snpsnames)
  {

    chr<-haplotypes[[snpname]][1,]$chr
    for(i in 1:nrow(haplotypes[[snpname]]))
    {
      combstr<-as.character(haplotypes[[snpname]][i,]$clusterComb)
      comb<-unlist(strsplit(as.character(combstr), "\\|"))
      snps<-unlist(strsplit(as.character(haplotypes[[snpname]][i,]$snps), "\\|"))
      snpsdata<-hapmap[[chr]][snps,]
      # replace values not exactly equal to comb value to NA
      for(i in 1:length(comb))
      {
        snpvalue<-comb[i]
        # in the i snp in data, replace all values not equal to comb[i] to NA
        snpsdata[i,]<-ifelse(snpsdata[i,]!=snpvalue, NA, snpsdata[i,])
      }
      # remove columns with one or more NA
      snpsdata<-snpsdata[,colSums(is.na(snpsdata))==0]
      # get the samples
      samples<-colnames(snpsdata)
      samples<-paste(samples, collapse="|")
      # create a vector of samples concatenated by |, snpname, and chr
      # snpname, chr, pos, combstr, samples
      hap_comb_samples<-rbind(hap_comb_samples, c(snpname, chr, combstr, samples))
    }

  }
  hap_comb_samples<-hap_comb_samples[-1,]
  colnames(hap_comb_samples)<-c("snp", "chr", "comb", "samples")
  return(hap_comb_samples)
}




# combstr<-as.character(haplotypes[[snpname]][1,]$clusterComb)
# comb<-unlist(strsplit(as.character(combstr), "\\|"))
# snps<-unlist(strsplit(as.character(haplotypes[[snpname]][1,]$snps), "\\|"))
# snpsdata<-hapmap[[chr]][snps,]
# # replace values not exactly equal to comb value to NA
# for(i in 1:length(comb))
# {
#   snpvalue<-comb[i]
#   # in the i snp in data, replace all values not equal to comb[i] to NA
#   snpsdata[i,]<-ifelse(snpsdata[i,]!=snpvalue, NA, snpsdata[i,])
# }
# # remove columns with one or more NA
# snpsdata<-snpsdata[,colSums(is.na(snpsdata))==0]
# # get the samples
# samples<-colnames(snpsdata)
# samples<-paste(samples, collapse="|")
# # create a vector of samples concatenated by |, snpname, and chr
# # snpname, chr, pos, combstr, samples
# hap_comb_samples<-rbind(hap_comb_samples, c(snpname, chr, combstr, samples))

