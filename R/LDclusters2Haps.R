#' Convert LD clusters to haplotypes and calculate the frequency of each haplotype combination
#'@param hapmap A hapmap object
#'@param clusterLDs A list of LD clusters
#'@param comb_freq_threshold The minimum frequency of a haplotype combination
#'@return A list of haplotype combinations with their frequencies
#'@export
convertLDclusters2Haps<-function(hapmap, clusterLDs, comb_freq_threshold)
{
  hap_comb<-list()
  # loop over the clusters
  for (i in 1:length(clusterLDs))
  {
    cluster<-clusterLDs[[i]]
    chr<-strsplit(cluster[[1]][1], ":")[[1]][1]
    # select the data of these snps
    cluster.data<-hapmap[[chr]][cluster$SNP,]
    # remove all non variant columns first 11 columns
    cluster.data<-cluster.data[,-c(1:11)]
    # transpose the data
    cluster.data<-t(cluster.data)
    # concatenate the data in one column
    clusterComb<-apply(cluster.data, 1, paste, collapse="|")
    # calculate the frequency of each combination
    snpsComb<-as.data.frame(table(clusterComb))
    totalsample<-sum(snpsComb$Freq)
    snpsComb$Freq<-snpsComb$Freq/totalsample
    snpsComb$snps<-paste(cluster$SNP, collapse="|")
    snpsComb$chr<-chr
    # remove rows with frequency less than comb_freq_threshold
    snpsComb<-snpsComb[snpsComb$Freq>comb_freq_threshold,]
    hap_comb[[i]]<-snpsComb
  }
  # add names to the list
  names(hap_comb)<-names(clusterLDs)
  return(hap_comb)
}
