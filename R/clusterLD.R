#' Cluster SNPs based on LD matrices
#' @param LDsInfo a list containing the LD matrices and the info of the clusters
#' @param ld_threshold The minimum LD value to consider two SNPs in LD
#' @importFrom igraph graph_from_adjacency_matrix
#' @return LD cluster information
clusterLD<-function(LDsInfo, ld_threshold,cls_count=3)
{
  out_info<-LDsInfo$out_info
  ld_folder<-LDsInfo$ld_folder
  cluster_info<-list()
  for(i in 1:nrow(out_info))
  {
    chr<-out_info[i,][1]
    cls<-out_info[i,][2]
    ld_matrix<-read.csv(file.path(ld_folder, paste(chr, cls, "ld_matrix.csv", sep="_")), row.names=1)
    ld_matrix<-as.matrix(ld_matrix)
    colnames(ld_matrix)<-rownames(ld_matrix)
    snps<-rownames(ld_matrix)
    # if the value is less than 0.3, set it to 0
    ld_matrix[ld_matrix<ld_threshold]<-0
    colnames(ld_matrix)<-1:ncol(ld_matrix)
    rownames(ld_matrix)<-1:ncol(ld_matrix)
    # convert the matrix to a graph
    g<-graph_from_adjacency_matrix(ld_matrix, mode="lower", weighted=TRUE)
    # get the clusters
    clusters<-components(g)
    # get cluster membership
    cluster_membership<-clusters$membership
    #convert the membership to a data frame
    cluster_membership<-data.frame(SNP=snps[1:ncol(ld_matrix)], cluster=cluster_membership)
    # remove clusters with less than 3 SNPs
    # count the number of SNPs in each cluster
    cluster_count<-table(cluster_membership$cluster)
    # select clusters with more than 3 SNPs
    cluster_membership<-cluster_membership[cluster_membership$cluster %in% names(cluster_count)[cluster_count>=cls_count],]
    cluster_info[[length(cluster_info)+1]]<-cluster_membership
  }
  names(cluster_info)<-paste(out_info[,1], out_info[,2], sep="_")
  return(cluster_info)
}
