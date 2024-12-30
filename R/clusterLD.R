#' Cluster SNPs based on LD matrices
#'
#' @param LDsInfo a list containing the LD matrices and the info of the clusters
#' @importFrom igraph graph_from_adjacency_matrix components
#' @export
#' @return LD cluster information
getLDclusters <- function(LDsInfo) {
  # Retrieve configuration parameters
  ld_threshold <- get_config("ld_threshold")
  cls_count <- get_config("dist_cluster_count")

  out_info<-LDsInfo$out_info
  ld_folder<-LDsInfo$ld_folder
  cluster_info<-list()
  for(i in 1:nrow(out_info))
  {
    chr<-out_info[i,][1]
    cls<-out_info[i,][2]
    ld_matrix<-read.csv(file.path(ld_folder, paste(cls, "ld_matrix.csv", sep="_")), row.names=1)
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
  #names(cluster_info)<-paste(out_info[,1], out_info[,2], sep="_")
  names(cluster_info)<-out_info[,2]
  # delete clusters where the significant snp is not there
  for(i in 1:length(cluster_info))
  {
    main_snp <-names(cluster_info)
    # which cluster is the main snp is
    main_snp_cluster<-cluster_info[[i]][cluster_info[[i]]$SNP==main_snp[i],]$cluster
    # remove all clusters where the main snp is not there
    cluster_info[[i]]<-cluster_info[[i]][cluster_info[[i]]$cluster==main_snp_cluster,]$SNP
  }
  return(cluster_info)
}
