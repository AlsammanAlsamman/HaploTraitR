#' Compute LD matrix for a cluster
#' @param hapmap a hapmap object containing the genotypic data
#' @param haplotype_clusters a list of haplotype clusters
#' @param chr the chromosome
#' @param cls the cluster
#' @param outfolder the output folder
computeLDclusterByChr<-function(hapmap, haplotype_clusters, chr, cls, outfolder)
{
  geno_cluster<-haplotype_clusters[[chr]][[cls]]
  # lets extract the gentypic data of this cluster
  geno_cluster_data<-hapmap[[chr]][paste(chr, geno_cluster, sep=":"),]
  # get the snps
  cluster.snps<-rownames(geno_cluster_data)
  # convert the data to numeric
  geno_cluster_data_numeric<-convertGenoBi2Numeric(geno_cluster_data)
  # calculate the LD matrix
  ld_matrix <- calculate_ld_matrix(geno_cluster_data_numeric)
  rownames(ld_matrix)<-cluster.snps
  colnames(ld_matrix)<-cluster.snps
  # save the LD matrix
  ld_matrix_file<-file.path(outfolder, paste(cls, "ld_matrix.csv", sep="_"))
  write.csv(ld_matrix, ld_matrix_file)
}


#' Compute LD matrix for a clusters
#' @param hapmap a hapmap object containing the genotypic data
#' @param haplotype_clusters a list of haplotype clusters
#' @param outfolder the output folder
#' @return a folder containing the LD matrices
#' @export

computeLDclusters<-function(hapmap, haplotype_clusters, outfolder)
{
  # create a temporary folder to store LD matrices
  if(is.null(outfolder))
  {
    outfolder<-tempdir()
    print("A temporary folder will be created to store the LD matrices")
  }

  # create a folder to store the LD matrices
  # create the output folder if it does not exist
  if(!dir.exists(outfolder))
  {
    dir.create(outfolder, showWarnings = FALSE)
  }
  ld_folder<-file.path(outfolder, "LD_matrices")
  dir.create(ld_folder, showWarnings = FALSE)
  out_info<-list()
  for(chr in names(haplotype_clusters))
  {
    for(cls in names(haplotype_clusters[[chr]]))
    {
      out_info[[length(out_info)+1]]<-c(chr, cls)
      computeLDclusterByChr(hapmap, haplotype_clusters, chr, cls, ld_folder)
    }
  }
  # convert the list to a data frame
  out_info<-do.call(rbind, out_info)
  # return list
  out_list<-list(ld_folder, out_info)
  names(out_list)<-c("ld_folder", "out_info")
  return(out_list)
}
