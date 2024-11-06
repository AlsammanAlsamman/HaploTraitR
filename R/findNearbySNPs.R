#' Find nearby SNPs to a given SNP list assuming they are in the same chromosome
#' @param snp_list1 A list of SNPs (Reference)
#' @param snp_list2 A list of SNPs (Query)
#' @param threshold The maximum distance between the SNPs
#' @return A list of SNPs from snp_list2 that are within the threshold distance from snp_list1
#' @export
#' @examples
#' list1 <- c(1000, 2000, 3000)
#' list2 <- c(1500, 2500, 3500, 5000)
#' threshold <- 1000
#' result <- find_nearby_snps(list1, list2, threshold)
#' print(result)
find_nearby_snps<-function(snp_list1, snp_list2, threshold)
{
  result <- .Call("find_nearby_snps_c", as.integer(snp_list1), as.integer(snp_list2), as.integer(threshold))
  return(result)
}
