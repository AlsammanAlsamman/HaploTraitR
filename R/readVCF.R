#' Read a VCF file
#' @param vcfpath Path to the VCF file
#' @return A data frame with the VCF data
#' @export
readVCF<-function(vcfpath)
{
 return(read.csv(vcfpath,sep="\t",header=FALSE, comment.char="#"))
}
