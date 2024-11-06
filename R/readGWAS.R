#' read gwas results file rs, chr, pos, p
#' @param gwasfile Path to the GWAS file
#' @return A data frame with the GWAS data
#' @export
#' @examples
#' gwasfile<-"SignificantSNP_GWAS.csv"
#' gwas<-readGWAS(gwasfile)
#' print(gwas)
readGWAS<-function(gwasfile,header=TRUE,sep="\t")
{
  gwas<-read.csv(gwasfile,header=header,sep=sep)
  colnames(gwas)<-c("rs","chr","pos","p")
  # split by chr
  gwas<-split(gwas, gwas$chr)
  return(gwas)
}
