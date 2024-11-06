#' Read hapmap file into R
#' @param path Path to hapmap file
#' @param sep Separator
#' @param header Logical indicating if the file has a header
#' @param ... Additional arguments to \code{\link{read.table}}
#' @return A data frame
#' @export
#' @examples
#' readHapmap("hapmap.txt")
#' readHapmap("hapmap.txt", header = FALSE)
#' readHapmap("hapmap.txt", header = FALSE, sep = "\t")
#' readHapmap("hapmap.txt", header = FALSE, sep = "\t", skip = 1)
#' readHapmap("hapmap.txt", header = FALSE, sep = "\t", skip = 1, nrows = 10)
readHapmap <- function(path, sep = "\t", header = TRUE, ...) {
  # remove the # character from the header
  hap<-read.csv(path, sep = sep, header = header, ...)
  # create rs column
  hap$rsid<-paste(hap$chrom, hap$pos, sep=":")
  # check for duplicated rsid if there are any remove them and print a message
  if(sum(duplicated(hap$rsid))>0)
  {
    cat("warning: duplicated rsid found, it will be removed!\n")
    hap<-hap[!duplicated(hap$rsid),]
  }
  # make the rsid the rownames
  rownames(hap)<-hap$rsid
  # remove the rsid column
  hap<-hap[,-ncol(hap)]
  # split by chromosome
  hap<-split(hap, hap$chrom)
  return(hap)
}
