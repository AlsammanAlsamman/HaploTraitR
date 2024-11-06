#' Convert a habmap object to a numeric matrix
#' @param habmap a habmap object
#' @return a numeric matrix
#' @export

convertGenoBi2Numeric<-function(habmap)
{
  # remove the first 11 columns
  habmap<-habmap[,-c(1:11)]
  # convert the data to matrix
  habmap<-as.matrix(habmap)
  binaryMatrix<-.Call("convertGenoBi2Numeric_c", habmap)
  return(binaryMatrix)
}
