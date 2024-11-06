#' save LD matrices to a folder
#' @param LDsInfo a list containing the LD matrices and the info
#' @param outfolder the folder to save the LD matrices
#' @export

saveLDs2folder<-function(LDsInfo, outfolder)
{
  ld_folder<-LDsInfo$ld_folder
  out_info<-LDsInfo$out_info
  # save to outfolder
  for(ldfile in list.files(ld_folder))
  {
    file.copy(file.path(ld_folder, ldfile), outfolder)
  }
  # save the info
  write.csv(out_info, file.path(outfolder, "ld_info.csv"))
}
