#' Barley leaf area example data
#'
#' GWAS results, genotypes and phenotypes for leaf area in a barley panel
#' (275 accessions, Barley 50K SNP array). The same data are available as
#' files in `system.file("extdata", package = "HaploTraitR")`.
#'
#' @format A list with three elements:
#' \describe{
#'   \item{gwas}{GWAS results (data frame) for leaf area}
#'   \item{hapmap}{Genotypes as returned by [readHapmap()] (list of data frames per chromosome)}
#'   \item{pheno}{Phenotypes: accession name (`Taxa`) and leaf area (`Area`, cm2)}
#' }
#' @source ICARDA barley breeding program.
"barley_area"
