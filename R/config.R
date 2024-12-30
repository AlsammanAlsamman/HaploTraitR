# Create a package environment to store configuration settings
.haplotraitr_env <- new.env(parent = emptyenv())

# Function to initialize the default configuration
initialize_config <- function() {
  .haplotraitr_env$outfolder <- NULL
  .haplotraitr_env$dist_threshold <- 1000000
  .haplotraitr_env$dist_cluster_count <- 5
  .haplotraitr_env$ld_threshold <- 0.3
  .haplotraitr_env$comb_freq_threshold <- 0.1
  .haplotraitr_env$fdr_threshold <- 0.05
  .haplotraitr_env$fdr_method <- "fdr"
  .haplotraitr_env$t_test_threshold <- 0.05
  .haplotraitr_env$t_test_method <- "t.test"
  .haplotraitr_env$rsid_col <- "rsid"
  .haplotraitr_env$pos_col <- "pos"
  .haplotraitr_env$chr_col <- "chr"
  .haplotraitr_env$pval_col <- "p"
  .haplotraitr_env$fdr_col <- "fdr"
  .haplotraitr_env$phenotypename <- "Phenotype"
  .haplotraitr_env$phenotypeunit <- NULL
}

# Initialize configuration when the package is loaded
.onLoad <- function(libname, pkgname) {
  initialize_config()
  print("HaploTraitR configuration initialized.")
  print("In order to change the configuration, use the set_config function.")
}

#' Set configuration parameters for HaploTraitR
#'
#' @param config A named list of configuration parameters.
#' @export
#' @examples
#' config <- list(dist_threshold = 2000000, dist_cluster_count = 10)
#' set_config(config)
set_config <- function(config = list()) {
  if (!is.list(config)) {
    stop("The config parameter must be a list.")
  }

  for (name in names(config)) {
    if (exists(name, envir = .haplotraitr_env)) {
      assign(name, config[[name]], envir = .haplotraitr_env)
    } else {
      warning(paste("Unknown configuration parameter:", name))
    }
  }
}

#' Get a configuration parameter
#'
#' @param name The name of the configuration parameter
#' @return The value of the configuration parameter
#' @export
#' @examples
#' get_config("dist_threshold")
get_config <- function(name) {
  if (exists(name, envir = .haplotraitr_env)) {
    get(name, envir = .haplotraitr_env)
  } else {
    stop(paste("Configuration parameter", name, "not found."))
  }
}
