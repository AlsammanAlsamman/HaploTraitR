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
  .haplotraitr_env$pheno_col <- "Phenotype"
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

#' list all the configuration parameters
#' @export
list_config <- function() {
  ls.str(.haplotraitr_env)
}

#' Save configuration to a file
#' @param filename The name of the file to save the configuration to
#' @export
#' @examples
#' save_config("config.txt")
save_config <- function(filename) {
  config <- as.list(.haplotraitr_env)
  # save as text file
  for (name in names(config)) {
    write(paste(name, "=", config[[name]], sep = ""), file = filename, append = TRUE)
  }
  # save as RDS file
  saveRDS(config, file = paste0(filename, ".rds"))
  print("Configuration saved.")
}

#' Load configuration from a file
#' @param filename The name of the file to load the configuration from
#' @export
#' @examples
#' load_config("config.rds")
#' load_config("config.txt")
load_config <- function(filename) {
  config<-list()
  if (grepl(".rds$", filename)) {
    config <- readRDS(filename)}
  else if (grepl(".txt$", filename))
  {
    config <- read.table(filename, sep = "=", col.names = c("name", "value"), stringsAsFactors = FALSE)
    config <- as.list(config$value)
  }
  set_config(config)
  print("Configuration loaded.")
}

#' reset the configuration to default
#' @export
#' @examples
#' reset_config()
reset_config <- function() {
  initialize_config()
  print("Configuration reset to default.")
}

