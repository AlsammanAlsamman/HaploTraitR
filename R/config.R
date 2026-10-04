# Package environment that stores the configuration settings
.haplotraitr_env <- new.env(parent = emptyenv())

# Default values of every configuration parameter
.config_defaults <- function() {
  list(
    outfolder = NULL,
    dist_threshold = 1000000,
    dist_cluster_count = 5,
    ld_threshold = 0.3,
    comb_freq_threshold = 0.1,
    hap_max_missing = 0.2,
    fdr_threshold = 0.05,
    fdr_method = "fdr",
    t_test_threshold = 0.05,
    t_test_method = "t.test",
    p_adjust_method = "holm",
    higher_is_better = TRUE,
    rsid_col = "rsid",
    pos_col = "pos",
    chr_col = "chr",
    pval_col = "p",
    fdr_col = "fdr",
    phenotypename = "Phenotype",
    pheno_col = "Phenotype",
    phenotypeunit = NULL,
    plot_format = "png",
    plot_dpi = 300
  )
}

# One-line description of every configuration parameter
.config_help <- c(
  outfolder = "Output folder for all results",
  dist_threshold = "Window (bp) around each GWAS SNP in which SNPs are collected",
  dist_cluster_count = "Minimum number of SNPs in a window / LD block",
  ld_threshold = "Minimum LD (r2) to link two SNPs into the same block",
  comb_freq_threshold = "Minimum frequency of a haplotype (rarer ones are pooled as 'Other')",
  hap_max_missing = "Maximum fraction of missing calls allowed when assigning an accession to a haplotype",
  fdr_threshold = "FDR threshold used to select significant GWAS SNPs",
  fdr_method = "Method used by p.adjust() to compute GWAS FDR",
  t_test_threshold = "Significance level (alpha) for haplotype tests",
  t_test_method = "Test used to compare haplotypes: 't.test' (Welch) or 'wilcox.test'",
  p_adjust_method = "Multiple-testing correction for pairwise haplotype tests",
  higher_is_better = "TRUE if a higher trait value is desirable (e.g. yield), FALSE otherwise (e.g. disease score)",
  rsid_col = "Marker name column in the GWAS file",
  pos_col = "Position column in the GWAS file",
  chr_col = "Chromosome column in the GWAS file",
  pval_col = "P-value column in the GWAS file",
  fdr_col = "FDR column in the GWAS file (computed if absent)",
  phenotypename = "Trait name used in titles and file names",
  pheno_col = "Trait column in the phenotype file",
  phenotypeunit = "Trait unit shown on plot axes",
  plot_format = "Plot file format: 'png' or 'pdf'",
  plot_dpi = "Resolution of PNG plots"
)

# Function to initialize the default configuration
initialize_config <- function() {
  rm(list = ls(.haplotraitr_env), envir = .haplotraitr_env)
  defaults <- .config_defaults()
  for (name in names(defaults)) {
    assign(name, defaults[[name]], envir = .haplotraitr_env)
  }
}

.onLoad <- function(libname, pkgname) {
  initialize_config()
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(cli::format_inline(
    "{.strong HaploTraitR} {utils::packageVersion('HaploTraitR')} \u00b7 start with {.fun run_haplotraitr}; ",
    "{.fun list_config} shows the settings."))
}

#' Set configuration parameters for HaploTraitR
#'
#' @param config A named list of configuration parameters. See [list_config()]
#'   for the available parameters.
#' @return Invisibly, the previous values of the changed parameters.
#' @export
#' @examples
#' set_config(list(dist_threshold = 2000000, dist_cluster_count = 10))
#' reset_config()
set_config <- function(config = list()) {
  if (!is.list(config)) {
    stop("The config parameter must be a list, e.g. set_config(list(ld_threshold = 0.5)).")
  }
  old <- list()
  for (name in names(config)) {
    if (exists(name, envir = .haplotraitr_env, inherits = FALSE)) {
      old[name] <- list(get(name, envir = .haplotraitr_env))
      assign(name, config[[name]], envir = .haplotraitr_env)
    } else {
      warning("Unknown configuration parameter: ", name,
              ". Run list_config() to see the valid parameters.")
    }
  }
  invisible(old)
}

#' Get a configuration parameter
#'
#' @param name The name of the configuration parameter
#' @return The value of the configuration parameter
#' @export
#' @examples
#' get_config("dist_threshold")
get_config <- function(name) {
  if (exists(name, envir = .haplotraitr_env, inherits = FALSE)) {
    get(name, envir = .haplotraitr_env)
  } else {
    stop("Configuration parameter '", name, "' not found. Run list_config() to see the valid parameters.")
  }
}

#' List all the configuration parameters
#'
#' @return A data frame with the parameter names, current values and descriptions
#'   (printed and returned invisibly).
#' @export
#' @examples
#' list_config()
list_config <- function() {
  out <- list_config_df()
  defaults <- vapply(.config_defaults(), format_config_value, character(1))
  changed <- out$value != defaults[out$parameter]
  name_w <- max(nchar(out$parameter)) + 2
  value_w <- min(max(nchar(out$value)), 30) + 2
  cli::cli_h2("HaploTraitR settings")
  for (i in seq_len(nrow(out))) {
    value <- if (changed[i]) cli::col_yellow(cli::style_bold(out$value[i])) else cli::style_bold(out$value[i])
    cli::cli_verbatim(paste0("  ", pad(cli::col_cyan(out$parameter[i]), name_w),
                             pad(value, value_w), muted(out$description[i])))
  }
  cli::cli_text(muted("Changed values are yellow. Use set_config(list(name = value)) to change a setting."))
  invisible(out)
}

format_config_value <- function(v) {
  if (is.null(v)) return("NULL")
  if (is.numeric(v)) return(paste(format(v, big.mark = ",", scientific = FALSE, trim = TRUE), collapse = ", "))
  paste(format(v), collapse = ", ")
}

list_config_df <- function() {
  names <- names(.config_defaults())
  values <- vapply(names, function(n) {
    v <- get(n, envir = .haplotraitr_env)
    format_config_value(v)
  }, character(1))
  data.frame(parameter = names, value = values, description = unname(.config_help[names]),
             row.names = NULL, stringsAsFactors = FALSE)
}

#' Save configuration to a file
#'
#' Writes a human-readable `name=value` text file and an `.rds` copy that
#' preserves the value types.
#' @param filename The name of the text file to save the configuration to
#' @return Invisibly, the file name.
#' @export
#' @examples
#' f <- tempfile(fileext = ".txt")
#' save_config(f)
save_config <- function(filename) {
  config <- mget(names(.config_defaults()), envir = .haplotraitr_env)
  lines <- vapply(names(config), function(n) {
    v <- config[[n]]
    paste0(n, "=", if (is.null(v)) "NULL" else paste(v, collapse = ","))
  }, character(1))
  writeLines(lines, filename)
  saveRDS(config, file = paste0(filename, ".rds"))
  cli::cli_alert_success("Configuration saved to {.file {filename}}")
  invisible(filename)
}

#' Load configuration from a file
#'
#' @param filename A file written by [save_config()] (`.txt` or `.rds`).
#' @return Invisibly, the loaded configuration list.
#' @export
#' @examples
#' f <- tempfile(fileext = ".txt")
#' save_config(f)
#' load_config(f)
load_config <- function(filename) {
  if (grepl("\\.rds$", filename, ignore.case = TRUE)) {
    config <- readRDS(filename)
  } else {
    lines <- readLines(filename)
    lines <- lines[nzchar(lines)]
    keys <- sub("=.*$", "", lines)
    vals <- sub("^[^=]*=", "", lines)
    config <- lapply(vals, function(v) {
      if (v == "NULL") NULL else utils::type.convert(v, as.is = TRUE)
    })
    names(config) <- keys
  }
  set_config(config)
  cli::cli_alert_success("Configuration loaded from {.file {filename}}")
  invisible(config)
}

#' Reset the configuration to default
#' @return Invisibly `NULL`.
#' @export
#' @examples
#' reset_config()
reset_config <- function() {
  initialize_config()
  cli::cli_alert_success("HaploTraitR configuration reset to default.")
  invisible(NULL)
}
