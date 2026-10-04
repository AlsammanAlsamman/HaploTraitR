# Internal helpers -------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

# Make a string safe for use as a file name on every operating system
safe_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

# Return the configured output folder, creating it if needed
get_outfolder <- function(required = TRUE) {
  outfolder <- get_config("outfolder")
  if (is.null(outfolder)) {
    if (required) {
      stop("No output folder is set. Use set_config(list(outfolder = \"my_results\")) ",
           "or create_unique_result_folder().")
    }
    return(NULL)
  }
  if (!dir.exists(outfolder)) dir.create(outfolder, recursive = TRUE, showWarnings = FALSE)
  outfolder
}

# Return (and create) a sub folder of the output folder
out_subdir <- function(...) {
  path <- file.path(get_outfolder(), ...)
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

# Write a data frame as CSV into the 'tables' folder of the output folder
save_table <- function(df, filename) {
  outfolder <- get_outfolder(required = FALSE)
  if (is.null(outfolder)) return(invisible(NULL))
  path <- file.path(out_subdir("tables"), filename)
  utils::write.csv(df, path, row.names = FALSE)
  invisible(path)
}

# Trait label with unit, used for axis titles
trait_label <- function() {
  unit <- get_config("phenotypeunit")
  name <- get_config("phenotypename")
  if (is.null(unit) || !nzchar(unit)) name else paste0(name, " (", unit, ")")
}

# Save a ggplot or haplotraitr_figure using the configured format and dpi
save_plot <- function(plot, path_no_ext, width, height) {
  format <- tolower(get_config("plot_format"))
  dpi <- get_config("plot_dpi")
  file <- paste0(path_no_ext, ".", format)
  if (format == "pdf") {
    grDevices::pdf(file, width = width, height = height)
  } else {
    type <- if (capabilities("cairo")) "cairo" else NULL
    if (is.null(type)) {
      grDevices::png(file, width = width, height = height, units = "in", res = dpi)
    } else {
      grDevices::png(file, width = width, height = height, units = "in", res = dpi, type = type)
    }
  }
  on.exit(grDevices::dev.off())
  print(plot)
  invisible(file)
}

# Order chromosome names naturally (1H, 2H, ..., 10 after 9)
chr_order <- function(chr) {
  chr <- unique(as.character(chr))
  num <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", chr)))
  chr[order(is.na(num), num, chr)]
}

# Order haplotype labels: H1, H2, ..., Other
hap_levels <- function(haps) {
  haps <- unique(as.character(haps))
  num <- suppressWarnings(as.integer(sub("^H", "", haps)))
  haps[order(is.na(num), num, haps)]
}

# Significance stars for adjusted p-values
signif_stars <- function(p) {
  out <- as.character(cut(p, c(-Inf, 1e-4, 1e-3, 1e-2, 0.05, Inf),
                          labels = c("****", "***", "**", "*", "ns")))
  out[is.na(out)] <- "ns"
  out
}

# Format a p-value for plot subtitles
format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-4) formatC(p, format = "e", digits = 1) else formatC(p, format = "g", digits = 2)
}

# Compact letter display (insert-and-absorb algorithm).
# groups: group names ordered best first; sig: two-column matrix of
# significantly different pairs. Groups sharing a letter do not differ.
cld_letters <- function(groups, sig) {
  cols <- list(groups)
  if (length(sig) > 0) {
    sig <- matrix(sig, ncol = 2)
    for (i in seq_len(nrow(sig))) {
      a <- sig[i, 1]
      b <- sig[i, 2]
      new <- list()
      for (col in cols) {
        if (a %in% col && b %in% col) {
          new <- c(new, list(setdiff(col, a)), list(setdiff(col, b)))
        } else {
          new <- c(new, list(col))
        }
      }
      # absorb columns that are contained in another column
      keep <- vapply(seq_along(new), function(j) {
        !any(vapply(seq_along(new), function(m) {
          m != j && all(new[[j]] %in% new[[m]]) &&
            (length(new[[m]]) > length(new[[j]]) || m < j)
        }, logical(1)))
      }, logical(1))
      cols <- new[keep & lengths(new) > 0]
    }
  }
  first <- vapply(cols, function(col) min(match(col, groups)), numeric(1))
  cols <- cols[order(first)]
  out <- vapply(groups, function(g) {
    paste(letters[which(vapply(cols, function(col) g %in% col, logical(1)))], collapse = "")
  }, character(1))
  names(out) <- groups
  out
}
