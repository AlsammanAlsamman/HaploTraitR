# Console output built on cli (always installed with ggplot2). All output is
# sent as messages, so it can be silenced with suppressMessages().

# Mid-grey that stays readable on light and dark consoles
muted <- function(x) cli::make_ansi_style("grey55")(x)

# Pad a possibly coloured string to a fixed display width
pad <- function(x, width, align = "left") cli::ansi_align(x, width = width, align = align)

# Colour a haplotype evidence level
style_evidence <- function(x) {
  ifelse(x == "Strong", cli::col_green(cli::style_bold(x)),
         ifelse(x == "Moderate", cli::col_yellow(x), muted(x)))
}

# Colour a percentage difference: green when favourable, red when unfavourable
style_effect <- function(pct, higher = isTRUE(get_config("higher_is_better"))) {
  txt <- ifelse(is.na(pct), "NA", sprintf("%+.1f%%", pct))
  good <- !is.na(pct) & ((pct > 0) == higher)
  ifelse(is.na(pct), muted(txt), ifelse(good, cli::col_green(txt), cli::col_red(txt)))
}

# Colour a p-value: significant values bold
style_p <- function(p, alpha = get_config("t_test_threshold")) {
  txt <- vapply(p, format_p, character(1))
  ifelse(!is.na(p) & p < alpha, cli::style_bold(txt), muted(txt))
}

# Print a table of block results with aligned, coloured columns
cli_block_table <- function(blocks) {
  cols <- list(
    Block = blocks$block,
    Haps = as.character(blocks$n_haplotypes),
    `p value` = style_p(blocks$p_value),
    Evidence = style_evidence(blocks$strength),
    Superior = ifelse(blocks$strength == "None", muted(blocks$superior),
                      cli::style_bold(blocks$superior)),
    `vs pop.` = style_effect(blocks$superior_vs_pop_pct)
  )
  widths <- vapply(names(cols), function(n) max(cli::ansi_nchar(c(n, cols[[n]]))), numeric(1)) + 2
  header <- paste(mapply(pad, names(cols), widths), collapse = "")
  cli::cli_verbatim(paste0("  ", cli::style_bold(header)))
  cli::cli_verbatim(paste0("  ", muted(strrep("\u2500", sum(widths) - 2))))
  for (i in seq_len(nrow(blocks))) {
    row <- paste(mapply(function(col, w) pad(col[i], w), cols, widths), collapse = "")
    cli::cli_verbatim(paste0("  ", row))
  }
}
