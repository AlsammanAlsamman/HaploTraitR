#' Test trait differences between haplotypes
#'
#' For every haplotype block, all pairs of haplotypes (with at least two
#' phenotyped accessions) are compared with a Welch t-test or a Wilcoxon test
#' (`t_test_method`), and p-values are adjusted with `p_adjust_method`.
#' @param comb_sample_Tables The accession table (from [getSNPcombTables()])
#' @param savecopy Save the results (`tables/<trait>_pairwise_tests.csv`) in the output folder
#' @return A named list (per block) of data frames with the columns `group1`,
#'   `group2`, `n1`, `n2`, `mean1`, `mean2`, `diff`, `statistic`, `df`, `p`,
#'   `p.adj`, `p.adj.signif`, `block`, `y.position`, `xmin` and `xmax`.
#'   `NULL` if no block could be tested.
#' @export
testSNPcombs <- function(comb_sample_Tables, savecopy = TRUE) {
  method <- get_config("t_test_method")
  adjust <- get_config("p_adjust_method")
  results <- list()
  for (block in unique(comb_sample_Tables$block)) {
    d <- comb_sample_Tables[comb_sample_Tables$block == block, ]
    levels <- hap_levels(d$hap)
    counts <- table(factor(d$hap, levels))
    groups <- names(counts)[counts >= 2]
    if (length(groups) < 2) next
    pairs <- utils::combn(groups, 2)
    res <- do.call(rbind, lapply(seq_len(ncol(pairs)), function(i) {
      x <- d$Pheno[d$hap == pairs[1, i]]
      y <- d$Pheno[d$hap == pairs[2, i]]
      test <- tryCatch(
        if (method == "wilcox.test") stats::wilcox.test(x, y, exact = FALSE) else stats::t.test(x, y),
        error = function(e) NULL
      )
      data.frame(group1 = pairs[1, i], group2 = pairs[2, i], n1 = length(x), n2 = length(y),
                 mean1 = mean(x), mean2 = mean(y), diff = mean(x) - mean(y),
                 statistic = if (is.null(test)) NA else unname(test$statistic),
                 df = if (is.null(test$parameter)) NA else unname(test$parameter),
                 p = if (is.null(test)) NA else test$p.value,
                 stringsAsFactors = FALSE)
    }))
    res$p.adj <- p.adjust(res$p, method = adjust)
    res$p.adj.signif <- signif_stars(res$p.adj)
    res$block <- block
    step <- 0.08 * diff(range(d$Pheno))
    res$y.position <- max(d$Pheno) + step * seq_len(nrow(res))
    res$xmin <- match(res$group1, levels)
    res$xmax <- match(res$group2, levels)
    results[[block]] <- res
  }
  if (length(results) == 0) {
    warning("No haplotype block had at least two haplotypes with phenotyped accessions.")
    return(NULL)
  }
  if (savecopy) saveTTestResultsToFile(results)
  results
}

#' Save the haplotype test results to a file
#'
#' @param t_test_snpComp The list returned by [testSNPcombs()]
#' @return Invisibly, the path of the CSV file.
#' @export
saveTTestResultsToFile <- function(t_test_snpComp) {
  df <- do.call(rbind, t_test_snpComp)
  rownames(df) <- NULL
  df$trait <- get_config("phenotypename")
  path <- save_table(df[, c("block", "group1", "group2", "n1", "n2", "mean1", "mean2", "diff",
                            "statistic", "df", "p", "p.adj", "p.adj.signif", "trait")],
                     paste0(safe_name(get_config("phenotypename")), "_pairwise_tests.csv"))
  if (!is.null(path)) cli::cli_alert_success("Haplotype tests saved to {.file {path}}")
  invisible(path)
}

#' Summarise haplotype effects for breeders
#'
#' Computes, for every block and haplotype, the trait mean with its 95%
#' confidence interval, the difference from the population mean, letter groups
#' (haplotypes sharing a letter do not differ significantly) and identifies the
#' superior haplotype according to `higher_is_better`.
#' @param haplotypes A haplotype table with samples (from [getHapCombSamples()])
#' @param SNPcombTables The accession table (from [getSNPcombTables()])
#' @param t_test_snpComp Pairwise tests (from [testSNPcombs()]); computed if `NULL`
#' @return A list with two data frames: `effects` (one row per haplotype) and
#'   `blocks` (one row per block, including a selection recommendation).
#' @export
summarize_haplotypes <- function(haplotypes, SNPcombTables, t_test_snpComp = NULL) {
  if (is.null(t_test_snpComp)) t_test_snpComp <- testSNPcombs(SNPcombTables, savecopy = FALSE)
  alpha <- get_config("t_test_threshold")
  higher <- isTRUE(get_config("higher_is_better"))
  method <- get_config("t_test_method")
  unit <- get_config("phenotypeunit") %||% ""

  effects <- list()
  blocks <- list()
  for (block in unique(haplotypes$block)) {
    h <- haplotypes[haplotypes$block == block, ]
    d <- SNPcombTables[SNPcombTables$block == block, ]
    if (nrow(d) == 0) next
    pop_mean <- mean(d$Pheno)
    groups <- hap_levels(h$hap)

    stats <- do.call(rbind, lapply(groups, function(g) {
      v <- d$Pheno[d$hap == g]
      n <- length(v)
      m <- if (n) mean(v) else NA
      se <- if (n > 1) stats::sd(v) / sqrt(n) else NA
      half <- if (n > 1) stats::qt(0.975, n - 1) * se else NA
      data.frame(hap = g, n_pheno = n, mean = m, sd = if (n > 1) stats::sd(v) else NA,
                 se = se, ci_low = m - half, ci_high = m + half,
                 median = if (n) stats::median(v) else NA,
                 min = if (n) min(v) else NA, max = if (n) max(v) else NA,
                 stringsAsFactors = FALSE)
    }))
    stats <- merge(h[, c("block", "chr", "start", "end", "n_snps", "lead_snps", "hap", "alleles", "n", "freq")],
                   stats, by = "hap", sort = FALSE)
    stats <- stats[match(groups, stats$hap), ]
    stats$diff_vs_pop <- stats$mean - pop_mean
    stats$pct_vs_pop <- 100 * stats$diff_vs_pop / pop_mean

    # overall test across haplotypes with >= 2 phenotyped accessions
    testable <- stats$hap[stats$n_pheno >= 2]
    dt <- d[d$hap %in% testable, ]
    block_p <- NA
    if (length(testable) >= 2) {
      block_p <- tryCatch(
        if (method == "wilcox.test") stats::kruskal.test(Pheno ~ factor(hap), data = dt)$p.value
        else stats::oneway.test(Pheno ~ factor(hap), data = dt, var.equal = FALSE)$p.value,
        error = function(e) NA)
    }

    # letter groups, best haplotype first
    ord <- order(if (higher) -stats$mean else stats$mean, na.last = TRUE)
    ranked <- stats$hap[ord][stats$hap[ord] %in% testable]
    pw <- t_test_snpComp[[block]]
    sig <- if (is.null(pw)) NULL else as.matrix(pw[!is.na(pw$p.adj) & pw$p.adj < alpha, c("group1", "group2")])
    stats$letters <- ""
    if (length(ranked)) stats$letters[match(ranked, stats$hap)] <- cld_letters(ranked, sig)
    stats$rank <- match(stats$hap, stats$hap[ord])

    # superior haplotype: best named haplotype with >= 2 phenotyped accessions
    named <- stats[stats$hap != "Other" & stats$n_pheno >= 2, ]
    named <- named[order(if (higher) -named$mean else named$mean), ]
    superior <- if (nrow(named)) named$hap[1] else NA
    stats$superior <- stats$hap %in% superior
    effects[[block]] <- stats

    significant <- !is.na(block_p) && block_p < alpha
    n_differs <- 0
    if (!is.na(superior) && !is.null(pw)) {
      others <- setdiff(named$hap, superior)
      hits <- pw[(pw$group1 == superior & pw$group2 %in% others) |
                   (pw$group2 == superior & pw$group1 %in% others), ]
      n_differs <- sum(hits$p.adj < alpha, na.rm = TRUE)
    }
    strength <- if (!significant || n_differs == 0) "None"
                else if (n_differs == nrow(named) - 1) "Strong" else "Moderate"
    sup <- stats[stats$hap %in% superior, ]
    worst <- if (nrow(named)) named[nrow(named), ] else NULL
    recommendation <- if (is.na(superior)) {
      "Not enough phenotyped accessions to compare haplotypes."
    } else if (significant && length(ranked) && ranked[1] == "Other") {
      other <- stats[stats$hap == "Other", ]
      sprintf("Significant effect (p = %s), but the best group is 'Other' (rare haplotypes, mean %.2f %s, %d accessions). Lower comb_freq_threshold to split it into separate haplotypes.",
              format_p(block_p), other$mean, unit, other$n)
    } else if (strength == "None") {
      sprintf("No significant haplotype effect (p = %s); low priority for selection.", format_p(block_p))
    } else {
      sprintf("Select %s: mean %.2f %s (%+.1f%% vs population); significantly %s than %d of %d other haplotype(s); carried by %d accessions (%.0f%%).",
              superior, sup$mean, unit, sup$pct_vs_pop, if (higher) "higher" else "lower",
              n_differs, nrow(named) - 1, sup$n, 100 * sup$freq)
    }
    blocks[[block]] <- data.frame(
      block = block, chr = h$chr[1], start = h$start[1], end = h$end[1], n_snps = h$n_snps[1],
      lead_snps = h$lead_snps[1], n_haplotypes = sum(h$hap != "Other"),
      n_accessions = nrow(d), test = if (method == "wilcox.test") "Kruskal-Wallis" else "Welch ANOVA",
      p_value = block_p, significant = significant, superior = superior,
      superior_mean = if (nrow(sup)) sup$mean else NA, population_mean = pop_mean,
      superior_vs_pop_pct = if (nrow(sup)) sup$pct_vs_pop else NA,
      best_worst_diff = if (!is.null(worst) && nrow(sup)) sup$mean - worst$mean else NA,
      superior_n = if (nrow(sup)) sup$n else NA, strength = strength,
      recommendation = recommendation, stringsAsFactors = FALSE)
  }
  effects <- do.call(rbind, effects)
  blocks <- do.call(rbind, blocks)
  rownames(effects) <- NULL
  rownames(blocks) <- NULL
  if (!is.null(blocks)) blocks$fdr <- p.adjust(blocks$p_value, method = "BH")
  list(effects = effects, blocks = blocks)
}

#' Find diagnostic SNPs for the superior haplotypes
#'
#' For each block with a significant effect, lists the SNPs whose allele in
#' the superior haplotype differs from the other haplotypes. Fully diagnostic
#' SNPs distinguish the superior haplotype from all others and are good
#' candidates for marker-assisted selection (e.g. KASP assays).
#' @param haplotypes A haplotype table (from [getHapCombSamples()])
#' @param summary The list returned by [summarize_haplotypes()]
#' @param significant_only Only report blocks with a significant effect
#' @return A data frame with one row per informative SNP.
#' @export
diagnostic_markers <- function(haplotypes, summary, significant_only = TRUE) {
  blocks <- summary$blocks
  if (significant_only) blocks <- blocks[blocks$strength != "None", ]
  out <- list()
  for (i in seq_len(nrow(blocks))) {
    block <- blocks$block[i]
    sup <- blocks$superior[i]
    h <- haplotypes[haplotypes$block == block & haplotypes$hap != "Other", ]
    if (is.na(sup) || nrow(h) < 2) next
    info <- block_info(haplotypes, block)
    al <- do.call(rbind, strsplit(h$alleles, "|", fixed = TRUE))
    rownames(al) <- h$hap
    others <- al[rownames(al) != sup, , drop = FALSE]
    fav <- al[sup, ]
    differs <- colSums(others != matrix(fav, nrow(others), length(fav), byrow = TRUE))
    keep <- differs > 0
    if (!any(keep)) next
    ld <- attr(info, "ld")
    leads <- intersect(attr(info, "leads"), info$snp)
    r2_lead <- if (length(leads)) apply(ld[info$snp, leads, drop = FALSE], 1, max, na.rm = TRUE) else NA
    other_alleles <- apply(others, 2, function(a) paste(unique(a), collapse = "/"))
    out[[block]] <- data.frame(
      block = block, superior = sup, marker = info$marker, snp = info$snp, pos = info$pos,
      favourable_allele = fav, other_alleles = other_alleles,
      distinguishes = paste0(differs, "/", nrow(others)),
      fully_diagnostic = differs == nrow(others),
      is_gwas_snp = info$snp %in% leads, r2_with_gwas_snp = round(r2_lead, 3),
      stringsAsFactors = FALSE)[keep, ]
  }
  out <- do.call(rbind, out)
  if (is.null(out)) return(data.frame())
  out <- out[order(out$block, -out$fully_diagnostic, -out$r2_with_gwas_snp), ]
  rownames(out) <- NULL
  out
}

#' Rank accessions by the superior haplotypes they carry
#'
#' Builds an accession x block matrix of haplotypes and counts, for every
#' accession, how many superior haplotypes of significant blocks it carries.
#' Accessions that stack many superior haplotypes are candidate parents.
#' Genotyped accessions without phenotype data are included.
#' @param haplotypes A haplotype table with samples (from [getHapCombSamples()])
#' @param pheno A phenotype data frame (accession, trait)
#' @param summary The list returned by [summarize_haplotypes()]
#' @return A data frame with one row per accession.
#' @export
rank_accessions <- function(haplotypes, pheno, summary) {
  pheno <- prepare_pheno(pheno)
  higher <- isTRUE(get_config("higher_is_better"))
  all_samples <- unique(unlist(strsplit(haplotypes$samples, "|", fixed = TRUE)))
  out <- data.frame(accession = all_samples, stringsAsFactors = FALSE)
  out$trait <- pheno$Pheno[match(all_samples, pheno$Sample)]
  blocks <- summary$blocks
  sig_blocks <- blocks$block[blocks$strength != "None"]
  out$n_superior <- 0L
  for (block in blocks$block) {
    h <- haplotypes[haplotypes$block == block, ]
    hap <- rep(NA_character_, length(all_samples))
    for (i in seq_len(nrow(h))) {
      s <- strsplit(h$samples[i], "|", fixed = TRUE)[[1]]
      hap[all_samples %in% s] <- h$hap[i]
    }
    out[[block]] <- hap
    if (block %in% sig_blocks) {
      sup <- blocks$superior[blocks$block == block]
      out$n_superior <- out$n_superior + (!is.na(hap) & hap == sup)
    }
  }
  out$n_significant_blocks <- length(sig_blocks)
  trait_order <- if (higher) -out$trait else out$trait
  out <- out[order(-out$n_superior, trait_order, na.last = TRUE), ]
  out$rank <- seq_len(nrow(out))
  rownames(out) <- NULL
  out[, c("rank", "accession", "trait", "n_superior", "n_significant_blocks", blocks$block)]
}
