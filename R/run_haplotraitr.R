#' Run the complete HaploTraitR analysis
#'
#' One call that reads the data, builds haplotype blocks around the
#' significant GWAS SNPs, tests the haplotype effects on the trait, and writes
#' all plots, CSV tables and a formatted Excel report to the output folder.
#'
#' Output folder layout:
#' * `HaploTraitR_<trait>.xlsx` - the Excel report (start here)
#' * `plots/genome_overview` - GWAS Manhattan plot with the haplotype blocks
#' * `plots/haplotype_effects_summary` - effect of every haplotype in every block
#' * `plots/candidate_accessions` - accessions carrying the most superior haplotypes
#' * `plots/haplotype_blocks/` - alleles, trait effects and LD of each block
#' * `plots/haplotype_trait_plots/` - trait distribution per haplotype
#' * `plots/gwas_snp_plots/` - trait distribution per genotype of each GWAS SNP
#' * `tables/` - CSV copies of all tables; `LD_matrices/` - LD (r2) matrices
#'
#' @param gwas GWAS results: a file path or a data frame from [read_gwas_file()]
#' @param genotypes Genotypes: a HapMap or VCF file path, or the list returned by
#'   [readHapmap()] / [readVCF()]
#' @param pheno Phenotypes: a file path or a data frame whose first column
#'   holds the accession names
#' @param trait Name of the trait column in `pheno` (default: the second column)
#' @param unit Trait unit shown on the plots (optional)
#' @param outfolder Output folder. Default: a new `haplotraitR_run` folder in the
#'   working directory
#' @param higher_is_better `TRUE` if higher trait values are desirable (yield),
#'   `FALSE` otherwise (disease score, days to heading...)
#' @param config Other configuration parameters as a named list, see [list_config()]
#' @param snp_plots Also draw one plot per significant GWAS SNP
#' @return Invisibly, a `haplotraitr_result` list with the elements
#'   `summary`, `haplotypes`, `accessions`, `pairwise`, `markers`, `ranking`,
#'   `gwas_significant`, `LDblocks` and `outfolder`.
#' @export
#' @examples
#' \donttest{
#' ex <- function(f) system.file("extdata", f, package = "HaploTraitR")
#' res <- run_haplotraitr(
#'   gwas = ex("gwas_area_ann19.csv.gz"),
#'   genotypes = ex("Barley_50K.tsv.gz"),
#'   pheno = ex("area_ann19.tsv"),
#'   trait = "Area", unit = "cm2",
#'   outfolder = file.path(tempdir(), "haplo_example"),
#'   config = list(fdr_threshold = 0.1)
#' )
#' res
#' }
run_haplotraitr <- function(gwas, genotypes, pheno, trait = NULL, unit = NULL, outfolder = NULL,
                            higher_is_better = TRUE, config = list(), snp_plots = TRUE) {
  t_start <- Sys.time()
  n_steps <- if (snp_plots) 8 else 7
  i_step <- 0
  step <- function(text) {
    i_step <<- i_step + 1
    counter <- cli::col_cyan(sprintf("[%d/%d]", i_step, n_steps))
    cli::cli_text("{counter} {.strong {text}}")
  }

  # settings
  if (length(config)) set_config(config)
  if (is.character(pheno)) pheno <- read_table_auto(pheno)
  pheno <- as.data.frame(pheno)
  if (is.null(trait)) trait <- colnames(pheno)[2]
  if (!trait %in% colnames(pheno)) {
    stop("Trait '", trait, "' not found. Phenotype columns: ", paste(colnames(pheno)[-1], collapse = ", "))
  }
  pheno <- pheno[, c(colnames(pheno)[1], trait)]
  set_config(list(phenotypename = trait, pheno_col = trait, phenotypeunit = unit,
                  higher_is_better = higher_is_better))
  if (is.null(outfolder)) outfolder <- create_unique_result_folder()
  set_config(list(outfolder = outfolder))
  get_outfolder()

  cli::cli_h1("HaploTraitR {utils::packageVersion('HaploTraitR')}")
  direction <- if (higher_is_better) "higher is better" else "lower is better"
  cli::cli_text("Trait: {.strong {trait_label()}} {cli::symbol$bullet} {direction} {cli::symbol$bullet} output: {.path {outfolder}}")
  cli::cli_text("")
  indent <- cli::cli_div(theme = list(".alert" = list("margin-left" = 6), ".bullets" = list("margin-left" = 6)))

  step("Reading data")
  if (is.character(gwas)) gwas <- read_gwas_file(gwas)
  if (is.character(genotypes)) {
    genotypes <- if (grepl("\\.vcf(\\.gz)?$", genotypes, ignore.case = TRUE)) readVCF(genotypes) else readHapmap(genotypes)
  }
  n_geno <- ncol(genotypes[[1]]) - 11
  big <- function(x) format(x, big.mark = ",")
  n_snps <- big(sum(vapply(genotypes, nrow, integer(1))))
  n_gwas <- big(nrow(gwas))
  n_pheno <- nrow(pheno)
  shared <- length(intersect(colnames(genotypes[[1]])[-(1:11)], pheno[[1]]))
  cli::cli_bullets(c(
    "*" = "GWAS: {.strong {n_gwas}} SNPs",
    "*" = "Genotypes: {.strong {n_snps}} SNPs {cli::symbol$times} {.strong {n_geno}} accessions",
    "*" = "Phenotypes: {.strong {n_pheno}} accessions ({shared} also genotyped)"))
  if (shared == 0) stop("No accession name is shared by the genotype and phenotype data.")

  step("Selecting significant GWAS SNPs")
  gwas_sig <- filter_gwas_data(gwas)
  if (length(gwas_sig) == 0) stop("No significant GWAS SNP. Increase fdr_threshold, e.g. config = list(fdr_threshold = 0.1).")

  step("Building LD blocks around the GWAS SNPs")
  windows <- getDistClusters(gwas_sig, genotypes)
  if (length(windows) == 0) stop("No GWAS SNP has enough genotyped SNPs nearby. Try a larger dist_threshold.")
  ld_info <- computeLDclusters(genotypes, windows)
  ld_blocks <- getLDclusters(ld_info)

  step("Identifying haplotypes")
  haps <- convertLDclusters2Haps(genotypes, ld_blocks)
  if (nrow(haps) == 0) stop("No haplotype block could be built. Try a lower ld_threshold or comb_freq_threshold.")
  haps <- getHapCombSamples(haps, genotypes)
  accessions <- getSNPcombTables(haps, pheno)

  step("Testing haplotype effects")
  pairwise <- testSNPcombs(accessions, savecopy = FALSE)
  summary <- summarize_haplotypes(haps, accessions, pairwise)
  markers <- diagnostic_markers(haps, summary)
  ranking <- rank_accessions(haps, pheno, summary)
  suppressMessages(saveTTestResultsToFile(pairwise))
  save_table(summary$blocks, "block_summary.csv")
  save_table(summary$effects, "haplotype_effects.csv")
  save_table(markers, "diagnostic_markers.csv")
  save_table(ranking, "candidate_accessions.csv")

  n_blocks <- nrow(summary$blocks)
  step("Drawing haplotype plots")
  plot_haplotypes_genome(gwas, haps, summary = summary)
  plot_effect_summary(summary)
  plot_candidate_accessions(ranking, summary)
  plotLDCombMatrix(haplotypes = haps, summary = summary)
  generateHapCombBoxPlots(accessions, pairwise)

  if (snp_plots) {
    step("Drawing GWAS SNP plots")
    boxplot_genotype_phenotype(pheno, gwas_sig, genotypes)
  }

  result <- structure(list(summary = summary, haplotypes = haps, accessions = accessions,
                           pairwise = pairwise, markers = markers, ranking = ranking,
                           gwas_significant = do.call(rbind, gwas_sig), LDblocks = ld_blocks,
                           outfolder = outfolder, trait = trait, unit = unit,
                           higher_is_better = higher_is_better,
                           excel = file.path(outfolder, paste0("HaploTraitR_", safe_name(trait), ".xlsx"))),
                      class = "haplotraitr_result")
  step("Writing Excel report")
  suppressMessages({
    write_excel_report(result, result$excel)
    save_config(file.path(outfolder, "parameters.txt"))
  })
  cli::cli_end(indent)
  elapsed <- format(round(difftime(Sys.time(), t_start, units = "secs"), 1))
  cli::cli_text("")
  cli::cli_alert_success("Analysis finished in {elapsed}")
  print(result)
  invisible(result)
}

#' @export
print.haplotraitr_result <- function(x, ...) {
  b <- x$summary$blocks
  old <- set_config(list(higher_is_better = x$higher_is_better %||% TRUE))
  on.exit(set_config(old))
  unit <- if (is.null(x$unit) || !nzchar(x$unit)) "" else paste0(" (", x$unit, ")")
  n_sig <- sum(b$strength != "None")

  cli::cli_h2("Results: {x$trait %||% get_config('phenotypename')}{unit}")
  cli::cli_alert_success("{nrow(b)} haplotype block{?s}, {.strong {n_sig}} with a significant haplotype effect")
  cli::cli_text("")
  cli_block_table(b)
  cli::cli_text("")

  sig <- b[b$strength != "None", ]
  for (i in seq_len(nrow(sig))) {
    block <- sig$block[i]
    rec <- sig$recommendation[i]
    cli::cli_alert("{.strong {block}}: {rec}")
  }
  other_best <- b[grepl("best group is 'Other'", b$recommendation), ]
  if (nrow(other_best)) {
    blocks <- other_best$block
    cli::cli_alert_warning("In {.val {blocks}} the best group is {.emph Other} (rare haplotypes); lower {.field comb_freq_threshold} to split it.")
  }

  m <- x$markers
  if (!is.null(m) && nrow(m)) {
    n_m <- nrow(m)
    n_full <- sum(m$fully_diagnostic)
    cli::cli_alert_info("{n_m} diagnostic SNP{?s} for the superior haplotypes ({.strong {n_full}} fully diagnostic, ready for KASP design)")
  }
  if (n_sig > 0 && !is.null(x$ranking)) {
    top <- head(x$ranking$accession[x$ranking$n_superior > 0], 5)
    cli::cli_alert_info("Top candidate accessions: {.val {top}}")
  }

  if (!is.null(x$outfolder)) {
    excel <- x$excel %||% file.path(x$outfolder, "HaploTraitR.xlsx")
    plots <- file.path(x$outfolder, "plots")
    n_plots <- length(list.files(plots, recursive = TRUE))
    tables <- file.path(x$outfolder, "tables")
    cli::cli_h2("Output")
    cli::cli_bullets(c(
      "*" = "Excel report: {.file {excel}}",
      "*" = "Plots ({n_plots}): {.path {plots}}",
      "*" = "CSV tables: {.path {tables}}"))
  }
  invisible(x)
}
