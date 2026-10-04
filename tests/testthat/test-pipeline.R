test_that("the full pipeline runs on the barley example", {
  skip_on_cran()
  ex <- function(f) system.file("extdata", f, package = "HaploTraitR")
  out <- file.path(tempdir(), "haplotraitr_test")
  unlink(out, recursive = TRUE)
  on.exit(reset_config())
  res <- suppressMessages(run_haplotraitr(
    gwas = ex("gwas_area_ann19.csv.gz"), genotypes = ex("Barley_50K.tsv.gz"),
    pheno = ex("area_ann19.tsv"), trait = "Area", unit = "cm2", outfolder = out,
    config = list(fdr_threshold = 0.1), snp_plots = FALSE))

  expect_s3_class(res, "haplotraitr_result")
  # the five GWAS SNPs collapse into three distinct haplotype blocks
  expect_equal(nrow(res$summary$blocks), 3)
  strong <- res$summary$blocks[res$summary$blocks$strength == "Strong", ]
  expect_equal(strong$superior, "H3")
  # every accession sample is phenotyped exactly once per block
  expect_false(any(duplicated(res$accessions[, c("Sample", "block")])))
  # phenotype values are attached to the right accession
  pheno <- read.delim(ex("area_ann19.tsv"))
  expect_equal(res$accessions$Pheno, pheno$Area[match(res$accessions$Sample, pheno$Taxa)])
  expect_true(file.exists(file.path(out, "HaploTraitR_Area.xlsx")))
  expect_true(file.exists(file.path(out, "plots", "genome_overview.png")))
  expect_length(list.files(file.path(out, "plots", "haplotype_blocks")), 3)
})
