library(HaploTraitR)

# Example data shipped with the package
ex <- function(f) system.file("extdata", f, package = "HaploTraitR")

# Run the whole analysis in one call ------------------------------------------
results <- run_haplotraitr(
  gwas      = ex("gwas_area_ann19.csv.gz"),
  genotypes = ex("Barley_50K.tsv.gz"),
  pheno     = ex("area_ann19.tsv"),
  trait     = "Area",
  unit      = "cm2",
  outfolder = create_unique_result_folder(location = "sampleout"),
  higher_is_better = TRUE,
  config    = list(fdr_threshold = 0.1)
)

# Main tables (also in the Excel report)
results$summary$blocks     # one row per haplotype block, with recommendation
results$summary$effects    # one row per haplotype
results$markers            # diagnostic SNPs for the superior haplotypes
head(results$ranking)      # candidate accessions

# Re-draw a single block figure, e.g. for a presentation
fig <- plot_haplotype_block("2H:582673984", results$haplotypes, results$summary)
print(fig)
