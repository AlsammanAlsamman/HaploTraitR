library(HaploTraitR)

# Step 0: Reset configuration
HaploTraitR::reset_config()

# Step 0.1: List available configuration parameters
# list_config()

# Step 0.2: Set the phenotype information
phenotype <- "Area"
phenotype_unit <- "cm2"
set_config(config = list(phenotypename = phenotype,
                         pheno_col = phenotype,
                         phenotypeunit = phenotype_unit))

# Step 1: Create a unique output folder for results
result_folder <- create_unique_result_folder(location = "sampleout")
set_config(config = list(outfolder = result_folder))

# Step 2: Read Data Files from Sample Data
gwas_file <- "sampledata/gwas_area_ann19.csv"
hapmap_file <- "sampledata/Barley_50K.tsv"
pheno_file <- "sampledata/area_ann19.tsv"

# Read GWAS data
gwas <- read_gwas_file(gwas_file, sep = ",")
# Read hapmap data
hapmap <- readHapmap(hapmap_file)
# Read phenotype data
pheno <- read.csv(pheno_file, sep = "\t", header = TRUE)

# Or you can use the data provided in the package
# data("barley_area", package = "HaploTraitR")
# gwas <- barley_area$gwas
# hapmap <- barley_area$hapmap
# pheno <- barley_area$pheno

# Step 3: Set FDR threshold and filter GWAS data
set_config(config = list(fdr_threshold = 0.1))
gwas_filtered <- filter_gwas_data(gwas)
print(names(gwas_filtered))
print(gwas_filtered)

# Step 4: Generate Boxplots for Genotype-Phenotype Distribution
bxplts <- boxplot_genotype_phenotype(pheno, gwas_filtered, hapmap)

# Step 5: Cluster SNPs by Chromosome
dist_clusters <- getDistClusters(gwas_filtered, hapmap)

# Step 6: Compute LD Matrices
LDsInfo <- computeLDclusters(hapmap, dist_clusters)

# Step 7: Cluster LD Values and Convert to Haplotypes
clusterLDs <- getLDclusters(LDsInfo)
haplotypes <- HaploTraitR::convertLDclusters2Haps(hapmap, clusterLDs)

# Step 8: Get Sample Haplotypes and Save
snpcomb_sample <- HaploTraitR::getHapCombSamples(haplotypes, hapmap)

# Step 9: Generate SNP Combination Tables
SNPcombTables <- getSNPcombTables(snpcomb_sample, pheno)

# Step 10: Perform t-tests on SNP Combinations
t_test_snpComp <- HaploTraitR::testSNPcombs(SNPcombTables)

# Step 11: Generate and Save Haplotype Combination Boxplots
generateHapCombBoxPlots(SNPcombTables, t_test_snpComp)

# Step 12: Plot Haplotype Distribution (Optional)
# lDplots <- plotLDForClusters(clusterLDs)

# Step 13: Plot LD Combination Matrix with LD
LDCombs <- plotLDCombMatrix(clusterLDs, haplotypes, gwas_filtered)
# Note: For a large number of SNPs, the plot may not be clear or may raise an error.

# Step 14: Plot Haplotypes on the Genome
plot_haplotypes_genome(gwas, haplotypes)
