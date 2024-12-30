library(HaploTraitR)

# This script performs a haplotype combination analysis using GWAS, hapmap, and phenotype data.
# The analysis consists of two parts: one that requires phenotype data and one that does not.

# Initialize configuration (uncomment if needed)
# initialize_config()

# Step 1: Create a unique output folder for results
result_folder <- create_unique_result_folder(location = "sampleout")
set_config(config = list(outfolder = result_folder))
get_config("outfolder")

# Step 2: Read Data Files

# GWAS data file containing significant SNPs
gwasfile <- "sampledata/SignificantSNP_GWAS.tsv"

# Hapmap file containing SNP genotypes
hapfile <- "sampledata/Barley_50K_KNNimp.hmp.txt"

# Phenotype data file
phenofile <- "sampledata/Pheno_ANN19.tsv"

# Read GWAS data
gwas <- readGWAS(gwasfile, sep = "\t")

# Read Hapmap data
hapmap <- readHapmap(hapfile)

# Read Phenotype data
pheno <- read.csv(phenofile, header = TRUE, sep = "\t")

# Set phenotype information in the configuration
set_config(config = list(phenotypename = "PH", phenotypeunit = "cm"))

# Step 3: Generate Boxplots for Genotype-Phenotype Distribution
# This step generates boxplots showing the distribution of haplotype combinations compared to each other.
bxplts <- boxplot_genotype_phenotype(pheno, gwas, hapmap)

# Step 4: Cluster SNPs by Chromosome
# This step clusters SNPs based on their chromosome and physical distance.
dist_clusters <- getDistClusters(gwas, hapmap)

# Step 5: Compute LD Matrices
# This step computes linkage disequilibrium (LD) matrices for each SNP cluster.
# The matrices will be saved in the outfolder/LD_matrices folder.
LDsInfo <- computeLDclusters(hapmap, dist_clusters)

# Step 6: Cluster LD Values and Convert to Haplotypes
# This step clusters LD values and converts them to haplotypes.
clusterLDs <- getLDclusters(LDsInfo)
haplotypes <- HaploTraitR::convertLDclusters2Haps(hapmap, clusterLDs)

# Step 7: Get Sample Haplotypes and Save
# This step extracts sample haplotypes from the hapmap data.
snpcombsample <- HaploTraitR::getHapCombSamples(haplotypes, hapmap)

# Step 8: Generate SNP Combination Tables
# This step generates tables of SNP combinations and their corresponding phenotypes.
SNPcombTables <- getSNPcombTables(snpcombsample, pheno)

# Step 9: Perform t-tests on SNP Combinations
# This step performs t-tests on the SNP combinations to identify significant differences.
t_test_snpComp <- HaploTraitR::testSNPcombs(SNPcombTables)

# Step 10: Generate and Save Haplotype Combination Boxplots
# This step generates boxplots showing the distribution of haplotype combinations and saves them in categorized directories.
generateHapCombBoxPlots(SNPcombTables, t_test_snpComp)

# Step 11: Plot Haplotype Distribution
# This step plots the distribution of haplotypes for each LD cluster.
lDplots <- plotLDForClusters(clusterLDs)

# Step 12: Plot LD Combination Matrix with LD
# This step plots the LD combination matrix, showing the LD values for different SNP combinations.
LDCombs <- plotLDCombMatrix(clusterLDs, haplotypes, gwas)
