#!/usr/bin/env Rscript

# outfolder: NULL (Output folder for results)
# dist_threshold: 1000000 (Distance threshold for clustering SNPs within 1Mb)
# dist_cluster_count: 5 (Minimum number of SNPs in a cluster)
# ld_threshold: 0.3 (Minimum LD value to consider two SNPs in LD)
# comb_freq_threshold: 0.1 (Minimum frequency of genotypes with the same combination)
# fdr_threshold: 0.05 (False discovery rate threshold for multiple testing correction)
# fdr_method: "fdr" (Method for FDR correction)
# t_test_threshold: 0.05 (Threshold for t-test significance)
# t_test_method: "t.test" (Method for t-test)
# rsid_col: "rsid" (Column name for rsid in GWAS file)
# pos_col: "pos" (Column name for position in GWAS file)
# chr_col: "chr" (Column name for chromosome in GWAS file)
# pval_col: "p" (Column name for p-value in GWAS file)
# fdr_col: "fdr" (Column name for FDR in GWAS file)
# phenotypename: "Phenotype" (Name of the phenotype)
# phenotypeunit: NULL (Unit of the phenotype)

# Load required libraries
library(HaploTraitR)


setwd("~/Documents/publishing/HaploTraitR/Testing/Barley")
# result folder
result_folder<-"~/Documents/publishing/HaploTraitR/Testing/Barley/results"


# gwas folder contains several gwas data
gwas_folder<-"~/Documents/publishing/HaploTraitR/Testing/Barley/gwas"
set_config(config = list( rsid_col = "Marker",
                          pos_col = "Pos",
                          chr_col = "Chr",
                          pval_col = "p"))

# phenotypes folder contains several phenotype data
pheno_folder<-"~/Documents/publishing/HaploTraitR/Testing/Barley/pheno"
# hapmap file contains the genotype data
hapmap_file<-"~/Documents/publishing/HaploTraitR/Testing/Barley/Barley_50K.txt"

# loop over the gwas files
gwas_files<-list.files(gwas_folder, full.names = TRUE)
phenotype_files<-list.files(pheno_folder, full.names = TRUE)

# read hapmap data
hapmap <- readHapmap(hapmap_file)

for(gw in gwas_files)
{
  #gw<-"~/Documents/publishing/HaploTraitR/Testing/Barley/gwas/Area_ANN19.csv"
  gw_basename<-basename(gw)
  gw_basename
  # remove .csv extension
  gw_basename<-substr(gw_basename,1,nchar(gw_basename)-4)
  trial<-strsplit(gw_basename,"_")[[1]][2]
  #remove csv extension
  phenotyepe<-strsplit(gw_basename,"_")[[1]][1]
  pheno_file<-paste0(pheno_folder,"/",trial,".tsv")

  # Read GWAS data
  gwas <- read_gwas_file(gw, sep = ",")
  # filter GWAS data by FDR threshold
  gwas_filtered <- filter_gwas_data(gwas)

  # if the list is empty, skip the rest of the loop
  if(length(gwas_filtered)==0)
  {
    next
  }
  # Read Phenotype data
  set_config(config = list(phenotypename = phenotyepe, phenotypeunit = NULL, pheno_col = phenotyepe))
  pheno <- read.csv(pheno_file, header = TRUE, sep = "\t")

  pheno <-pheno[,c("Taxa",get_config("pheno_col"))]

  # create a folder for the results for each gwas_trial
  gw_trial_result_folder <-paste0(result_folder,"/",gw_basename)
  if (!dir.exists(gw_trial_result_folder)) {
    dir.create(gw_trial_result_folder, showWarnings = FALSE, recursive = TRUE)
  }else{
    # if the folder already exists, skip the rest of the loop
    next
  }
  set_config(config = list(outfolder = gw_trial_result_folder))
  # Step 4: Cluster SNPs by Chromosome
  # This step clusters SNPs based on their chromosome and physical distance.
  dist_clusters <- getDistClusters(gwas_filtered, hapmap)
  # Step 5: Compute LD Matrices
  # This step computes linkage disequilibrium (LD) matrices for each SNP cluster.
  # The matrices will be saved in the outfolder/LD_matrices folder.
  LDsInfo <- computeLDclusters(hapmap, dist_clusters)
  # if no clusters are found, skip the rest of the loop
  if(length(LDsInfo$out_info)==0)
  {
    next
  }
  # Step 6: Cluster LD Values and Convert to Haplotypes
  # This step clusters LD values and converts them to haplotypes.
  clusterLDs <- getLDclusters(LDsInfo)

  haplotypes <- HaploTraitR::convertLDclusters2Haps(hapmap, clusterLDs)
  # if no haplotypes
  if(nrow(haplotypes)==0)
  {
    next
  }
  # Step 7: Get Sample Haplotypes and Save
  # This step extracts sample haplotypes from the hapmap data.
  snpcombsample <- HaploTraitR::getHapCombSamples(haplotypes, hapmap)
  # Step 8: Generate SNP Combination Tables
  # This step generates tables of SNP combinations and their corresponding phenotypes.
  SNPcombTables <- getSNPcombTables(snpcombsample, pheno)
  # Step 9: Perform t-tests on SNP Combinations
  # This step performs t-tests on the SNP combinations to identify significant differences.
  t_test_snpComp <-  HaploTraitR::testSNPcombs(SNPcombTables)
  # if no t-test is significant, skip the rest of the loop
  if(length(t_test_snpComp)==0)
  {
    next
  }
  # Step 10: Generate and Save Haplotype Combination Boxplots
  # This step generates boxplots showing the distribution of haplotype combinations and saves them in categorized directories.
  generateHapCombBoxPlots(SNPcombTables, t_test_snpComp)
  # Step 11: Plot Haplotype Distribution
  # This step plots the distribution of haplotypes for each LD cluster.
  #lDplots <- plotLDForClusters(clusterLDs)
  # Step 12: Plot LD Combination Matrix with LD
  # This step plots the LD combination matrix, showing the LD values for different SNP combinations.
  tryCatch({
  LDCombs <- plotLDCombMatrix(clusterLDs, haplotypes, gwas_filtered)
  }, error = function(e) {
    message("Error in plotting LD Combination Matrix with LD")
  })
  # Step 13: Generate Boxplots for Genotype-Phenotype Distribution
  plot_haplotypes_genome(gwas, haplotypes)
}






# Step 1: Create a unique output folder for results
#result_folder <- create_unique_result_folder(location = "/home/samman/Documents/publishing/HaploTraitR/Testing/Barley/results")

# initialize the configuration
#initialize_config()
result_folder <- create_unique_result_folder(location = "~/Documents/publishing/HaploTraitR/Testing/Barley/results/haplotraitR_run")
set_config(config = list(outfolder = result_folder))



# Step 2: Read Data Files

# GWAS data file containing significant SNPs
gwasfile <- "/home/samman/Documents/publishing/HaploTraitR/Testing/Barley/gwas/Area_ANN19.csv"

# Hapmap file containing SNP genotypes
hapfile <- "/home/samman/Documents/publishing/HaploTraitR/Testing/Barley/Barley_50K.txt"

# Read Hapmap data
hapmap <- readHapmap(hapfile)
# Read GWAS data
gwas <- read_gwas_file(gwasfile, sep = ",")


# filter GWAS data by FDR threshold
gwas_filtered <- filter_gwas_data(gwas)


# Phenotype data file
phenofile <- "/home/samman/Documents/publishing/HaploTraitR/Testing/Barley/pheno/ANN19.tsv"
# Read Phenotype data
pheno <- read.csv(phenofile, header = TRUE, sep = "\t")


set_config(config = list(phenotypename = "Area", phenotypeunit = NULL, pheno_col = "Area"))
pheno <-pheno[,c("Taxa",get_config("pheno_col"))]




# Step 4: Cluster SNPs by Chromosome
# This step clusters SNPs based on their chromosome and physical distance.
dist_clusters <- getDistClusters(gwas_filtered, hapmap)

# Step 5: Compute LD Matrices
# This step computes linkage disequilibrium (LD) matrices for each SNP cluster.
# The matrices will be saved in the outfolder/LD_matrices folder.
LDsInfo <- computeLDclusters(hapmap, dist_clusters)
#LDsInfo <- retrieveLDMatricesFromFolder("~/Documents/publishing/HaploTraitR/Testing/Barley/results/haplotraitR_run/haplotraitR_run_2/LD_matrices")

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




############# Plots
# Step 3: Generate Boxplots for Genotype-Phenotype Distribution
# This step generates boxplots showing the distribution of haplotype combinations compared to each other.
bxplts <- boxplot_genotype_phenotype(pheno, gwas_filtered, hapmap)


plot_haplotypes_genome(gwas, haplotypes)
# Step 10: Generate and Save Haplotype Combination Boxplots
# This step generates boxplots showing the distribution of haplotype combinations and saves them in categorized directories.
# somthing wrong here
generateHapCombBoxPlots(SNPcombTables, t_test_snpComp)
# Step 11: Plot Haplotype Distribution
# This step plots the distribution of haplotypes for each LD cluster.
#lDplots <- plotLDForClusters(clusterLDs)

# Step 12: Plot LD Combination Matrix with LD
# This step plots the LD combination matrix, showing the LD values for different SNP combinations.
LDCombs <- plotLDCombMatrix(clusterLDs, haplotypes, gwas_filtered)
