library(HaploTraitR)
setwd("/home/samman/Documents/MyGitHub/HaploTraitR")

# Define file paths and output folder
hapfile <- "sampledata/Barley_50K_KNNimp.hmp.txt"
gwasfile <- "sampledata/SignificantSNP_GWAS.tsv"
phenofile <- "sampledata/Pheno_ANN19.tsv"
outfolder <- "sampleout"

### Setting Thresholds

dist_threshold <- 1000000 # Distance threshold (1Mb)
dist_cluster_count <- 5   # Minimum SNPs in a cluster
ld_threshold <- 0.3       # Minimum LD value for two SNPs to be in LD
comb_freq_threshold <- 0.1 # Minimum genotype frequency for a combination


## Step 1: Reading Data Files

gwas <- readGWAS(gwasfile, sep = "\t")

hapmap <- readHapmap(hapfile)

pheno <- read.csv(phenofile, header = TRUE, sep = "\t")

## Step 2: Looking at the the variation phenotypic data across the genotypes
# The analysis uses t-tests to compare the phenotypic data across the genotypes

# select significant snps
subhapmap<-extract_hapmap(hapmap, gwas)
# get pheno and geno
geno_pheno_table<-get_pheno_geno(subhapmap, pheno)

# plot the boxplot
# create a subfolder to store the plots
boxplot_geno_pheno_folder<-file.path(outfolder, "boxplot_geno_pheno")
if (!dir.exists(boxplot_geno_pheno_folder)) {
  dir.create(boxplot_geno_pheno_folder, showWarnings = FALSE)
}


boxplot_genotype_phenotype(genotype_phenotype_data = geno_pheno_table,
                           outfolder = boxplot_geno_pheno_folder,
                           method = "t.test")

############################################################################################################
## Step 2: Cluster SNPs by Chromosome

dist_clusters <- getDistClusters(gwas, hapmap, dist_threshold, dist_cluster_count)

## Step 3: Compute LD Matrices

LDsInfo <- computeLDclusters(hapmap, dist_clusters, outfolder)

## Step 4: Save LD Matrices

saveLDs2folder(LDsInfo, outfolder)

## Step 5: Cluster LD Values and Convert to Haplotypes

clusterLDs <- getLDclusters(LDsInfo,ld_threshold, cls_count = 3)

haplotypes <- HaploTraitR::convertLDclusters2Haps(hapmap, clusterLDs, comb_freq_threshold)

## Step 6: Get Sample Haplotypes and Save

snpcombsample <- HaploTraitR::getHapCombSamples(haplotypes, hapmap)

snpcombsample <- as.data.frame(snpcombsample)
write.csv(snpcombsample, file = file.path(outfolder, "haplotype_combinations.csv"), row.names = FALSE)

## Step 8: Generate SNP Combination Tables

SNPcombTables <- getSNPcombTables(snpcombsample, pheno)

## Step 9: Perform t-tests on SNP Combinations


t_test_snpComp <-  HaploTraitR::testSNPcombs(SNPcombTables)

## Step 10: Visualize Results


for (cls_snp in names(t_test_snpComp)) {
       HaploTraitR::plotHapCombBoxPlot(cls_snp, SNPcombTables, t_test_snpComp, outfolder= outfolder)
       ggplot2::ggsave(file.path(outfolder, paste0(cls_snp, "_LD_boxplot.png")), width = 8, height = 6)
}

### Plot Haplotype Distribution
#library(ggplot2)
lDplots<-plotLDForClusters(ld_matrix_folder=LDsInfo[[1]], clusterLDs, outfolder=outfolder)

### Plot LD Combination Matrix with LD

LDCombs<-plotLDCombMatrix(LDsInfo[[1]], clusterLDs, haplotypes, gwas, outfolder=outfolder)

