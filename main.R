
library(HaploTraitR)
library(ggplot2)
library(ggpubr)  # For statistical annotations
library(rstatix)
library(reshape2)

# Define file paths and output folder
hapfile <- "sampledata/Barley_50K_KNNimp.hmp.txt"
gwasfile <- "sampledata/SignificantSNP_GWAS.tsv"
phenofile <- "sampledata/Pheno_ANN19.tsv"
outfolder <- "sampleout"

### Setting Thresholds

dist_threshold <- 1000000 # Distance threshold (1Mb)
dist_cluster_count <- 5   # Minimum SNPs in a cluster
ld_threshold <- 0.3       # Minimum LD value
comb_freq_threshold <- 0.1 # Minimum genotype frequency for a combination


## Step 1: Reading Data Files

gwas <- readGWAS(gwasfile, sep = "\t")


hapmap <- readHapmap(hapfile)


## New addition
pheno <- read.csv(phenofile, header = TRUE, sep = "\t")


# select sig SNPs
subhapmap<-extract_hapmap(hapmap, gwas)
# get pheno and geno
geno_pheno<-get_pheno_geno(subhapmap, pheno)
# plot the boxplot
boxplot_geno_pheno(geno_pheno, outfolder)




























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
haplotypes
## Step 6: Get Sample Haplotypes and Save

snpcombsample <- HaploTraitR::getHapCombSamples(haplotypes, hapmap)
snpcombsample <- as.data.frame(snpcombsample)
write.csv(snpcombsample, file = file.path(outfolder, "haplotype_combinations.csv"), row.names = FALSE)

## Step 7: Read Phenotype Data



## Step 8: Generate SNP Combination Tables

SNPcombTables <- getSNPcombTables(snpcombsample, pheno)

## Step 9: Perform t-tests on SNP Combinations


t_test_snpComp <-  HaploTraitR::testSNPcombs(SNPcombTables)

## Step 10: Visualize Results


for (cls_snp in names(t_test_snpComp)) {
       HaploTraitR::plotHapCombBoxPlot(cls_snp, SNPcombTables, t_test_snpComp, outfolder= outfolder)
       ggplot2::ggsave(file.path(outfolder, paste0(cls_snp, "_boxplot.png")), width = 8, height = 6)
}

### Plot Haplotype Distribution

lDplots<-plotLDForClusters(ld_matrix_folder=LDsInfo[[1]], clusterLDs, outfolder=outfolder)

### Plot LD Combination Matrix with LD

LDCombs<-plotLDCombMatrix(LDsInfo[[1]], clusterLDs, haplotypes, gwas, outfolder=outfolder)
