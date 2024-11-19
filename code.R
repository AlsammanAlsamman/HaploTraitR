rm(list = ls())
library(HaploTraitR)

# Define file paths and output folder
hapfile <- "sampledata/Barley_50K_KNNimp.hmp.txt"
gwasfile <- "sampledata/SignificantSNP_GWAS.tsv"
phenofile <- "sampledata/Pheno_ANN19.tsv"
outfolder <- "sampleout"

#Set thresholds for distance, clustering, LD, and combination frequency.

dist_threshold <- 1000000 # Distance threshold (1Mb)
dist_cluster_count <- 5   # Minimum SNPs in a cluster
ld_threshold <- 0.3       # Minimum LD value
comb_freq_threshold <- 0.1 # Minimum genotype frequency for a combination

#First, load the GWAS data to identify significant SNPs for clustering.
gwas <- readGWAS(gwasfile, sep = "\t")


#Load the hapmap file to access genotype data for the analysis.

hapmap <- readHapmap(hapfile)

#Cluster SNPs within the specified distance threshold around significant SNPs from the GWAS data.

dist_clusters <- getDistClusters(gwas, hapmap, dist_threshold, dist_cluster_count)


############ Its beeter to select the cluster that have the significant SNP and remove others


## Step 3: Compute LD Matrices

#Calculate linkage disequilibrium (LD) values for each SNP cluster.

LDsInfo <- computeLDclusters(hapmap, dist_clusters, outfolder)

## Step 4: Save LD Matrices

#Save the LD matrices to the specified output folder if no folder is provided, a temporary folder will be created.

saveLDs2folder(LDsInfo, outfolder)

## Step 5: Cluster LD Values and Convert to Haplotypes

#Cluster SNPs based on LD values and apply frequency thresholding to derive haplotypes.


clusterLDs <- getLDclusters(LDsInfo,ld_threshold, cls_count = 3)

# calculate haplotypes Linked SNPs from the clusters
haplotypes <- convertLDclusters2Haps(hapmap, clusterLDs, comb_freq_threshold)

#write.csv(haplotypes, file = file.path(outfolder, "haplotypes.csv"), row.names = FALSE)
## Step 6: Get Sample Haplotypes and Save

#Extract the haplotype combinations for samples and save to a CSV file.
snpcombsample <- getHapCombSamples(haplotypes, hapmap)


############# Working updating
#write.csv(snpcombsample, file = file.path(outfolder, "haplotype_combinations.csv"), row.names = FALSE)


## Step 7: Read Phenotype Data

#Read the phenotype data for associating haplotypes with phenotypic traits.

pheno <- read.csv(phenofile, header = TRUE, sep = "\t")

## Step 8: Generate SNP Combination Tables

#Generate tables to link SNP combinations with phenotype data.

SNPcombTables <- getSNPcombTables(snpcombsample, pheno)

## Step 9: Perform t-tests on SNP Combinations

#Perform statistical testing on SNP combinations to see if there are significant associations with phenotypes.


#library(magrittr)
#library(rstatix)

t_test_snpComp <-  testSNPcombs(SNPcombTables)

## Step 10: Visualize Results

### Plot Haplotypes as Boxplot

#Select an SNP and create a boxplot to compare phenotypic values among different haplotype combinations.


#library(ggplot2)
#library(ggpubr)


for (cls_snp in names(t_test_snpComp)) {
       plotHapCombBoxPlot(cls_snp, SNPcombTables, t_test_snpComp)
       ggplot2::ggsave(file.path(outfolder, paste0(cls_snp, "_boxplot.png")), width = 8, height = 6)
}

plotLDCombMatrix(LDsInfo[[1]], clusterLDs, haplotypes, gwas, outfolder=outfolder)


lDplots<-plotLDForClusters(ld_matrix_folder=LDsInfo[[1]], clusterLDs, outfolder=outfolder)


