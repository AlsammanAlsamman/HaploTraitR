rm(list = ls())
library(HaploTraitR)

# Define file paths and output folder
hapfile <- "/home/samman/Documents/MyGitHub/HaploTraitR/sampledata/Barley_50K_KNNimp.hmp.txt"
gwasfile <- "/home/samman/Documents/MyGitHub/HaploTraitR/sampledata/SignificantSNP_GWAS.tsv"
phenofile <- "/home/samman/Documents/MyGitHub/HaploTraitR/sampledata/Pheno_ANN19.tsv"
outfolder <- "/home/samman/Documents/MyGitHub/HaploTraitR/outfolder"

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
dist_clusters

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
snpcombsample
ncol(snpcombsample)
colnames(snpcombsample)
haplotypes
############# Working updating
#write.csv(snpcombsample, file = file.path(outfolder, "haplotype_combinations.csv"), row.names = FALSE)


## Step 7: Read Phenotype Data

#Read the phenotype data for associating haplotypes with phenotypic traits.

pheno <- read.csv(phenofile, header = TRUE, sep = "\t")

## Step 8: Generate SNP Combination Tables

#Generate tables to link SNP combinations with phenotype data.

SNPcombTables <- getSNPcombTables(snpcombsample, pheno)


# only for the documentation
# plot the haplotypes as boxplot
library(ggplot2)
library(ggpubr)
# Sample    Pheno                                                                         ldcls comb         SNP
# 1     G2 70.56720 CC|TT|AA|TT|GG|AA|GG|AA|AA|GG|TT|AA|CC|AA|GG|GG|AA|AA|TT|GG|CC|CC|AA|GG|AA|CC    1 2H:48081081
# 2    G17 78.18191 CC|TT|AA|TT|GG|AA|GG|AA|AA|GG|TT|AA|CC|AA|GG|GG|AA|AA|TT|GG|CC|CC|AA|GG|AA|CC    1 2H:48081081
# 3    G25 80.10286 CC|TT|AA|TT|GG|AA|GG|AA|AA|GG|TT|AA|CC|AA|GG|GG|AA|AA|TT|GG|CC|CC|AA|GG|AA|CC    1 2H:48081081
# 4    G26 66.19045 CC|TT|AA|TT|GG|AA|GG|AA|AA|GG|TT|AA|CC|AA|GG|GG|AA|AA|TT|GG|CC|CC|AA|GG|AA|CC    1 2H:48081081
# 5    G33 84.09275 CC|TT|AA|TT|GG|AA|GG|AA|AA|GG|TT|AA|CC|AA|GG|GG|AA|AA|TT|GG|CC|CC|AA|GG|AA|CC    1 2H:48081081
# 6    G36 78.57061 CC|TT|AA|TT|GG|AA|GG|AA|AA|GG|TT|AA|CC|AA|GG|GG|AA|AA|TT|GG|CC|CC|AA|GG|AA|CC    1 2H:48081081




## Step 9: Perform t-tests on SNP Combinations

#Perform statistical testing on SNP combinations to see if there are significant associations with phenotypes.


#library(magrittr)
#library(rstatix)

t_test_snpComp <-  testSNPcombs(SNPcombTables)
t_test_snpComp


library(ggplot2)
library(ggpubr)
ggplot(SNPcombTables, aes(x = comb, y = Pheno)) +
  geom_violin(aes(fill = comb), alpha = 0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = paste("Phenotype Distribution by Haplotype Combination for SNP", "2H:48081081"),
       x = "Haplotype Combination",
       y = "Phenotype Value") +
  theme(legend.position = "none")+
  facet_wrap(~SNP)

# compare shared genotypes between the haplotypes



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
