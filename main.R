
library(HaploTraitR)


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
head(gwas)

hapmap <- readHapmap(hapfile)


## New addition
pheno <- read.csv(phenofile, header = TRUE, sep = "\t")
head(pheno)


# select sig SNPs
subhapmap <- list()
for (i in 1:length(gwas)) {
  chr <- names(gwas)[i]
  sigSNPs <- gwas[[i]]$rs
  subhapmap[[i]] <- hapmap[[chr]][sigSNPs,]
}
# rbinding the hapmap
subhapmap <- do.call(rbind, subhapmap)
# remove the first 11 columns
subhapmap[,1:11]<-NULL
# transpose the data
subhapmap<-t(subhapmap)
snps_target<-colnames(subhapmap)
# merge the pheno data with the hapmap
geno_pheno<-merge(subhapmap, pheno, by.x="row.names", by.y="Taxa")
# transfere row.names to the first column
rownames(geno_pheno)<-geno_pheno[,1]
geno_pheno<-geno_pheno[,-1]
#melting the data and using the SNP and Taxa as the id variables
library(reshape2)
Trait_name<-"PH"


geno_pheno_melt<-melt(geno_pheno, id.vars=c(Trait_name))




#################################
library(ggplot2)
library(ggpubr)  # For statistical annotations
library(rstatix)

# Assuming `geno_pheno_melt` is the data frame and `snps_target` contains SNPs to plot
head(geno_pheno_melt)

# Extract unique SNP values
snp_var <- unique(geno_pheno_melt$value)

# Generate consistent colors for SNP values
snp_colors <- rainbow(length(snp_var))
names(snp_colors) <- snp_var  # Assign names for consistent mapping

# Loop through the SNP targets
for (snp in snps_target) {
  # Subset data for the SNP
  sub_data <- geno_pheno_melt[geno_pheno_melt$variable %in% snp, ]
  colnames(sub_data) <- c("Trait", "snp_var", "value")
  sub_data$value <- as.factor(sub_data$value)
  # remove the NN values
  sub_data<-sub_data[!sub_data$value=="NN",]
  stat.test <- NULL
  # Perform the t-test using compare_means()
  tryCatch({
    stat.test <- compare_means(
      formula = Trait ~ value,  # Define the formula for comparison
      data = sub_data,          # Use the subset data
      method = "t.test"         # Specify Welch Two Sample t-test
    )
    # Add y.position for annotations
    stat.test <- stat.test %>%
      mutate(y.position = max(sub_data$Trait) + 1)  # Position above the boxplot
  }, error = function(e) {
    print(paste("Error in", snp, ":", e))
  })

  # Create the boxplot with annotations
  p <- ggplot(sub_data, aes(x = value, y = Trait)) +
    geom_violin(aes(fill = value, alpha = 0.01)) +
    geom_jitter(width = 0.2, alpha = 0.3) +
    scale_fill_manual(values = snp_colors)+
    # no legend
    theme(legend.position = "none")+
    # labs title
    labs(title = paste("Phenotype Distribution by SNP", snp),
         x = "Haplotype Combination",
         y = "Phenotype Value")
  if (!is.null(stat.test)) {
    p<- p + stat_pvalue_manual(stat.test[1, ], label = "p.signif", tip.length = 0.01)
  }
  ggsave(paste0(snp, "_boxplot.png"), plot = p, path = outfolder, width = 6, height = 6)
}










##############################################
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
