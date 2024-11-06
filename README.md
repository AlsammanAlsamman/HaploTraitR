# Haplotype Analysis Project

This project performs haplotype analysis and analyzing and visualizing phenotypic data across multiple genetic / haplotypes combinations 
to understand the assocaition between haplotype (SNPs linked by LD) and variation in traits. 

## Sample Files

- **Haplotype File:** `Barley_50K_KNNimp.hmp.txt`
- **GWAS File:** `SignificantSNP_GWAS.csv`
- **Phenotype File:** `Pheno_ANN19.tsv`


## test Parameters

- **Distance Threshold:** 1,000,000 (for clustering all SNPs within 1Mb of the significant SNPs)
- **Cluster Count Threshold:** 5 (minimum number of SNPs in a cluster)
- **LD Threshold:** 0.3 (minimum LD value to consider two SNPs in LD)
- **Combination Frequency Threshold:** 0.1 (minimum genotypes with the same combination to consider the combination)

## Example code Steps

1. **Read GWAS Data:**
    ```r
    gwas <- readGWAS(gwasfile, sep = ",")
    ```

2. **Read Haplotype Data:**
    ```r
    hapmap <- readHapmap(hapfile)
    ```

3. **Cluster SNPs by Chromosome:**
    ```r
    haplotype_clusters <- getHaplotypeClusters(gwas, hapmap, dist_threshold, dist_cluster_count)
    ```

4. **Compute LD Matrices:**
    ```r
    LDsInfo <- computeLDclusters(hapmap, haplotype_clusters)
    ```

5. **Save LD Matrices to Output Folder:**
    ```r
    saveLDs2folder(LDsInfo, outfolder)
    ```

6. **Cluster LDs:**
    ```r
    clusterLDs <- clusterLD(LDsInfo, ld_threshold, cls_count = 3)
    ```

7. **Convert LD Clusters to Haplotypes:**
    ```r
    haplotypes <- convertLDclusters2Haps(hapmap, clusterLDs, comb_freq_threshold)
    ```

8. **Get Haplotype Combination Samples:**
    ```r
    snpcombsample <- getHapCombSamples(haplotypes, hapmap)
    snpcombsample <- as.data.frame(snpcombsample)
    ```

9. **Save Haplotype Combinations:**
    ```r
    write.csv(snpcombsample, file = file.path(outfolder, "haplotype_combinations.csv"), row.names = FALSE)
    ```

10. **Read Phenotype Data:**
    ```r
    pheno <- read.csv(phenofile, header = TRUE, sep = "\t")
    ```

11. **Get SNP Combination Tables:**
    ```r
    SNPcombTables <- getSNPcombTables(snpcombsample, pheno)
    ```

12. **Perform T-Test on SNP Combinations:**
    ```r
    t_test_snpComp <- testSNPcombs(SNPcombTables)
    ```

13. **Plot Haplotype Combination Box Plot:**
    ```r
    snp <- "2H_JHI-Hv50k-2016-79696"
    plotHapCombBoxPlot(snp, SNPcombTables, t_test_snpComp)
    ```

14. **Plot Haplotype Combination Distribution:**
    ```r
    plotHapCombDistribution(snp, SNPcombTables)
    ```
