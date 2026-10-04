# Cluster SNPs based on LD matrices

SNPs are linked when their r2 is at least \`ld_threshold\`; the
connected group of SNPs that contains the GWAS SNP forms its LD block.
Blocks smaller than \`dist_cluster_count\` SNPs are discarded.

## Usage

``` r
getLDclusters(LDsInfo)
```

## Arguments

- LDsInfo:

  A list returned by \[computeLDclusters()\] or
  \[retrieveLDMatricesFromFolder()\]

## Value

A named list (per GWAS SNP) of the SNP identifiers in its LD block.
