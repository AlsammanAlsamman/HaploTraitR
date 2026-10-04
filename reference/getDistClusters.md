# Get distance clusters around the significant GWAS SNPs

For every significant GWAS SNP, collects the genotyped SNPs within
\`dist_threshold\` bp. Windows with fewer than \`dist_cluster_count\`
SNPs are dropped.

## Usage

``` r
getDistClusters(gwas, hapmap)
```

## Arguments

- gwas:

  A list of GWAS data frames split by chromosome (from
  \[filter_gwas_data()\])

- hapmap:

  A list of HapMap data frames split by chromosome (from
  \[readHapmap()\])

## Value

A list (per chromosome) of named lists (per GWAS SNP) of SNP positions.
