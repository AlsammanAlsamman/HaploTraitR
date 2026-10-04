# Compute LD matrices for the distance clusters

Computes the r2 (squared correlation of allele dosages) between all SNPs
of each window and saves the matrices in the \`LD_matrices\` folder.

## Usage

``` r
computeLDclusters(hapmap, haplotype_clusters)
```

## Arguments

- hapmap:

  A list of HapMap data frames (from \[readHapmap()\])

- haplotype_clusters:

  A list of distance clusters (from \[getDistClusters()\])

## Value

A list with the LD folder (\`ld_folder\`), cluster information
(\`out_info\`) and the LD matrices (\`matrices\`).
