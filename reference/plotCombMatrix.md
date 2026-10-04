# Plot the haplotype allele matrix of a block

Plot the haplotype allele matrix of a block

## Usage

``` r
plotCombMatrix(
  cluster_id,
  haplotypes,
  gwas = NULL,
  snps = NULL,
  outfolder = NULL
)
```

## Arguments

- cluster_id:

  A block id or GWAS SNP

- haplotypes:

  A haplotype table (from \[convertLDclusters2Haps()\] or
  \[getHapCombSamples()\])

- gwas:

  Not used; kept for compatibility

- snps:

  Not used; kept for compatibility

- outfolder:

  Folder in which to save the plot (\`NULL\`: do not save)

## Value

A ggplot object.
