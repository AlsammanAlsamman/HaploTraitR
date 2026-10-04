# Plot LD and haplotype matrices of all blocks

Saves one \[plot_haplotype_block()\] figure per block in
\`plots/haplotype_blocks\`.

## Usage

``` r
plotLDCombMatrix(clusterLDs = NULL, haplotypes, gwas = NULL, summary = NULL)
```

## Arguments

- clusterLDs:

  Not used; kept for compatibility with HaploTraitR \< 0.2

- haplotypes:

  A haplotype table (from \[getHapCombSamples()\])

- gwas:

  Not used; kept for compatibility with HaploTraitR \< 0.2

- summary:

  Optional output of \[summarize_haplotypes()\]

## Value

Invisibly, a named list of \`haplotraitr_figure\` objects.
