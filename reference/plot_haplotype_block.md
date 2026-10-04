# Plot a haplotype block: alleles, trait effects and LD

Produces a three-part figure for one haplotype block: \* top: the
alleles of every haplotype at every SNP (major allele in blue, minor
allele in orange; red triangles mark the GWAS SNPs); \* right: the trait
mean and 95 letter groups (green = superior haplotype, dashed line =
population mean); \* bottom: the physical position of the SNPs and their
pairwise LD (r2).

## Usage

``` r
plot_haplotype_block(block, haplotypes, summary = NULL)
```

## Arguments

- block:

  A block id or GWAS SNP of the block

- haplotypes:

  A haplotype table (from \[getHapCombSamples()\] or
  \[convertLDclusters2Haps()\])

- summary:

  Optional output of \[summarize_haplotypes()\]; adds the trait panel

## Value

A \`haplotraitr_figure\` object (print it to draw it). Its \`width\` and
\`height\` elements give a suitable size in inches.
