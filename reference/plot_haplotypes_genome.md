# Plot the haplotype blocks across the genome

Manhattan plot of the GWAS results in which the haplotype blocks are
highlighted (green: significant haplotype effect, orange: not
significant) and labelled.

## Usage

``` r
plot_haplotypes_genome(
  gwas_data,
  haplotype_data,
  use_facet = TRUE,
  pwidth = 12,
  pheight = 5,
  pdpi = 300,
  summary = NULL
)
```

## Arguments

- gwas_data:

  The full GWAS results (data frame from \[read_gwas_file()\], or a list
  of data frames)

- haplotype_data:

  A haplotype table (from \[getHapCombSamples()\])

- use_facet:

  If \`TRUE\` a genome-wide plot is drawn; if \`FALSE\` one zoomed plot
  per chromosome carrying a block

- pwidth, pheight, pdpi:

  Size (inches) and resolution of the saved plot(s)

- summary:

  Optional output of \[summarize_haplotypes()\] used to colour the
  blocks

## Value

A named list of ggplot objects (also saved in \`plots/\`).
