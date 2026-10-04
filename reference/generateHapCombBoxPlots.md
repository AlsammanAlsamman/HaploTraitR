# Generate and save the haplotype trait plots of all blocks

Plots are saved in \`plots/haplotype_trait_plots/significant\` or
\`.../not_significant\` according to the overall test of the block.

## Usage

``` r
generateHapCombBoxPlots(
  SNPcombTables,
  t_test_snpComp,
  pwidth = NULL,
  pheight = NULL
)
```

## Arguments

- SNPcombTables:

  The accession table (from \[getSNPcombTables()\])

- t_test_snpComp:

  Pairwise tests (from \[testSNPcombs()\])

- pwidth, pheight:

  Not used anymore; the size adapts to the number of haplotypes.

## Value

Invisibly, a named list of ggplot objects.
