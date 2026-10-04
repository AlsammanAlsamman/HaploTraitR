# Plot the trait distribution of each haplotype

Violin + box plot of the trait for every haplotype of a block, with the
number of accessions, the mean (diamond), the population mean (dashed
line) and letter groups. The superior haplotype is shown in green.

## Usage

``` r
plotHapCombBoxPlot(
  cls_snp = NULL,
  SNPcombTables,
  t_test_snpComp,
  with_dots = TRUE,
  outfolder = NULL,
  show_pairs = NULL
)
```

## Arguments

- cls_snp:

  A block id or GWAS SNP. If \`NULL\`, the first block is plotted.

- SNPcombTables:

  The accession table (from \[getSNPcombTables()\])

- t_test_snpComp:

  Pairwise tests (from \[testSNPcombs()\])

- with_dots:

  Show the individual accessions

- outfolder:

  Folder in which to save the plot (\`NULL\`: do not save)

- show_pairs:

  Draw brackets for significant pairs (\`NULL\`: only when there are at
  most four haplotypes)

## Value

A ggplot object.
