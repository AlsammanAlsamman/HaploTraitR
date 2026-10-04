# Summarise haplotype effects for breeders

Computes, for every block and haplotype, the trait mean with its 95
confidence interval, the difference from the population mean, letter
groups (haplotypes sharing a letter do not differ significantly) and
identifies the superior haplotype according to \`higher_is_better\`.

## Usage

``` r
summarize_haplotypes(haplotypes, SNPcombTables, t_test_snpComp = NULL)
```

## Arguments

- haplotypes:

  A haplotype table with samples (from \[getHapCombSamples()\])

- SNPcombTables:

  The accession table (from \[getSNPcombTables()\])

- t_test_snpComp:

  Pairwise tests (from \[testSNPcombs()\]); computed if \`NULL\`

## Value

A list with two data frames: \`effects\` (one row per haplotype) and
\`blocks\` (one row per block, including a selection recommendation).
