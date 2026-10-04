# Find diagnostic SNPs for the superior haplotypes

For each block with a significant effect, lists the SNPs whose allele in
the superior haplotype differs from the other haplotypes. Fully
diagnostic SNPs distinguish the superior haplotype from all others and
are good candidates for marker-assisted selection (e.g. KASP assays).

## Usage

``` r
diagnostic_markers(haplotypes, summary, significant_only = TRUE)
```

## Arguments

- haplotypes:

  A haplotype table (from \[getHapCombSamples()\])

- summary:

  The list returned by \[summarize_haplotypes()\]

- significant_only:

  Only report blocks with a significant effect

## Value

A data frame with one row per informative SNP.
