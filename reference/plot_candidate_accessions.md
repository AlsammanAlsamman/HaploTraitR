# Plot the best candidate accessions

Heat map of the haplotypes carried by the top-ranked accessions (see
\[rank_accessions()\]) in the significant blocks; green cells mark the
superior haplotypes. Useful to choose parents that stack favourable
haplotypes.

## Usage

``` r
plot_candidate_accessions(ranking, summary, top = 30)
```

## Arguments

- ranking:

  The data frame returned by \[rank_accessions()\]

- summary:

  The list returned by \[summarize_haplotypes()\]

- top:

  Number of accessions to show

## Value

A ggplot object (also saved as \`plots/candidate_accessions\`).
