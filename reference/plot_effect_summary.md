# Plot the haplotype effects of all blocks

One row per block; each point is a haplotype placed at its difference
from the population mean (in proportional to its frequency. The superior
haplotype is green.

## Usage

``` r
plot_effect_summary(summary)
```

## Arguments

- summary:

  The list returned by \[summarize_haplotypes()\]

## Value

A ggplot object (also saved as \`plots/haplotype_effects_summary\`).
