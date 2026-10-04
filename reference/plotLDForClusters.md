# Plot LD heat maps for all LD blocks

Plot LD heat maps for all LD blocks

## Usage

``` r
plotLDForClusters(clusterLDs, LDsInfo = NULL)
```

## Arguments

- clusterLDs:

  A list of LD blocks (from \[getLDclusters()\])

- LDsInfo:

  Optional output of \[computeLDclusters()\]; otherwise the matrices are
  read from the \`LD_matrices\` folder of the output folder

## Value

Invisibly, a list of ggplot objects (also saved in
\`plots/LD_heatmaps\`).
