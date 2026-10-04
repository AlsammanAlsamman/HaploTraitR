# Extract HapMap data for significant SNPs

Extract HapMap data for significant SNPs

## Usage

``` r
extract_hapmap(hapmap, gwas)
```

## Arguments

- hapmap:

  A list of HapMap data frames (from \[readHapmap()\])

- gwas:

  A list of GWAS data frames (from \[filter_gwas_data()\])

## Value

A list of HapMap data frames containing only the GWAS SNPs.
