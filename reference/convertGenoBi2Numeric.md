# Convert HapMap genotypes to allele dosages

Convert HapMap genotypes to allele dosages

## Usage

``` r
convertGenoBi2Numeric(habmap)
```

## Arguments

- habmap:

  A HapMap data frame (11 annotation columns followed by genotype
  columns), or a character matrix of genotype calls.

## Value

An integer matrix (SNPs x accessions) counting copies of the minor
allele (0, 1, 2); missing calls are \`NA\`.

## Examples

``` r
g <- data.frame(matrix(NA, 2, 11), G1 = c("AA", "CC"), G2 = c("AG", "CT"), G3 = c("GG", "NN"))
convertGenoBi2Numeric(g)
#>      G1 G2 G3
#> [1,]  0  1  2
#> [2,]  0  1 NA
```
