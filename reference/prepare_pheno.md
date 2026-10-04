# Prepare phenotype data

Selects the trait column, removes missing values and averages replicated
measurements of the same accession.

## Usage

``` r
prepare_pheno(pheno)
```

## Arguments

- pheno:

  A data frame whose first column holds accession names and other
  columns hold traits. If it has more than two columns the trait is
  taken from the \`pheno_col\` configuration parameter.

## Value

A data frame with columns \`Sample\` and \`Pheno\`.

## Examples

``` r
pheno <- data.frame(Taxa = c("G1", "G1", "G2"), Yield = c(1, 3, 5))
prepare_pheno(pheno)
#> ℹ Replicated measurements found; the mean per accession is used.
#>   Sample Pheno
#> 1     G1     2
#> 2     G2     5
```
