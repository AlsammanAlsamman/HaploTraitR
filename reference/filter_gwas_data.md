# Filter GWAS data by FDR threshold

Filter GWAS data by FDR threshold

## Usage

``` r
filter_gwas_data(gwas)
```

## Arguments

- gwas:

  Data frame with the GWAS data (from \[read_gwas_file()\])

## Value

A list of data frames, one per chromosome, with the significant SNPs.

## Examples

``` r
gwas <- read_gwas_file(system.file("extdata", "gwas_area_ann19.csv.gz", package = "HaploTraitR"))
set_config(list(fdr_threshold = 0.1))
sig <- filter_gwas_data(gwas)
#> ℹ FDR column not found; FDR computed with method "fdr".
#> ℹ 5 significant GWAS SNPs (FDR < 0.1).
reset_config()
#> ✔ HaploTraitR configuration reset to default.
```
