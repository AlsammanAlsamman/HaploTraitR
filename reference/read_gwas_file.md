# Read GWAS results file with specified columns

The column names are taken from the configuration (\`rsid_col\`,
\`chr_col\`, \`pos_col\`, \`pval_col\`); see \[set_config()\].

## Usage

``` r
read_gwas_file(gwasfile, header = TRUE, sep = NULL)
```

## Arguments

- gwasfile:

  Path to the GWAS file (tab, comma or semicolon separated; may be
  gzipped)

- header:

  Logical, whether the file contains a header row

- sep:

  Separator used in the file. \`NULL\` (default) detects it
  automatically.

## Value

A data frame with the standardised columns \`rsid\`, \`chr\`, \`pos\`,
\`p\` and \`rs\` (\`chr:pos\` identifier).

## Examples

``` r
gwas <- read_gwas_file(system.file("extdata", "gwas_area_ann19.csv.gz", package = "HaploTraitR"))
head(gwas)
#>        Trait               rsid chr    pos       p MarkerR2        LOD Allel
#> 1 Area_ANN19     SCRI_RS_169369  1H 142369 0.57641 1.15e-03 0.23926849     C
#> 2 Area_ANN19 JHI-Hv50k-2016-358  1H 145558 0.88809 7.27e-05 0.05154302     C
#> 3 Area_ANN19 JHI-Hv50k-2016-343  1H 146097 0.77031 3.13e-04 0.11333446     A
#> 4 Area_ANN19 JHI-Hv50k-2016-342  1H 146155 0.73977 4.05e-04 0.13090329     C
#> 5 Area_ANN19     SCRI_RS_149838  1H 146218 0.54869 1.32e-03 0.26067295     A
#> 6 Area_ANN19     SCRI_RS_227491  1H 147494 0.56155 1.24e-03 0.25061157     A
#>     Effect        rs
#> 1 -0.21060 1H:142369
#> 2  0.07235 1H:145558
#> 3  0.13457 1H:146097
#> 4  0.18553 1H:146155
#> 5 -0.22529 1H:146218
#> 6  0.21533 1H:147494
```
