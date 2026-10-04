# Generate boxplots of the trait for each significant GWAS SNP

Generate boxplots of the trait for each significant GWAS SNP

## Usage

``` r
boxplot_genotype_phenotype(
  pheno,
  gwas,
  hapmap,
  plot_width = 5,
  plot_height = 5.5,
  plot_all = TRUE
)
```

## Arguments

- pheno:

  A data frame containing phenotype data (accession, trait)

- gwas:

  A list of significant GWAS SNPs (from \[filter_gwas_data()\])

- hapmap:

  A list of HapMap data frames (from \[readHapmap()\])

- plot_width:

  The width of the saved plots (inches)

- plot_height:

  The height of the saved plots (inches)

- plot_all:

  If \`TRUE\` plot all SNPs, otherwise only SNPs with a significant
  genotype effect

## Value

A named list of ggplot objects.
