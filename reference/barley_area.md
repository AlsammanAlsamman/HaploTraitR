# Barley leaf area example data

GWAS results, genotypes and phenotypes for leaf area in a barley panel
(275 accessions, Barley 50K SNP array). The same data are available as
files in \`system.file("extdata", package = "HaploTraitR")\`.

## Usage

``` r
barley_area
```

## Format

A list with three elements:

- gwas:

  GWAS results (data frame) for leaf area

- hapmap:

  Genotypes as returned by \[readHapmap()\] (list of data frames per
  chromosome)

- pheno:

  Phenotypes: accession name (\`Taxa\`) and leaf area (\`Area\`, cm2)

## Source

ICARDA barley breeding program.
