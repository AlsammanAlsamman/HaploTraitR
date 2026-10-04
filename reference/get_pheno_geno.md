# Get the phenotype and genotype data

Get the phenotype and genotype data

## Usage

``` r
get_pheno_geno(subhapmap, pheno)
```

## Arguments

- subhapmap:

  A list of HapMap data frames (from \[extract_hapmap()\])

- pheno:

  A phenotype data frame (accession, trait)

## Value

A long data frame with columns \`Phenotype\`, \`variable\` (SNP),
\`value\` (genotype call) and \`Sample\`.
