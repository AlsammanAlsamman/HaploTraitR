# Rank accessions by the superior haplotypes they carry

Builds an accession x block matrix of haplotypes and counts, for every
accession, how many superior haplotypes of significant blocks it
carries. Accessions that stack many superior haplotypes are candidate
parents. Genotyped accessions without phenotype data are included.

## Usage

``` r
rank_accessions(haplotypes, pheno, summary)
```

## Arguments

- haplotypes:

  A haplotype table with samples (from \[getHapCombSamples()\])

- pheno:

  A phenotype data frame (accession, trait)

- summary:

  The list returned by \[summarize_haplotypes()\]

## Value

A data frame with one row per accession.
