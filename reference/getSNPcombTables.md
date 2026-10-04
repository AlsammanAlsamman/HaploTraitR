# Build the accession-level haplotype and phenotype table

Build the accession-level haplotype and phenotype table

## Usage

``` r
getSNPcombTables(haplotypes, pheno, savecopy = TRUE)
```

## Arguments

- haplotypes:

  A haplotype table with samples (from \[getHapCombSamples()\])

- pheno:

  A phenotype data frame (accession, trait); see \[prepare_pheno()\]

- savecopy:

  Save a copy (\`tables/accession_haplotype_phenotype.csv\`) in the
  output folder

## Value

A long data frame with columns \`Sample\`, \`Pheno\`, \`block\`, \`hap\`
and \`lead_snps\` (phenotyped accessions only).
