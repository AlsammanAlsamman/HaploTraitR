# Assign accessions to haplotypes

Accessions whose alleles match a haplotype exactly are assigned to it.
Accessions with some missing calls (up to \`hap_max_missing\`) are
assigned to the single haplotype compatible with their observed alleles.
All other accessions are pooled as \`"Other"\`.

## Usage

``` r
getHapCombSamples(haplotypes, hapmap, savecopy = TRUE)
```

## Arguments

- haplotypes:

  A haplotype table (from \[convertLDclusters2Haps()\])

- hapmap:

  A list of HapMap data frames (from \[readHapmap()\])

- savecopy:

  Save a copy (\`tables/haplotype_accessions.csv\`) in the output folder

## Value

The haplotype table with an extra \`"Other"\` row per block and the
columns \`n\` (number of accessions), \`freq\` and \`samples\`
(\`\|\`-separated).
