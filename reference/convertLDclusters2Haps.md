# Convert LD clusters to haplotypes

GWAS SNPs whose LD blocks contain exactly the same SNPs are merged into
one haplotype block. Within each block, every distinct combination of
alleles carried by at least \`comb_freq_threshold\` of the accessions
becomes a haplotype (H1 = most frequent).

## Usage

``` r
convertLDclusters2Haps(hapmap, clusterLDs, savecopy = TRUE)
```

## Arguments

- hapmap:

  A list of HapMap data frames (from \[readHapmap()\])

- clusterLDs:

  A list of LD blocks (from \[getLDclusters()\])

- savecopy:

  Save a copy (\`tables/haplotypes_alleles.csv\`) in the output folder

## Value

A data frame with one row per haplotype: \`block\`, \`chr\`, \`start\`,
\`end\`, \`n_snps\`, \`lead_snps\`, \`snps\`, \`hap\`, \`alleles\` and
\`freq\`. SNP details and LD of each block are kept in the \`"blocks"\`
attribute.
