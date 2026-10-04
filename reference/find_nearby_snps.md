# Find nearby SNPs to a given SNP list assuming they are in the same chromosome

Find nearby SNPs to a given SNP list assuming they are in the same
chromosome

## Usage

``` r
find_nearby_snps(snp_list1, snp_list2, threshold)
```

## Arguments

- snp_list1:

  A vector of SNP positions (reference)

- snp_list2:

  A vector of SNP positions (query)

- threshold:

  The maximum distance between the SNPs

## Value

A list with, for each position in \`snp_list1\`, the positions of
\`snp_list2\` within \`threshold\`.

## Examples

``` r
find_nearby_snps(c(1000, 2000, 3000), c(1500, 2500, 3500, 5000), 1000)
#> [[1]]
#> [1] 1500
#> 
#> [[2]]
#> [1] 1500 2500
#> 
#> [[3]]
#> [1] 2500 3500
#> 
```
