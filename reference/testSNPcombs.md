# Test trait differences between haplotypes

For every haplotype block, all pairs of haplotypes (with at least two
phenotyped accessions) are compared with a Welch t-test or a Wilcoxon
test (\`t_test_method\`), and p-values are adjusted with
\`p_adjust_method\`.

## Usage

``` r
testSNPcombs(comb_sample_Tables, savecopy = TRUE)
```

## Arguments

- comb_sample_Tables:

  The accession table (from \[getSNPcombTables()\])

- savecopy:

  Save the results (\`tables/\<trait\>\_pairwise_tests.csv\`) in the
  output folder

## Value

A named list (per block) of data frames with the columns \`group1\`,
\`group2\`, \`n1\`, \`n2\`, \`mean1\`, \`mean2\`, \`diff\`,
\`statistic\`, \`df\`, \`p\`, \`p.adj\`, \`p.adj.signif\`, \`block\`,
\`y.position\`, \`xmin\` and \`xmax\`. \`NULL\` if no block could be
tested.
