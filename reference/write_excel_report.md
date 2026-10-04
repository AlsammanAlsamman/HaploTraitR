# Write the HaploTraitR results to a formatted Excel workbook

The workbook contains the sheets \*Read me\*, \*Block summary\*,
\*Haplotype effects\*, \*Diagnostic markers\*, \*Candidate accessions\*,
\*Pairwise tests\*, \*GWAS SNPs\* and \*Parameters\*. Superior
haplotypes are highlighted in green.

## Usage

``` r
write_excel_report(results, file = NULL)
```

## Arguments

- results:

  The list returned by \[run_haplotraitr()\] (or a list with the
  elements \`summary\`, \`haplotypes\`, \`pairwise\`, \`markers\`,
  \`ranking\`, \`gwas_significant\`)

- file:

  Path of the \`.xlsx\` file. Default:
  \`\<outfolder\>/HaploTraitR\_\<trait\>.xlsx\`

## Value

Invisibly, the file path.
