# Read a VCF file into HapMap format

Converts the \`GT\` field of a (possibly gzipped) VCF file into the
HapMap structure used by HaploTraitR, so VCF data can be used directly
in the pipeline. Only bi- and multi-allelic SNPs are kept (indels are
skipped); heterozygous and missing calls are preserved (\`AG\`, \`NN\`).

## Usage

``` r
readVCF(vcfpath)
```

## Arguments

- vcfpath:

  Path to the VCF file

## Value

A list of data frames, one per chromosome, as returned by
\[readHapmap()\].
