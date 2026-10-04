# Read a HapMap genotype file

Reads a HapMap file (11 annotation columns followed by one column per
accession). Genotypes can be two-letter (\`AA\`, \`AG\`) or
single-letter IUPAC codes (\`A\`, \`R\`); missing calls (\`NN\`, \`N\`,
\`-\`) are standardised to \`NN\`.

## Usage

``` r
readHapmap(path, sep = "\t", header = TRUE, ...)
```

## Arguments

- path:

  Path to the HapMap file (may be gzipped)

- sep:

  Separator (default tab)

- header:

  Logical indicating if the file has a header

- ...:

  Additional arguments passed to \[utils::read.table()\]

## Value

A list of data frames, one per chromosome, with row names \`chr:pos\`.

## Examples

``` r
hapmap <- readHapmap(system.file("extdata", "Barley_50K.tsv.gz", package = "HaploTraitR"))
#> ℹ 38 SNPs share a position with another SNP; only the first is kept.
names(hapmap)
#> [1] "2H" "6H" "7H"
```
