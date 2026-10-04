# Run the complete HaploTraitR analysis

One call that reads the data, builds haplotype blocks around the
significant GWAS SNPs, tests the haplotype effects on the trait, and
writes all plots, CSV tables and a formatted Excel report to the output
folder.

## Usage

``` r
run_haplotraitr(
  gwas,
  genotypes,
  pheno,
  trait = NULL,
  unit = NULL,
  outfolder = NULL,
  higher_is_better = TRUE,
  config = list(),
  snp_plots = TRUE
)
```

## Arguments

- gwas:

  GWAS results: a file path or a data frame from \[read_gwas_file()\]

- genotypes:

  Genotypes: a HapMap or VCF file path, or the list returned by
  \[readHapmap()\] / \[readVCF()\]

- pheno:

  Phenotypes: a file path or a data frame whose first column holds the
  accession names

- trait:

  Name of the trait column in \`pheno\` (default: the second column)

- unit:

  Trait unit shown on the plots (optional)

- outfolder:

  Output folder. Default: a new \`haplotraitR_run\` folder in the
  working directory

- higher_is_better:

  \`TRUE\` if higher trait values are desirable (yield), \`FALSE\`
  otherwise (disease score, days to heading...)

- config:

  Other configuration parameters as a named list, see \[list_config()\]

- snp_plots:

  Also draw one plot per significant GWAS SNP

## Value

Invisibly, a \`haplotraitr_result\` list with the elements \`summary\`,
\`haplotypes\`, \`accessions\`, \`pairwise\`, \`markers\`, \`ranking\`,
\`gwas_significant\`, \`LDblocks\` and \`outfolder\`.

## Details

Output folder layout: \* \`HaploTraitR\_\<trait\>.xlsx\` - the Excel
report (start here) \* \`plots/genome_overview\` - GWAS Manhattan plot
with the haplotype blocks \* \`plots/haplotype_effects_summary\` -
effect of every haplotype in every block \*
\`plots/candidate_accessions\` - accessions carrying the most superior
haplotypes \* \`plots/haplotype_blocks/\` - alleles, trait effects and
LD of each block \* \`plots/haplotype_trait_plots/\` - trait
distribution per haplotype \* \`plots/gwas_snp_plots/\` - trait
distribution per genotype of each GWAS SNP \* \`tables/\` - CSV copies
of all tables; \`LD_matrices/\` - LD (r2) matrices

## Examples

``` r
# \donttest{
ex <- function(f) system.file("extdata", f, package = "HaploTraitR")
res <- run_haplotraitr(
  gwas = ex("gwas_area_ann19.csv.gz"),
  genotypes = ex("Barley_50K.tsv.gz"),
  pheno = ex("area_ann19.tsv"),
  trait = "Area", unit = "cm2",
  outfolder = file.path(tempdir(), "haplo_example"),
  config = list(fdr_threshold = 0.1)
)
#> 
#> ── HaploTraitR 0.2.0 ───────────────────────────────────────────────────────────
#> Trait: Area (cm2) • higher is better • output: /tmp/RtmpuvFqob/haplo_example
#> 
#> [1/8] Reading data
#>       ℹ 38 SNPs share a position with another SNP; only the first is kept.
#>       • GWAS: 36,534 SNPs
#>       • Genotypes: 6,704 SNPs × 275 accessions
#>       • Phenotypes: 275 accessions (275 also genotyped)
#> [2/8] Selecting significant GWAS SNPs
#>       ℹ FDR column not found; FDR computed with method "fdr".
#>       ℹ 5 significant GWAS SNPs (FDR < 0.1).
#> [3/8] Building LD blocks around the GWAS SNPs
#> [4/8] Identifying haplotypes
#>       ℹ 3 haplotype blocks from 5 GWAS SNPs.
#>       ℹ 275 of 275 genotyped accessions have phenotype data.
#> [5/8] Testing haplotype effects
#> [6/8] Drawing haplotype plots
#> [7/8] Drawing GWAS SNP plots
#> [8/8] Writing Excel report
#> 
#> ✔ Analysis finished in 12 secs
#> 
#> ── Results: Area (cm2) ──
#> 
#> ✔ 3 haplotype blocks, 1 with a significant haplotype effect
#> 
#>   Block               Haps  p value  Evidence  Superior  vs pop.  
#>   ──────────────────────────────────────────────────────────────
#>   2H:579.59-581.37Mb  2     4.0e-05  None      H2        -2.1%    
#>   2H:580.27-581.68Mb  2     1.4e-05  None      H2        -1.8%    
#>   2H:581.68-583.66Mb  3     1.6e-06  Strong    H3        +9.4%    
#> 
#> → 2H:581.68-583.66Mb: Select H3: mean 25.19 cm2 (+9.4% vs population); significantly higher than 2 of 2 other haplotype(s); carried by 30 accessions (11%).
#> ! In "2H:579.59-581.37Mb" and "2H:580.27-581.68Mb" the best group is Other (rare haplotypes); lower comb_freq_threshold to split it.
#> ℹ 17 diagnostic SNPs for the superior haplotypes (10 fully diagnostic, ready for KASP design)
#> ℹ Top candidate accessions: "G182", "G183", "G111", "G255", and "G180"
#> 
#> ── Output ──
#> 
#> • Excel report: /tmp/RtmpuvFqob/haplo_example/HaploTraitR_Area.xlsx
#> • Plots (14): /tmp/RtmpuvFqob/haplo_example/plots
#> • CSV tables: /tmp/RtmpuvFqob/haplo_example/tables
res
#> 
#> ── Results: Area (cm2) ──
#> 
#> ✔ 3 haplotype blocks, 1 with a significant haplotype effect
#> 
#>   Block               Haps  p value  Evidence  Superior  vs pop.  
#>   ──────────────────────────────────────────────────────────────
#>   2H:579.59-581.37Mb  2     4.0e-05  None      H2        -2.1%    
#>   2H:580.27-581.68Mb  2     1.4e-05  None      H2        -1.8%    
#>   2H:581.68-583.66Mb  3     1.6e-06  Strong    H3        +9.4%    
#> 
#> → 2H:581.68-583.66Mb: Select H3: mean 25.19 cm2 (+9.4% vs population); significantly higher than 2 of 2 other haplotype(s); carried by 30 accessions (11%).
#> ! In "2H:579.59-581.37Mb" and "2H:580.27-581.68Mb" the best group is Other (rare haplotypes); lower comb_freq_threshold to split it.
#> ℹ 17 diagnostic SNPs for the superior haplotypes (10 fully diagnostic, ready for KASP design)
#> ℹ Top candidate accessions: "G182", "G183", "G111", "G255", and "G180"
#> 
#> ── Output ──
#> 
#> • Excel report: /tmp/RtmpuvFqob/haplo_example/HaploTraitR_Area.xlsx
#> • Plots (14): /tmp/RtmpuvFqob/haplo_example/plots
#> • CSV tables: /tmp/RtmpuvFqob/haplo_example/tables
# }
```
