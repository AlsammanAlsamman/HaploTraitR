# List all the configuration parameters

List all the configuration parameters

## Usage

``` r
list_config()
```

## Value

A data frame with the parameter names, current values and descriptions
(printed and returned invisibly).

## Examples

``` r
list_config()
#> 
#> ── HaploTraitR settings ──
#> 
#>   outfolder            NULL       Output folder for all results
#>   dist_threshold       1,000,000  Window (bp) around each GWAS SNP in which SNPs are collected
#>   dist_cluster_count   5          Minimum number of SNPs in a window / LD block
#>   ld_threshold         0.3        Minimum LD (r2) to link two SNPs into the same block
#>   comb_freq_threshold  0.1        Minimum frequency of a haplotype (rarer ones are pooled as 'Other')
#>   hap_max_missing      0.2        Maximum fraction of missing calls allowed when assigning an accession to a haplotype
#>   fdr_threshold        0.05       FDR threshold used to select significant GWAS SNPs
#>   fdr_method           fdr        Method used by p.adjust() to compute GWAS FDR
#>   t_test_threshold     0.05       Significance level (alpha) for haplotype tests
#>   t_test_method        t.test     Test used to compare haplotypes: 't.test' (Welch) or 'wilcox.test'
#>   p_adjust_method      holm       Multiple-testing correction for pairwise haplotype tests
#>   higher_is_better     TRUE       TRUE if a higher trait value is desirable (e.g. yield), FALSE otherwise (e.g. disease score)
#>   rsid_col             rsid       Marker name column in the GWAS file
#>   pos_col              pos        Position column in the GWAS file
#>   chr_col              chr        Chromosome column in the GWAS file
#>   pval_col             p          P-value column in the GWAS file
#>   fdr_col              fdr        FDR column in the GWAS file (computed if absent)
#>   phenotypename        Phenotype  Trait name used in titles and file names
#>   pheno_col            Phenotype  Trait column in the phenotype file
#>   phenotypeunit        NULL       Trait unit shown on plot axes
#>   plot_format          png        Plot file format: 'png' or 'pdf'
#>   plot_dpi             300        Resolution of PNG plots
#> Changed values are yellow. Use set_config(list(name = value)) to change a
#> setting.
```
