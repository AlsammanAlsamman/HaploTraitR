% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_haplotypes_genome.R
\name{plot_haplotypes_genome}
\alias{plot_haplotypes_genome}
\title{Plot Haplotypes SNPs Across the Genome}
\usage{
plot_haplotypes_genome(
  gwas_data,
  haplotype_data,
  use_facet = TRUE,
  pwidth = 10,
  pheight = 5,
  pdpi = 300
)
}
\arguments{
\item{gwas_data}{A data frame or list of data frames containing the GWAS results.}

\item{haplotype_data}{A data frame containing the haplotype information.}

\item{use_facet}{Logical, whether to use facet wrapping for chromosomes.}
}
\value{
A list of ggplot2 objects representing the plots.
}
\description{
Plot Haplotypes SNPs Across the Genome
}
