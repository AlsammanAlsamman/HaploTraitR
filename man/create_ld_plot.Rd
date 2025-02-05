% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/snplinkage_subCodes.R
\name{create_ld_plot}
\alias{create_ld_plot}
\title{Function to create LD plot and integrate with SNP positions}
\usage{
create_ld_plot(
  df_ld,
  df_snp,
  biplot_subset = NULL,
  labels_colname = NULL,
  diamonds = length(unique(df_ld$SNP_A)) < 40,
  point_size = ifelse(is.null(biplot_subset), 120, 80)/sqrt(nrow(df_ld)),
  title = "",
  title_biplot = "",
  ...
)
}
\arguments{
\item{df_ld}{Data frame with columns SNP_A, SNP_B, and R2}

\item{df_snp}{Data frame with columns snpID, position, and labels_colname}

\item{biplot_subset}{Subset of SNPs to highlight in the biplot}

\item{labels_colname}{Column name for SNP labels}

\item{diamonds}{Should the values be displayed as diamonds or points ?}

\item{point_size}{Size for geom_point}

\item{title}{Title for the plot}

\item{title_biplot}{Title for the biplot}
}
\value{
ggplot
}
\description{
Function to create LD plot and integrate with SNP positions
}
