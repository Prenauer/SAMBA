% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Samba_functions.R
\name{Analyze_Samba_Genes}
\alias{Analyze_Samba_Genes}
\title{Analyze_Samba_Genes}
\usage{
Analyze_Samba_Genes(
  sgRes,
  control.gene = NULL,
  score.method = "GeneScore",
  fdr.threshold = 0.1,
  file.prefix = NULL,
  verbose = T
)
}
\arguments{
\item{sgRes}{data.frame of the results from the guide-level screen enrichment analysis. This is the output from the "Analyze_Samba_Guides" function.}

\item{control.gene}{character, gene label of the internal control guides (e.g. "NTC").}

\item{score.method}{character, indicating whether to use the 'KS' (Kolmogorov-Smirnov) or 'GeneScore' method to perform gene-level analysis.}

\item{fdr.threshold}{numeric, indicating the proportion of the control guides, or total guides if control.gene = NULL, will be used as a threshold for the weighted sum score.}

\item{file.prefix}{(optional) file.path and/or character prefix to use for exporting the result tables}

\item{verbose}{logical, indicating whether to display verbose output information.}
}
\value{
data.frame of the results from the gene-level screen enrichment analysis.
}
\description{
This function performs the gene-level screen enrichment analysis of SAMBA, using guide-level SAMBA results as an input.
}
