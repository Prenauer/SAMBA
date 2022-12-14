% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SAMBA_Functions.R
\name{Samba}
\alias{Samba}
\title{Samba}
\usage{
Samba(
  data,
  design,
  coefficient = NULL,
  contrast = NULL,
  ntc.as.null.dist = T,
  score.method = "MetaAnalysis",
  test.method = "QLF",
  GuideMap = NULL,
  file.prefix = NULL
)
}
\arguments{
\item{data}{data.frame with read counts. Must have "sgRNA" and "Gene" columns, followed by columns with raw sample counts.}

\item{design}{design matrix for the samples in "data". Note: The order of the rows must match that of the data read-counts.}

\item{coefficient}{character name or vector, indicating which coefficient of the linear model is tested to be equal to zero.}

\item{contrast}{contrast matrix, indicating the contrast(s) of the linear model to be tested as equal to zero (optional).}

\item{ntc.as.null.dist}{logical, indicating whether to ese NTCs or all guides to generate a null distribution for gene-level analyses.}

\item{score.method}{character, indicating whether to use the 'MetaAnalysis' or 'GeneScore' method to perform gene-level analysis.}

\item{test.method}{character, indicating whether to use 'QLF' or 'LRT' tests.}

\item{GuideMap}{data.frame that maps the guides of the library to their respective genes (optional). Must have "sgRNA" and "Gene" columns.}

\item{file.prefix}{Prefix for writing output files.}
}
\value{
DGEList object of the edgeR package.
}
\description{
This is an all-in-one function to preprocesses and analyze read count data from a CRISPR
screen. The structure of the input data is assumed to have columns labeled "sgRNA" and
"Gene", followed by read counts from each sample in appropriately labeled columns.
}
