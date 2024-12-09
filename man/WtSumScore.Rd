% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Samba_functions.R
\name{WtSumScore}
\alias{WtSumScore}
\title{WtSumScore}
\usage{
WtSumScore(data, fdr.thresh)
}
\arguments{
\item{scores}{numeric vector of the log-FC values for all guides of a single gene.}

\item{sg.threshold}{numeric, giving the logFC value of the null data that represents the FDR threshold}
}
\value{
numeric score of a gene.
}
\description{
This function calculates a single gene score for a given set of guide log-FC values.
}
