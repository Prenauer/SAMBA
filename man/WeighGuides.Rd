% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SAMBA_Functions.R
\name{WeighGuides}
\alias{WeighGuides}
\title{WeighGuides}
\usage{
WeighGuides(dge)
}
\arguments{
\item{dge}{DGEList object with raw count data.}
}
\value{
DGEList object that contains guide weights as a numeric vector.
}
\description{
This function calculates a weight for each guide, based on the number of screen samples that have detected counts.
}
