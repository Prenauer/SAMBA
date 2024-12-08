% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Samba_functions.R
\name{CreateDesignMatrix}
\alias{CreateDesignMatrix}
\title{CreateDesignMatrix}
\usage{
CreateDesignMatrix(hdr.names, screen.names, ctrl.names)
}
\arguments{
\item{hdr.names}{character vector of the count data column names.}

\item{screen.names}{character vector of the screen sample names.}

\item{ctrl.names}{character vector of the control sample names.}
}
\value{
data.frame of the design matrix.
}
\description{
This function generates a simple design matrix when the user inputs a count file and a
list of names for the screen and control samples.
}
