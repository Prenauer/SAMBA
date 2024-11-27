#' CRISPR Screen analysis with moderated Bayesian statistics and adaptive gene aggregation scoring
#'
#' This analysis algorithm employs generalized linear models (GLMs) to fit CRISPR
#' screen data, and it analyzes guide-level enrichment/depletion using moderated Bayesian
#' statistics and quasi-likelihood F tests. Gene-level scores are then calculated and compared
#' with a null distribution to calculate enrichment and depletion for the screen.
#'
#' @docType package
#'
#' @author Paul Renauer \email{paul.renauer@yale.edu}
#'
#' @name SAMBA
#'
#' @import dplyr
#' @import stringr
#' @import edgeR
#' @import limma
#' @import dqrng
#' @import utils
#' @import stats
NULL
