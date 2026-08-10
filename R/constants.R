# Estimation of selection and GC-biased gene conversion intensity from site frequency spectrum data by a least-square approach
# Sylvain Glemin
# sylvain.glemin@univ-rennes.fr

#' @title Internal Constants
#'
#' @description Constants used internally by the package.
#' @noRd
# Constants used when functions are not defined in 0
ZERO <- 10^(-8)
ONE <- 1 - ZERO
# Value below which the snp category is removed
SNPTHRESH <- 0.1
# Default parameters for optimization
# gBGC coefficient
BMIN <- -100
BMAX <- 100
# Log of the mutation bias
MMIN <- -2
MMAX <- 2
# Error rates
EMAX <- 0.49
EINIT <- 0.01
# Choice of the correction to applied to SFS
COR <- TRUE
# Option parameters of the optim function
MAXIT <- 100
FACTR <- 10^7
LMM <- 20
VERBOSE <- 0
USEGR <- FALSE
