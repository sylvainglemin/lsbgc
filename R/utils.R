# Estimation of selection and GC-biased gene conversion intensity from site frequency spectrum data by a least-square approach
# Sylvain Glemin
# sylvain.glemin@univ-rennes.fr

# Various functions to manipulate SFS or compute statistics


# Harmonic number

#' @title Harmonic number
#'
#' @param n an integer
#'
#' @returns the harmonic number of n
#'
#' @export
#'
#' @examples
#' hn(10)
hn <- function(n) {
  sum(1/c(1:n))
}



#' @title Skewness of the SFS
#' @description
#' Function to compute the skweness of a frequency spectrum and testing the asymmetry
#' Typically for testing the asymmetry of the GC unfolded spectrum
#'
#' @param sfs a SFS
#'
#' @returns the skewness of the sfs
#'
#' @importFrom stats pnorm
#'
#' @export
#'
#' @examples
#' sfs <- c(100,50,30,20,30,40,80)
#' skewness_sfs(sfs)
skewness_sfs <- function(sfs) {
  Nclass <- length(sfs)+1
  Nobs <- sum(sfs)
  freq <- c(1:(Nclass-1))/Nclass
  mean <- sum(sfs*freq)/Nobs
  skew <- (sqrt(Nobs*(Nobs-1))/(Nobs-2))*sum(sfs*(freq-mean)^3/Nobs)/(sum(sfs*(freq-mean)^2)/Nobs)^(3/2)
  SES <- sqrt(6*Nobs*(Nobs-1)/((Nobs-2)*(Nobs+1)*(Nobs+3))) # Standard error of skewness
  pval <- 2*pnorm(abs(skew/SES),0,1,lower.tail=F)
  return(list(skewness=skew,p.value=pval))
}



#' @title Projection of the SFS
#'
#' @param sfs initial SFS of size n = 1 + length(sfs)
#' @param m size of the projected sfs with m ≤ n
#'
#' @returns the projected sfs
#'
#' @export
#'
#' @examples
#' sfs <- c(100,50,30,20,30,40,80)
#' project_sfs(sfs,5)
project_sfs <- function(sfs,m){
  n <- length(sfs)
  output <- rep(0,m)
  index<-c(1:m)
  for(k in 1:n){
    output <- output + sfs[k]*choose(k,index)*choose(n-k,m-index)/choose(n,m)
  }
  return(output)
}
