# Estimation of selection and GC-biased gene conversion intensity from site frequency spectrum data by a least-square approach
# Sylvain Glemin
# sylvain.glemin@univ-rennes.fr


# Functions giving the expected number of SNP (ens) in frequency j/n under different models


#' @title Expectation of the neutral model
#'
#' @param theta population scale mutation rate = 4Neu
#' @param j number of copies of the derived allele in the sample
#'
#' @returns The expected number of SNPs in copy j in the sample
#'
#' @importFrom rlang abort
#'
#' @export
#'
#' @examples
#' theta <- 100
#' n <- 10
#' sfs <- ens_neutral(theta,j = 1:(n-1))
#' barplot(sfs)
ens_neutral <- function(theta,j) {
  #Error checking
  if(!is.numeric(theta) || theta <=0) {
    abort("theta must be a positive numeric value")
  }
  if(!is.numeric(j) || j <=0) {
    abort("j must be a positive integer")
  }
  #Main
  return(theta/j)
}
ens_neutral <- Vectorize(ens_neutral,vectorize.args = "j")



#' @title Expectation of constant gBGC with intensity B
#'
#' @description
#'  B can be positive (W->S mutations) or negative (S->W) mutations
#'
#' @param theta population scale mutation rate = 4Neu
#' @param B population scale biased gene conversion coefficient = 4Neb
#' @param n sample size
#' @param j number of copies of the derived allele
#'
#' @returns The expected number of SNPs in copy j in the sample
#'
#' @export
#'
#' @examples
#' theta <- 100
#' B <- 1.2
#' n <- 10
#' sfs <- ens_constant(theta,B,n,j = 1:(n-1))
#' barplot(sfs)
ens_constant <- function(theta,B,n,j) {
  #Error checking
  if(!is.numeric(theta) || theta <=0) {
    abort("theta must be a positive numeric value")
  }
  if(!is.numeric(B)) {
    abort("B must be a positive numeric value")
  }
  if(!is.numeric(n) || n<=1) {
    abort("n must be an integer higher than 1")
  }
  if(!is.numeric(j) || j<=0 || j>=n) {
    abort("j must be a positive integer lower than n")
  }
  #Main
  if (B==0) return(ens_neutral(theta,j))
  ens <- ((n*theta)/(j*(n-j))) * ((1-exp(-B*(1-j/n)))/(1-exp(-B)))
  return(ens)
}
ens_constant <- Vectorize(ens_constant,vectorize.args = "j")


#' @title Expectation of the first hotspot model
#'
#' @description
#' gBGC with intensity B in hotspot (fraction f) and 0 otherwise (fraction 1-f)
#'
#' @param theta population scale mutation rate = 4Neu
#' @param B population scale biased gene conversion coefficient in hotspots = 4Neb
#' @param f fraction of hotspots
#' @param n sample size
#' @param j number of copies of the derived allele
#'
#' @returns The expected number of SNPs in copy j in the sample
#'
#' @export
#'
#' @examples
#' theta <- 100
#' B <- 1.2
#' f <- 0.2
#' n <- 10
#' sfs <- ens_hotspot1(theta,B,f,n,j = 1:(n-1))
#' barplot(sfs)
ens_hotspot1 <- function(theta,B,f,n,j) {
  #Error checking
  if(!is.numeric(theta) || theta <=0) {
    abort("theta must be a positive numeric value")
  }
  if(!is.numeric(B)) {
    abort("B must be a positive numeric value")
  }
  if(!is.numeric(f) || f<0 || f>1) {
    abort("f must be a numeric value between 0 and 1")
  }
  if(!is.numeric(n) || n<=1) {
    abort("n must be an integer higher than 1")
  }
  if(!is.numeric(j) || j<=0 || j>=n) {
    abort("j must be a positive integer lower than n")
  }
  #Main
  if (B==0) return(ens_neutral(theta,j))
  ens <- f * ens_constant(theta,B,n,j) + (1-f) * ens_neutral(theta,j)
  return(ens)
}
ens_hotspot1 <- Vectorize(ens_hotspot1,vectorize.args = "j")


#' @title Expectation of the second hotspot model
#' @description
#' gBGC with intensity B1 in hotspot (fraction f) and B0 otherwise (fraction 1-f)
#'
#' @param theta population scale mutation rate = 4Neu
#' @param B0 population scale biased gene conversion coefficient in the background = 4Neb0
#' @param B1 population scale biased gene conversion coefficient in hotspots = 4Neb1
#' @param f fraction of hotspots
#' @param n sample size
#' @param j number of copies of the derived allele
#'
#' @returns The expected number of SNPs in copy j in the sample
#'
#' @export
#'
#' @examples
#' theta <- 100
#' B0 <- 0.2
#' B1 <- 3.5
#' f <- 0.2
#' n <- 10
#' sfs <- ens_hotspot2(theta,B0,B1,f,n,j = 1:(n-1))
#' barplot(sfs)
ens_hotspot2 <- function(theta,B0,B1,f,n,j) {
  #Error checking
  if(!is.numeric(theta) || theta <=0) {
    abort("theta must be a positive numeric value")
  }
  if(!is.numeric(B0)) {
    abort("B0 must be a positive numeric value")
  }
  if(!is.numeric(B1)) {
    abort("B1 must be a positive numeric value")
  }
  if(!is.numeric(f) || f<0 || f>1) {
    abort("f must be a numeric value between 0 and 1")
  }
  if(!is.numeric(n) || n<=1) {
    abort("n must be an integer higher than 1")
  }
  if(!is.numeric(j) || j<=0 || j>=n) {
    abort("j must be a positive integer lower than n")
  }
  #Main
  if (B0==0) return(ens_hotspot1(theta,B1,f,n,j))
  if (B1==0) return(ens_hotspot1(theta,B0,1-f,n,j))
  if (B0==0 & B1==0) return(ens_neutral(theta,j))
  ens <-f * ens_constant(theta,B1,n,j) + (1-f) * ens_constant(theta,B0,n,j)
  return(ens)
}
ens_hotspot2 <- Vectorize(ens_hotspot2,vectorize.args = "j")



# EXPECTATIONS WITH ORIENTATION ERRORS ####

#' @title Expectation of the neutral model with polarization error
#'
#' @param theta1 population scale mutation rate of the focal category = 4Neu1
#' @param theta2 population scale mutation rate of the complementary category = 4Neu2
#' @param e1 polarization error rate of the focal category
#' @param e2 polarization error rate of the complementary category
#' @param n sample size
#' @param j number of copies of the derived allele
#'
#' @returns The expected number of SNPs in copy j in the sample
#'
#' @export
#'
#' @examples
#' theta1 <- 100
#' theta2 <- 200
#' e1 <- 0.01
#' e2 <- 0.02
#' n <- 10
#' sfs <- ens_neutral_err(theta1,theta2,e1,e2,n,j = 1:(n-1))
#' barplot(sfs)
ens_neutral_err <- function(theta1,theta2,e1,e2,n,j) {
  #Error checking
  if(!is.numeric(theta1) || theta1 <=0) {
    abort("theta1 must be a positive numeric value")
  }
  if(!is.numeric(theta2) || theta2 <=0) {
    abort("theta2 must be a positive numeric value")
  }
  if(!is.numeric(e1) || e1<0 || e1>1) {
    abort("e1 must be a numeric value between 0 and 1")
  }
  if(!is.numeric(e2) || e2<0 || e2>1) {
    abort("e2 must be a numeric value between 0 and 1")
  }
  if(!is.numeric(n) || n<=1) {
    abort("n must be an integer higher than 1")
  }
  if(!is.numeric(j) || j<=0 || j>=n) {
    abort("j must be a positive integer lower than n")
  }
  #Main
  ens <- (1-e1)*ens_neutral(theta1,j) + e2*ens_neutral(theta2,n-j)
  return(ens)
}
ens_neutral_err <- Vectorize(ens_neutral_err,vectorize.args = "j")



#' @title Expectation of constant gBGC with intensity B
#'
#' @description
#' B can be positive (W->S mutations) or negative (S->W) mutations
#'
#' @param theta1 population scale mutation rate of the focal category = 4Neu1
#' @param theta2 population scale mutation rate of the complementary category = 4Neu2
#' @param B population scale biased gene conversion coefficient in hotspots = 4Neb
#' @param e1 polarization error rate of the focal category
#' @param e2 polarization error rate of the complementary category
#' @param n sample size
#' @param j number of copies of the derived allele
#'
#' @returns The expected number of SNPs in copy j in the sample
#'
#' @export
#'
#' @examples
#' theta1 <- 100
#' theta2 <- 200
#' B <- 1.2
#' e1 <- 0.01
#' e2 <- 0.02
#' n <- 10
#' sfs <- ens_constant_err(theta1,theta2,B,e1,e2,n,j = 1:(n-1))
#' barplot(sfs)
ens_constant_err <- function(theta1,theta2,B,e1,e2,n,j) {
  #Error checking
  if(!is.numeric(theta1) || theta1 <=0) {
    abort("theta1 must be a positive numeric value")
  }
  if(!is.numeric(theta2) || theta2 <=0) {
    abort("theta2 must be a positive numeric value")
  }
  if(!is.numeric(B)) {
    abort("B0 must be a positive numeric value")
  }
  if(!is.numeric(e1) || e1<0 || e1>1) {
    abort("e1 must be a numeric value between 0 and 1")
  }
  if(!is.numeric(e2) || e2<0 || e2>1) {
    abort("e2 must be a numeric value between 0 and 1")
  }
  if(!is.numeric(n) || n<=1) {
    abort("n must be an integer higher than 1")
  }
  if(!is.numeric(j) || j<=0 || j>=n) {
    abort("j must be a positive integer lower than n")
  }
  #Main
  ens <- (1-e1)*ens_constant(theta1,B,n,j) + e2*ens_constant(theta2,-B,n,n-j)
  return(ens)
}
ens_constant_err <- Vectorize(ens_constant_err,vectorize.args = "j")



#' @title Expectation of the first hotspot model
#'
#' @description
#' gBGC with intensity B in hotspot (fraction f) and 0 otherwise (fraction 1-f)
#'
#' @param theta1 population scale mutation rate of the focal category = 4Neu1
#' @param theta2 population scale mutation rate of the complementary category = 4Neu2
#' @param B population scale biased gene conversion coefficient in hotspots = 4Neb
#' @param f fraction of hotpsots
#' @param e1 polarization error rate of the focal category
#' @param e2 polarization error rate of the complementary category
#' @param n sample size
#' @param j number of copies of the derived allele
#'
#' @returns The expected number of SNPs in copy j in the sample
#'
#' @export
#'
#' @examples
#' theta1 <- 100
#' theta2 <- 200
#' B <- 1.2
#' f <- 0.2
#' e1 <- 0.01
#' e2 <- 0.02
#' n <- 10
#' sfs <- ens_hotspot1_err(theta1,theta2,B,f,e1,e2,n,j = 1:(n-1))
#' barplot(sfs)
ens_hotspot1_err <- function(theta1,theta2,B,f,e1,e2,n,j) {
  #Error checking
  if(!is.numeric(theta1) || theta1 <=0) {
    abort("theta1 must be a positive numeric value")
  }
  if(!is.numeric(theta2) || theta2 <=0) {
    abort("theta2 must be a positive numeric value")
  }
  if(!is.numeric(B)) {
    abort("B must be a positive numeric value")
  }
  if(!is.numeric(f) || f<0 || f>1) {
    abort("f must be a numeric value between 0 and 1")
  }
  if(!is.numeric(e1) || e1<0 || e1>1) {
    abort("e1 must be a numeric value between 0 and 1")
  }
  if(!is.numeric(e2) || e2<0 || e2>1) {
    abort("e2 must be a numeric value between 0 and 1")
  }
  if(!is.numeric(n) || n<=1) {
    abort("n must be an integer higher than 1")
  }
  if(!is.numeric(j) || j<=0 || j>=n) {
    abort("j must be a positive integer lower than n")
  }
  #Main
  ens <- (1-e1)*ens_hotspot1(theta1,B,f,n,j) + e2*ens_hotspot1(theta2,-B,f,n,n-j)
  return(ens)
}
ens_hotspot1_err <- Vectorize(ens_hotspot1_err,vectorize.args = "j")



#' @title Expectation of the second hotspot model
#'
#' @description
#' gBGC with intensity B1 in hotspot (fraction f) and B0 otherwise (fraction 1-f)
#'
#' @param theta1 population scale mutation rate of the focal category = 4Neu1
#' @param theta2 population scale mutation rate of the complementary category = 4Neu2
#' @param B0 population scale biased gene conversion coefficient in the background = 4Neb0
#' @param B1 population scale biased gene conversion coefficient in hotspots = 4Neb1
#' @param f fraction of hotpsots
#' @param e1 polarization error rate of the focal category
#' @param e2 polarization error rate of the complementary category
#' @param n sample size
#' @param j number of copies of the derived allele
#'
#' @returns The expected number of SNPs in copy j in the sample
#'
#' @export
#'
#' @examples
#' theta1 <- 100
#' theta2 <- 200
#' B0 <- 0.2
#' B1 <- 3.5
#' f <- 0.2
#' e1 <- 0.01
#' e2 <- 0.02
#' n <- 10
#' sfs <- ens_hotspot2_err(theta1,theta2,B0,B1,f,e1,e2,n,j = 1:(n-1))
#' barplot(sfs)
ens_hotspot2_err <- function(theta1,theta2,B0,B1,f,e1,e2,n,j) {
  #Error checking
  if(!is.numeric(theta1) || theta1 <=0) {
    abort("theta1 must be a positive numeric value")
  }
  if(!is.numeric(theta2) || theta2 <=0) {
    abort("theta2 must be a positive numeric value")
  }
  if(!is.numeric(B0)) {
    abort("B0 must be a positive numeric value")
  }
  if(!is.numeric(B1)) {
    abort("B1 must be a positive numeric value")
  }
  if(!is.numeric(f) || f<0 || f>1) {
    abort("f must be a numeric value between 0 and 1")
  }
  if(!is.numeric(e1) || e1<0 || e1>1) {
    abort("e1 must be a numeric value between 0 and 1")
  }
  if(!is.numeric(e2) || e2<0 || e2>1) {
    abort("e2 must be a numeric value between 0 and 1")
  }
  if(!is.numeric(n) || n<=1) {
    abort("n must be an integer higher than 1")
  }
  if(!is.numeric(j) || j<=0 || j>=n) {
    abort("j must be a positive integer lower than n")
  }
  #Main
  ens <- (1-e1)*ens_hotspot2(theta1,B0,B1,f,n,j) + e2*ens_hotspot2(theta2,-B0,-B1,f,n,n-j)
  return(ens)
}
ens_hotspot2_err <- Vectorize(ens_hotspot2_err,vectorize.args = "j")
