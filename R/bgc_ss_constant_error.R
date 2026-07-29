# Estimation of selection and GC-biased gene conversion intensity from site frequency spectrum data by a least-square approach
# Sylvain Glemin
# sylvain.glemin@univ-rennes.fr


# Series of SS functions and their corresponding gradient where parameters are constant
# Include polarization errors




######################################### #
# NO BGC, MUTATION BIAS ###################
######################################### #

#' @title Sum of squares of model M with errors
#'
#' @description Function that return the weighted sum of squares between the function of the SFSs and error rates and the linear predictor.
#' The theory predict that the expectation of log(Tws(j)/Tsw(j)) = -log(mutbias) + log(GC) - log(1 - GC)
#' where Tws and Tsw are the "True" SFSs and j/n the frequency class j.
#' Tws and Tsw can be expressed as a function of the Observed SFSs (Ows and Osw) and error rates (e1 and e2):
#' Tws = Ows(1-e2) - rev(Osw)e2 and Tsw = Osw(1-e1) - rev(Ows)*e1
#' The function return the weighted least-square as afunction of: sum(w(j) * (log(Tws(j)/Tsw(j)) - B * j/n + log(mut_bias) )^2)
#' where Tws and Tsw are expressed as a function of Ows, Osw, e1 and e2
#' w(j) = Ows(j)*Osw(j) / (Ows(j)+Osw(j)) is the weight used in the least square
#'
#' @param par a vector with the for parameters of the model.
#' par(1) = M (i.e. log(mut_bias))
#' par(2) = e1
#' par(3) = e2
#' @param WS the WS observed SFS
#' @param SW the SW observed SFS
#' @param GC GC content
#' @param cor a Boolean to add a correction to the transformed SFS (default = TRUE)
#'
#' @returns The weighted sum of squares
#'
#' @export
#'
#' @examples
#' sfsWS <- c(100,50,30,15,10)
#' sfsSW <- c(200,80,30,10,5)
#' param <- c(2,0.02,0.01)
#' sum_of_squares_M_err(param,sfsWS,sfsSW,0.5)
sum_of_squares_M_err <- function(par,WS,SW,GC,cor=COR) {
  #Error checking
  if(!is.numeric(c(WS,SW)) || length(which(c(WS,SW)<0))>0 ) {
    abort("The two SFSs must have positive numeric values")
  }
  if(length(WS)!=length(SW)) {
    abort("The two SFSs, WS and SW, must have the same length")
  }
  if(length(par)!=3) {
    abort("A three values vector must be given as par")
  }
  if(GC<=0 | GC>=1) {
    abort("GC content must be strictly between 0 and 1")
  }
  #Main
  M <- par[1]
  e1 <- par[2]
  e2 <- par[3]
  n <- length(WS)
  # True SFS as a function of observed one.
  WSt <- ((1 - e2)*WS - e2*rev(SW))/(1 - e1 - e2)
  SWt <- ((1 - e1)*SW - e1*rev(WS))/(1 - e1 - e2)
  # Same expression without 1-e1-e2 that simplifies in ratios
  #WSt2 <- ((1 - e2)*WS - e2*rev(SW))
  #SWt2 <- ((1 - e1)*SW - e1*rev(WS))
  # y <- log(WSt2/SWt2) + 1/(WSt) - 1/(SWt)
  removeNA <- which(WSt>=1 & SWt>=1)
  WSt <- WSt[removeNA]
  SWt <- SWt[removeNA]
  w <- 1/(t_variance(WSt,cor) + t_variance(SWt,cor))
  x <- c(1:n)/(n+1)
  x <- x[removeNA]
  y <- t_sfs(WSt,cor) - t_sfs(SWt,cor)
  ypred <- rep(- M - log(GC) + log(1 - GC),n)
  ypred <- ypred[removeNA]
  return( sum(w*(y-ypred)^2)/sum(w) )
}


#' @title Gradient of Sum of squares model M with errors
#'
#' @description Gradient of the sum of squares function
#'
#' @param par a vector with the for parameters of the model.
#' par(1) = M (i.e. log(mut_bias))
#' par(2) = e1
#' par(3) = e2
#' @param WS the WS observed SFS
#' @param SW the SW observed SFS
#' @param GC the GC content
#' @param cor a Boolean to add a correction to the transformed SFS (default = TRUE)
#'
#' @returns The gradient function
#'
#' @noRd
gr_sum_of_squares_M_err <- function(par,WS,SW,GC,cor=COR) {
  #Error checking
  if(!is.numeric(c(WS,SW)) || length(which(c(WS,SW)<0))>0 ) {
    abort("The two SFSs must have positive numeric values")
  }
  if(length(WS)!=length(SW)) {
    abort("The two SFSs, WS and SW, must have the same length")
  }
  if(length(par)!=3) {
    abort("A three values vector must be given as par")
  }
  if(GC<=0 | GC>=1) {
    abort("GC content must be strictly between 0 and 1")
  }
  #Main
  M <- par[1]
  e1 <- par[2]
  e2 <- par[3]
  n <- length(WS)
  # True SFS as a function of observed one.
  WSt <- ((1 - e2)*WS - e2*rev(SW))/(1 - e1 - e2)
  SWt <- ((1 - e1)*SW - e1*rev(WS))/(1 - e1 - e2)
  removeNA <- which(WSt>=1 & SWt>=1)
  WSt <- WSt[removeNA]
  SWt <- SWt[removeNA]
  w <- 1/(t_variance(WSt,cor) + t_variance(SWt,cor))
  x <- c(1:n)/(n+1)
  x <- x[removeNA]
  y <- t_sfs(WSt,cor) - t_sfs(SWt,cor)
  ypred <- rep(- M - log(GC) + log(1 - GC),n)
  # Derivative of y as a function of error rates
  # Approximated version, derivative of log(WSt2/SWt2) without additional terms
  dy1 <- d_expected_log_ratio(WS,SW,e1,e2)$d1
  dy2 <- d_expected_log_ratio(WS,SW,e1,e2)$d2
  grM <- sum(w*(2*(y - ypred))/sum(w))
  gre1 <- sum(w*(2*(y - ypred)*dy1)/sum(w))
  gre2 <- sum(w*(2*(y - ypred)*dy2)/sum(w))
  return(c(grM,gre1,gre2))
}



######################################### #
# CONSTANT BGC, NO MUTATION BIAS ##########
######################################### #


#' @title Sum of squares of model B with errors
#'
#' @description Function that return the weighted sum of squares between the function of the SFSs and error rates and the linear predictor.
#' The theory predict that the expectation of log(Tws(j)/Tsw(j)) = B * j/n + log(GC) - log(1 - GC)
#' where Tws and Tsw are the "True" SFSs, B = 4Neb is the population-scaled gBGC, and j/n the frequency class j.
#' Tws and Tsw can be expressed as a function of the Observed SFSs (Ows and Osw) and error rates (e1 and e2):
#' Tws = Ows(1-e2) - rev(Osw)e2 and Tsw = Osw(1-e1) - rev(Ows)*e1
#' The function return the weighted least-square as afunction of: sum(w(j) * (log(Tws(j)/Tsw(j)) - B * j/n + log(mut_bias) )^2)
#' where Tws and Tsw are expressed as a function of Ows, Osw, e1 and e2
#' w(j) = Ows(j)*Osw(j) / (Ows(j)+Osw(j)) is the weight used in the least square
#'
#' @param par a vector with the three parameters of the model.
#' par(1) = B
#' par(2) = e1
#' par(3) = e2
#' @param WS the WS observed SFS
#' @param SW the SW observed SFS
#' @param GC the GC content
#' @param cor a Boolean to add a correction to the transformed SFS (default = TRUE)
#'
#' @returns The weighted sum of squares
#'
#' @export
#'
#' @examples
#' sfsWS <- c(100,50,30,15,10)
#' sfsSW <- c(200,80, 30, 10, 5)
#' param <- c(1,0.02,0.01)
#' sum_of_squares_B_err(param,sfsWS,sfsSW,0.5)
sum_of_squares_B_err <- function(par,WS,SW,GC,cor=COR) {
  #Error checking
  if(!is.numeric(c(WS,SW)) || length(which(c(WS,SW)<0))>0 ) {
    abort("The two SFSs must have positive numeric values")
  }
  if(length(WS)!=length(SW)) {
    abort("The two SFSs, WS and SW, must have the same length")
  }
  if(length(par)!=3) {
    abort("A three values vector must be given as par")
  }
  if(GC<=0 | GC>=1) {
    abort("GC content must be strictly between 0 and 1")
  }
  #Main
  B <- par[1]
  e1 <- par[2]
  e2 <- par[3]
  n <- length(WS)
  # True SFS as a function of observed one.
  WSt <- ((1 - e2)*WS - e2*rev(SW))/(1 - e1 - e2)
  SWt <- ((1 - e1)*SW - e1*rev(WS))/(1 - e1 - e2)
  removeNA <- which(WSt>=1 & SWt>=1)
  WSt <- WSt[removeNA]
  SWt <- SWt[removeNA]
  w <- 1/(t_variance(WSt,cor) + t_variance(SWt,cor))
  x <- c(1:n)/(n+1)
  x <- x[removeNA]
  y <- t_sfs(WSt,cor) - t_sfs(SWt,cor)
  ypred <- B*x - log(GC) + log(1 - GC)
  return( sum(w*(y-ypred)^2)/sum(w) )
}


#' @title Gradient of Sum of squares of model B with errors
#'
#' @description Gradient of the sum of squares function
#'
#' @param par a vector with the four parameters of the model.
#' par(1) = B
#' par(2) = e1
#' par(3) = e2
#' @param WS the WS observed SFS
#' @param SW the SW observed SFS
#' @param GC the GC content
#' @param cor a Boolean to add a correction to the transformed SFS (default = TRUE)
#'
#' @returns The gradient function
#'
#' @noRd
gr_sum_of_squares_B_err <- function(par,WS,SW,GC,cor=COR) {
  #Error checking
  if(!is.numeric(c(WS,SW)) || length(which(c(WS,SW)<0))>0 ) {
    abort("The two SFSs must have positive numeric values")
  }
  if(length(WS)!=length(SW)) {
    abort("The two SFSs, WS and SW, must have the same length")
  }
  if(length(par)!=3) {
    abort("A three values vector must be given as par")
  }
  if(GC<=0 | GC>=1) {
    abort("GC content must be strictly between 0 and 1")
  }
  #Main
  B <- par[1]
  e1 <- par[2]
  e2 <- par[3]
  n <- length(WS)
  # True SFS as a function of observed one.
  WSt <- ((1-e2)*WS - e2*rev(SW))/(1 - e1 - e2)
  SWt <- ((1-e1)*SW - e1*rev(WS))/(1 - e1 - e2)
  removeNA <- which(WSt>=1 & SWt>=1)
  WSt <- WSt[removeNA]
  SWt <- SWt[removeNA]
  w <- 1/(t_variance(WSt,cor) + t_variance(SWt,cor))
  x <- c(1:n)/(n+1)
  x <- x[removeNA]
  y <- t_sfs(WSt,cor) - t_sfs(SWt,cor)
  dy1 <- d_expected_log_ratio(WS,SW,e1,e2)$d1
  dy2 <- d_expected_log_ratio(WS,SW,e1,e2)$d2
  grB <- sum(w*(-2*x*(y - ypred))/sum(w))
  gre1 <- sum(w*(2*(y - ypred)*dy1)/sum(w))
  gre2 <- sum(w*(2*(y - ypred)*dy2)/sum(w))
  return(c(grB,gre1,gre2))
}


######################################### #
# CONSTANT BGC, MUTATION BIAS #############
######################################### #

#' @title Sum of squares of model BM with errors
#'
#' @description Function that return the weighted sum of squares between the function of the SFSs and error rates and the linear predictor.
#' The theory predict that the expectation of log(Tws(j)/Tsw(j)) = B * j/n - log(mut_bias) + log(GC) - log(1 - GC)
#' where Tws and Tsw are the "True" SFSs, B = 4Neb is the population-scaled gBGC, mut_bias the mutation bias towards AT, and j/n the frequency class j.
#' Tws and Tsw can be expressed as a function of the Observed SFSs (Ows and Osw) and error rates (e1 and e2):
#' Tws = Ows(1-e2) - rev(Osw)e2 and Tsw = Osw(1-e1) - rev(Ows)*e1
#' The function return the weighted least-square as afunction of: sum(w(j) * (log(Tws(j)/Tsw(j)) - B * j/n + log(mut_bias) )^2)
#' where Tws and Tsw are expressed as a function of Ows, Osw, e1 and e2
#' w(j) = Ows(j)*Osw(j) / (Ows(j)+Osw(j)) is the weight used in the least square
#'
#' @param par a vector with the four parameters of the model.
#' par(1) = B
#' par(2) = M (i.e. log(mut_bias))
#' par(3) = e1
#' par(4) = e2
#' @param WS the WS observed SFS
#' @param SW the SW observed SFS
#' @param GC the GC content
#' @param cor a Boolean to add a correction to the transformed SFS (default = TRUE)
#'
#' @returns The weighted sum of squares
#'
#' @export
#'
#' @examples
#' sfsWS <- c(100,50,30,15,10)
#' sfsSW <- c(200,80, 30, 10, 5)
#' param <- c(1,2,0.02,0.01)
#' sum_of_squares_BM_err(param,sfsWS,sfsSW,0.5)
sum_of_squares_BM_err <- function(par,WS,SW,GC,cor=COR) {
  #Error checking
  if(!is.numeric(c(WS,SW)) || length(which(c(WS,SW)<0))>0 ) {
    abort("The two SFSs must have positive numeric values")
  }
  if(length(WS)!=length(SW)) {
    abort("The two SFSs, WS and SW, must have the same length")
  }
  if(length(par)!=4) {
    abort("A four values vector must be given as par")
  }
  if(GC<=0 | GC>=1) {
    abort("GC content must be strictly between 0 and 1")
  }
  #Main
  B <- par[1]
  M <- par[2]
  e1 <- par[3]
  e2 <- par[4]
  n <- length(WS)
  # True SFS as a function of observed one.
  WSt <- ((1 - e2)*WS - e2*rev(SW))/(1 - e1 - e2)
  SWt <- ((1 - e1)*SW - e1*rev(WS))/(1 - e1 - e2)
  removeNA <- which(WSt>=1 & SWt>=1)
  WSt <- WSt[removeNA]
  SWt <- SWt[removeNA]
  w <- 1/(t_variance(WSt,cor) + t_variance(SWt,cor))
  x <- c(1:n)/(n+1)
  x <- x[removeNA]
  y <- t_sfs(WSt,cor) - t_sfs(SWt,cor)
  ypred <- B*x - M - log(GC) + log(1 - GC)
  return( sum(w*(y-ypred)^2)/sum(w) )
}


#' @title Gradient of Sum of squares of model BM with errors
#'
#' @description Gradient of the sum of squares function
#'
#' @param par a vector with the four parameters of the model.
#' par(1) = B
#' par(2) = M (i.e. log(mut_bias))
#' par(3) = e1
#' par(4) = e2
#' @param WS the WS observed SFS
#' @param SW the SW observed SFS
#' @param GC the GC content
#' @param cor a Boolean to add a correction to the transformed SFS (default = TRUE)
#'
#' @returns The gradient function
#'
#' @noRd
gr_sum_of_squares_BM_err <- function(par,WS,SW,GC,cor=COR) {
  #Error checking
  if(!is.numeric(c(WS,SW)) || length(which(c(WS,SW)<0))>0 ) {
    abort("The two SFSs must have positive numeric values")
  }
  if(length(WS)!=length(SW)) {
    abort("The two SFSs, WS and SW, must have the same length")
  }
  if(length(par)!=4) {
    abort("A four values vector must be given as par")
  }
  if(GC<=0 | GC>=1) {
    abort("GC content must be strictly between 0 and 1")
  }
  #Main
  B <- par[1]
  M <- par[2]
  e1 <- par[3]
  e2 <- par[4]
  n <- length(WS)
  # True SFS as a function of observed one.
  WSt <- ((1-e2)*WS - e2*rev(SW))/(1 - e1 - e2)
  SWt <- ((1-e1)*SW - e1*rev(WS))/(1 - e1 - e2)
  removeNA <- which(WSt>=1 & SWt>=1)
  WSt <- WSt[removeNA]
  SWt <- SWt[removeNA]
  w <- 1/(t_variance(WSt,cor) + t_variance(SWt,cor))
  x <- c(1:n)/(n+1)
  x <- x[removeNA]
  y <- t_sfs(WSt,cor) - t_sfs(SWt,cor)
  ypred <- B*x - M - log(GC) + log(1 - GC)
  # Derivative of y as a function of error rates
  # Approximated version, derivative of log(WSt2/SWt2) without additional terms
  dy1 <- d_expected_log_ratio(WS,SW,e1,e2)$d1
  dy2 <- d_expected_log_ratio(WS,SW,e1,e2)$d2
  grB <- sum(w*(-2*x*(y - ypred))/sum(w))
  grM <- sum(w*(2*(y - ypred))/sum(w))
  gre1 <- sum(w*(2*(y - ypred)*dy1)/sum(w))
  gre2 <- sum(w*(2*(y - ypred)*dy2)/sum(w))
  return(c(grB,grM,gre1,gre2))
}
