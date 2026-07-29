# Estimation of selection and GC-biased gene conversion intensity from site frequency spectrum data by a least-square approach
# Sylvain Glemin
# sylvain.glemin@univ-rennes.fr

# Series of SS functions and their corresponding gradient for gBGC models with hotspots


######################################### #
# SIMPLE HOTSPOT MODEL HOTSPOT 1 ##########
######################################### #

#' @title Sum of squares of the hotspot 1 model: no background gBGC + hotspot
#'
#' @description Function that return the weighted sum of squares between the function of the SFSs and error rates and the predictor.
#' The theory predict the expectation of log(Tws(j)/Tsw(j))
#' where Tws and Tsw are the "True" SFSs
#' Tws and Tsw can be expressed as a function of the Observed SFSs (Ows and Osw) and error rates (e1 and e2):
#' Tws = Ows(1-e2) - rev(Osw)e2 and Tsw = Osw(1-e1) - rev(Ows)*e1
#' The function return the weighted least-square as afunction of: sum(w(j) * (log(Tws(j)/Tsw(j)) - B * j/n + log(mut_bias) )^2)
#' where Tws and Tsw are expressed as a function of Ows, Osw, e1 and e2
#' w(j) = Ows(j)*Osw(j) / (Ows(j)+Osw(j)) is the weight used in the least square
#'
#' @param par a vector with the six parameters of the model.
#' par(1) = B (gBGC background)
#' par(2) = f (proportion of hotspots)
#' par(3) = M (i.e. log(mut_bias))
#' par(4) = e1
#' par(5) = e2
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
#' param <- c(1,0.1,2,0.02,0.01)
#' sum_of_squares_hotspot1_err(param,sfsWS,sfsSW,0.5)
sum_of_squares_hotspot1_err <- function(par,WS,SW,GC,cor=COR) {
  #Error checking
  if(!is.numeric(c(WS,SW)) || length(which(c(WS,SW)<0))>0 ) {
    abort("The two SFSs must have positive numeric values")
  }
  if(length(par)!=5) {
    abort("A five values vector must be given as par")
  }
  if(length(WS)!=length(SW)) {
    abort("The two SFSs, WS and SW, must have the same length")
  }
  if(GC<=0 | GC>=1) {
    abort("GC content must be strictly between 0 and 1")
  }
  #Main
  B <- par[1]
  f <- par[2]
  M <- par[3]
  e1 <- par[4]
  e2 <- par[5]
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
  ypred <- - M - log(GC) + log(1 - GC) +
    log((1 - f)*ratio_B(0,x) + f*ratio_B(B,x)) -
    log((1 - f)*ratio_B(0,x) + f*ratio_B(-B,x))
  return( sum(w*(y-ypred)^2)/sum(w) )
}



#' @title Gradient of Sum of squares of model hotspot2
#'
#' @description Gradient of the sum of squares function
#'
#' @param par a vector with the six parameters of the model.
#' par(1) = B
#' par(2) = f
#' par(3) = M (i.e. log(mut_bias))
#' par(4) = e1
#' par(5) = e2
#' @param WS the WS observed SFS
#' @param SW the SW observed SFS
#' @param GC the GC content
#' @param cor a Boolean to add a correction to the transformed SFS (default = TRUE)
#'
#' @returns The gradient function
#'
#' @noRd
gr_sum_of_squares_hotspot1_err <- function(par,WS,SW,GC,cor=COR) {
  #Error checking
  if(!is.numeric(c(WS,SW)) || length(which(c(WS,SW)<0))>0 ) {
    abort("The two SFSs must have positive numeric values")
  }
  if(length(par)!=5) {
    abort("A five values vector must be given as par")
  }
  if(length(WS)!=length(SW)) {
    abort("The two SFSs, WS and SW, must have the same length")
  }
  if(GC<=0 | GC>=1) {
    abort("GC content must be strictly between 0 and 1")
  }
  #Main
  B <- par[1]
  f <- par[2]
  M <- par[3]
  e1 <- par[4]
  e2 <- par[5]
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
  ypred <- - M - log(GC) + log(1 - GC) +
    log((1 - f)*ratio_B(0,x) + f*ratio_B(B,x)) -
    log((1 - f)*ratio_B(0,x) + f*ratio_B(-B,x))
  # Derivative of the log ratio of expectation
  dB <- d_hotspot1(B,f,x)$dB
  df <- d_hotspot1(B,f,x)$df
  # Derivative of y as a function of error rates
  # Approximated version, derivative of log(WSt2/SWt2) without additional terms
  dy1 <- d_expected_log_ratio(WS,SW,e1,e2)$d1
  dy2 <- d_expected_log_ratio(WS,SW,e1,e2)$d2
  grB <- sum(w*(2*(y - ypred)*dB)/sum(w))
  grf <- sum(w*(2*(y - ypred)*df)/sum(w))
  grM <- sum(w*(2*(y - ypred))/sum(w))
  gre1 <- sum(w*(2*(y - ypred)*dy1)/sum(w))
  gre2 <- sum(w*(2*(y - ypred)*dy2)/sum(w))
  return(c(grB,grf,grM,gre1,gre2))
}



######################################### #
# FULL HOTSPOT MODEL HOTSPOT 2 ############
######################################### #

#' @title Sum of squares of the hotspot 2 model: background gBGC + hotspot
#'
#' @description Function that return the weighted sum of squares between the function of the SFSs and error rates and the predictor.
#' The theory predict the expectation of log(Tws(j)/Tsw(j))
#' where Tws and Tsw are the "True" SFSs
#' Tws and Tsw can be expressed as a function of the Observed SFSs (Ows and Osw) and error rates (e1 and e2):
#' Tws = Ows(1-e2) - rev(Osw)e2 and Tsw = Osw(1-e1) - rev(Ows)*e1
#' The function return the weighted least-square as afunction of: sum(w(j) * (log(Tws(j)/Tsw(j)) - B * j/n + log(mut_bias) )^2)
#' where Tws and Tsw are expressed as a function of Ows, Osw, e1 and e2
#' w(j) = Ows(j)*Osw(j) / (Ows(j)+Osw(j)) is the weight used in the least square
#'
#' @param par a vector with the six parameters of the model.
#' par(1) = B0 (gBGC background)
#' par(2) = B1 (gBGC in hotspots)
#' par(3) = f (proportion of hotspots)
#' par(4) = M (i.e. log(mut_bias))
#' par(5) = e1
#' par(6) = e2
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
#' param <- c(0.1,5,0.1,2,0.02,0.01)
#' sum_of_squares_hotspot2_err(param,sfsWS,sfsSW,0.5)
sum_of_squares_hotspot2_err <- function(par,WS,SW,GC,cor=COR) {
  #Error checking
  if(!is.numeric(c(WS,SW)) || length(which(c(WS,SW)<0))>0 ) {
    abort("The two SFSs must have positive numeric values")
  }
  if(length(par)!=6) {
    abort("A six values vector must be given as par")
  }
  if(length(WS)!=length(SW)) {
    abort("The two SFSs, WS and SW, must have the same length")
  }
  if(GC<=0 | GC>=1) {
    abort("GC content must be strictly between 0 and 1")
  }
  #Main
  B0 <- par[1]
  B1 <- par[2]
  f <- par[3]
  M <- par[4]
  e1 <- par[5]
  e2 <- par[6]
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
  ypred <- - M - log(GC) + log(1 - GC) +
    log((1 - f)*ratio_B(B0,x) + f*ratio_B(B1,x)) -
    log((1 - f)*ratio_B(-B0,x) + f*ratio_B(-B1,x))
  removeNA <- !is.na(y)
  w <- w[removeNA]
  y <- y[removeNA]
  ypred <- ypred[removeNA]
  return( sum(w*(y-ypred)^2)/sum(w) )
}



#' @title Gradient of Sum of squares of model hotspot2
#'
#' @description Gradient of the sum of squares function
#'
#' @param par a vector with the six parameters of the model.
#' par(1) = B0
#' par(2) = B1
#' par(3) = f
#' par(4) = M (i.e. log(mut_bias))
#' par(5) = e1
#' par(6) = e2
#' @param WS the WS observed SFS
#' @param SW the SW observed SFS
#' @param GC the GC content
#' @param cor a Boolean to add a correction to the transformed SFS (default = TRUE)
#'
#' @returns The gradient function
#'
#' @noRd
gr_sum_of_squares_hotspot2_err <- function(par,WS,SW,GC,cor=COR) {
  #Error checking
  if(!is.numeric(c(WS,SW)) || length(which(c(WS,SW)<0))>0 ) {
    abort("The two SFSs must have positive numeric values")
  }
  if(length(par)!=6) {
    abort("A six values vector must be given as par")
  }
  if(length(WS)!=length(SW)) {
    abort("The two SFSs, WS and SW, must have the same length")
  }
  if(GC<=0 | GC>=1) {
    abort("GC content must be strictly between 0 and 1")
  }
  #Main
  B0 <- par[1]
  B1 <- par[2]
  f <- par[3]
  M <- par[4]
  e1 <- par[5]
  e2 <- par[6]
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
  ypred <- - M - log(GC) + log(1 - GC) +
    log((1 - f)*ratio_B(B0,x) + f*ratio_B(B1,x)) -
    log((1 - f)*ratio_B(-B0,x) + f*ratio_B(-B1,x))
  # Derivative of the log ratio of expectation
  dB0 <- d_hotspot2(B0,B1,f,x)$dB0
  dB1 <- d_hotspot2(B0,B1,f,x)$dB1
  df <- d_hotspot2(B0,B1,f,x)$df
  # Derivative of y as a function of error rates
  # Approximated version, derivative of log(WSt2/SWt2) without additional terms
  dy1 <- d_expected_log_ratio(WS,SW,e1,e2)$d1
  dy2 <- d_expected_log_ratio(WS,SW,e1,e2)$d2
  grB0 <- sum(w*(2*(y - ypred)*dB0)/sum(w))
  grB1 <- sum(w*(2*(y - ypred)*dB1)/sum(w))
  grf <- sum(w*(2*(y - ypred)*df)/sum(w))
  grM <- sum(w*(2*(y - ypred))/sum(w))
  gre1 <- sum(w*(2*(y - ypred)*dy1)/sum(w))
  gre2 <- sum(w*(2*(y - ypred)*dy2)/sum(w))
  return(c(grB0,grB1,grf,grM,gre1,gre2))
}


############################################## #
# FIXED FRACTION OF HOTSPOT MODEL HOTSPOT 2BIS #
############################################## #

#' @title Sum of squares of the hotspot 2bis model: background gBGC + fixed proportion of hotspots
#'
#' @description Function that return the weighted sum of squares between the function of the SFSs and error rates and the predictor.
#' The theory predict the expectation of log(Tws(j)/Tsw(j))
#' where Tws and Tsw are the "True" SFSs
#' Tws and Tsw can be expressed as a function of the Observed SFSs (Ows and Osw) and error rates (e1 and e2):
#' Tws = Ows(1-e2) - rev(Osw)e2 and Tsw = Osw(1-e1) - rev(Ows)*e1
#' The function return the weighted least-square as afunction of: sum(w(j) * (log(Tws(j)/Tsw(j)) - B * j/n + log(mut_bias) )^2)
#' where Tws and Tsw are expressed as a function of Ows, Osw, e1 and e2
#' w(j) = Ows(j)*Osw(j) / (Ows(j)+Osw(j)) is the weight used in the least square
#'
#' @param par a vector with the five parameters of the model.
#' par(1) = B0 (gBGC background)
#' par(2) = B1 (gBGC in hotspots)
#' par(3) = M (i.e. log(mut_bias))
#' par(4) = e1
#' par(5) = e2
#' @param WS the WS observed SFS
#' @param SW the SW observed SFS
#' @param GC the GC content
#' @param f fraction of hotspots, fixed, not estimated (0<f<1/2)
#' @param cor a Boolean to add a correction to the transformed SFS (default = TRUE)
#'
#' @returns The weighted sum of squares
#'
#' @export
#'
#' @examples
#' sfsWS <- c(100,50,30,15,10)
#' sfsSW <- c(200,80, 30, 10, 5)
#' param <- c(0.1,5,2,0.02,0.01)
#' sum_of_squares_hotspot2bis_err(param,sfsWS,sfsSW,0.5,0.1)
sum_of_squares_hotspot2bis_err <- function(par,WS,SW,GC,f,cor=COR) {
  #Error checking
  if(!is.numeric(c(WS,SW)) || length(which(c(WS,SW)<0))>0 ) {
    abort("The two SFSs must have positive numeric values")
  }
  if(length(par)!=5) {
    abort("A five values vector must be given as par")
  }
  if(length(WS)!=length(SW)) {
    abort("The two SFSs, WS and SW, must have the same length")
  }
  if(GC<=0 | GC>=1) {
    abort("GC content must be strictly between 0 and 1")
  }
  if(f<=0 | f>=1/2) {
    abort("f must be strictly between 0 and 1/2")
  }
  #Main
  B0 <- par[1]
  B1 <- par[2]
  M <- par[3]
  e1 <- par[4]
  e2 <- par[5]
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
  ypred <- - M - log(GC) + log(1 - GC) +
    log((1 - f)*ratio_B(B0,x) + f*ratio_B(B1,x)) -
    log((1 - f)*ratio_B(-B0,x) + f*ratio_B(-B1,x))
  return( sum(w*(y-ypred)^2)/sum(w) )
}



#' @title Gradient of Sum of squares of model hotspot2bis
#'
#' @description Gradient of the sum of squares function
#'
#' @param par a vector with the five parameters of the model.
#' par(1) = B0
#' par(2) = B1
#' par(3) = M (i.e. log(mut_bias))
#' par(4) = e1
#' par(5) = e2
#' @param WS the WS observed SFS
#' @param SW the SW observed SFS
#' @param GC the GC content
#' @param f fraction of hotpsots
#' @param cor a Boolean to add a correction to the transformed SFS (default = TRUE)
#'
#' @returns The gradient function
#'
#' @noRd
gr_sum_of_squares_hotspot2bis_err <- function(par,WS,SW,GC,f,cor=COR) {
  #Error checking
  if(!is.numeric(c(WS,SW)) || length(which(c(WS,SW)<0))>0 ) {
    abort("The two SFSs must have positive numeric values")
  }
  if(length(par)!=5) {
    abort("A five values vector must be given as par")
  }
  if(length(WS)!=length(SW)) {
    abort("The two SFSs, WS and SW, must have the same length")
  }
  if(GC<=0 | GC>=1) {
    abort("GC content must be strictly between 0 and 1")
  }
  if(f<=0 | f>=1/2) {
    abort("f must be strictly between 0 and 1/2")
  }
  #Main
  B0 <- par[1]
  B1 <- par[2]
  M <- par[3]
  e1 <- par[4]
  e2 <- par[5]
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
  ypred <- - M - log(GC) + log(1 - GC) +
    log((1 - f)*ratio_B(B0,x) + f*ratio_B(B1,x)) -
    log((1 - f)*ratio_B(-B0,x) + f*ratio_B(-B1,x))
  # Derivative of the log ratio of expectation
  dB0 <- d_hotspot2(B0,B1,f,x)$dB0
  dB1 <- d_hotspot2(B0,B1,f,x)$dB1
  # Derivative of y as a function of error rates
  # Approximated version, derivative of log(WSt2/SWt2) without additional terms
  dy1 <- d_expected_log_ratio(WS,SW,e1,e2)$d1
  dy2 <- d_expected_log_ratio(WS,SW,e1,e2)$d2
  grB0 <- sum(w*(2*(y - ypred)*dB0)/sum(w))
  grB1 <- sum(w*(2*(y - ypred)*dB1)/sum(w))
  grM <- sum(w*(2*(y - ypred))/sum(w))
  gre1 <- sum(w*(2*(y - ypred)*dy1)/sum(w))
  gre2 <- sum(w*(2*(y - ypred)*dy2)/sum(w))
  return(c(grB0,grB1,grM,gre1,gre2))
}
