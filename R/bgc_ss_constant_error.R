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
sum_of_squares_M_err <- function(par,WS,SW,GC) {
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
  # w <- WSt2*SWt2/(WSt2 + SWt2)
  # Here the weight is written as a function of corrected SFSs
  # This complexifies the whole equation as e1 and e2 appear in w
  # The gradient function is also more complicated
  # Instead we use the weight as a function of observed SFSs
  w <- 1/(t_variance(WS) + t_variance(SW))
  y <- t_sfs(WSt) - t_sfs(SWt)
  # y <- log(WSt2/SWt2) + 1/(WSt) - 1/(SWt)
  ypred <- rep(- M - log(GC) + log(1 - GC),n)
  removeNA <- which(WS!=0 & SW!=0)
  w <- w[removeNA]
  y <- y[removeNA]
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
#'
#' @returns The gradient function
#'
#' @noRd
gr_sum_of_squares_M_err <- function(par,WS,SW,GC) {
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
  w <- 1/(t_variance(WS) + t_variance(SW))
  y <- t_sfs(WS) - t_sfs(SW)
  ypred <- rep(- M - log(GC) + log(1 - GC),n)
  # Derivative of y as a function of error rates
  # Approximated version, derivative of log(WSt2/SWt2) without additional terms
  dy1 <- d_expected_log_ratio(WS,SW,e1,e2)$d1
  dy2 <- d_expected_log_ratio(WS,SW,e1,e2)$d2
  removeNA <- which(WS!=0 & SW!=0)
  w <- w[removeNA]
  y <- y[removeNA]
  ypred <- ypred[removeNA]
  dy1 <- dy1[removeNA]
  dy2 <- dy2[removeNA]
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
sum_of_squares_B_err <- function(par,WS,SW,GC) {
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
  # Same expression without 1-e1-e2 that simplifies in ratios
  WSt2 <- ((1 - e2)*WS - e2*rev(SW))
  SWt2 <- ((1 - e1)*SW - e1*rev(WS))
  w <- WSt2*SWt2/(WSt2 + SWt2)
  # Here the weight is written as a function of corrected SFSs
  # This complexifies the whole equation as e1 and e2 appear in w
  # The gradient function is also more complicated
  # Instead we use the weight as a function of observed SFSs
  #w <- 1/(t_variance(WS) + t_variance(SW))
  x <- c(1:n)/(n+1)
  #y <- t_sfs(WS) - t_sfs(SW)
  y <- log(WSt2/SWt2) + 1/(WSt) - 1/(SWt)
  ypred <- B*x - log(GC) + log(1 - GC)
  removeNA <- which(WS!=0 & SW!=0)
  w <- w[removeNA]
  y <- y[removeNA]
  ypred <- ypred[removeNA]
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
#'
#' @returns The gradient function
#'
#' @noRd
gr_sum_of_squares_B_err <- function(par,WS,SW,GC) {
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
  # Same expression without 1-e1-e2 that simplifies in ratios
  #WSt2 <- ((1 - e2)*WS - e2*rev(SW))
  #SWt2 <- ((1 - e1)*SW - e1*rev(WS))
  w <- 1/(t_variance(WS) + t_variance(SW))
  x <- c(1:n)/(n+1)
  y <- t_sfs(WS) - t_sfs(SW)
  #y <- log(WSt2/SWt2) + 1/(WSt) - 1/(SWt)
  ypred <- B*x - log(GC) + log(1 - GC)
  # Derivative of y as a function of error rates
  # Approximated version, derivative of log(WSt2/SWt2) without additional terms
  dy1 <- d_expected_log_ratio(WS,SW,e1,e2)$d1
  dy2 <- d_expected_log_ratio(WS,SW,e1,e2)$d2
  removeNA <- which(WS!=0 & SW!=0)
  w <- w[removeNA]
  x <- x[removeNA]
  y <- y[removeNA]
  ypred <- ypred[removeNA]
  dy1 <- dy1[removeNA]
  dy2 <- dy2[removeNA]
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
sum_of_squares_BM_err <- function(par,WS,SW,GC) {
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
  # Same expression without 1-e1-e2 that simplifies in ratios
  WSt2 <- ((1 - e2)*WS - e2*rev(SW))
  SWt2 <- ((1 - e1)*SW - e1*rev(WS))
  w <- 1/(t_variance(WSt) + t_variance(SWt))
  x <- c(1:n)/(n+1)
  y <- t_sfs(WSt) - t_sfs(SWt)
  #w <- WSt2*SWt2/(WSt2 + SWt2)
  # Here the weight is written as a function of corrected SFSs
  # This complexifies the whole equation as e1 and e2 appear in w
  # The gradient function is also more complicated
  # Instead we use the weight as a function of observed SFSs
  # w <- WS*SW/(WS + SW)
  #w <- 1/( (-1+3*WS-5*WS^2+4*WS^3)/(4*WS^4) + (-1+3*SW-5*SW^2+4*SW^3)/(4*SW^4) )
  #w <- 1/((1+2*WS)^2/(4*WS*(1+WS)^2) + (1+2*SW)^2/(4*SW*(1+SW)^2))
  #w <- 4/(log((2 + 6*WS + WS^2)/(WS*(2 + WS))) + log((2 + 6*SW + SW^2)/(SW*(2 + SW))))
  #w <- 1/(2*Log((1 + 1/SW)*(1 + 1/WS)) -
  # Log(1 + 1/(SW*WS) + Log((1 + 1/SW)*(1 + 1/WS))))
  #w <- 1/(
  #  (-1 + 16*WS^2 + 44*WS^3 + 44*WS^4 + 16*WS^5)/((16*WS^2)* (1 + WS)^4) +
  #     (-1 + 16*SW^2 + 44*SW^3 + 44*SW^4 + 16*SW^5)/((16*SW^2)* (1 + SW)^4)
  #)
  #y <- log(WSt2/SWt2) - 1/(2*WSt) + 1/(2*SWt)
  #y <- log(WSt2/SWt2) -0.05/WSt + 4.6/(WSt^2) -5.4/(WSt^3) +0.05/SWt - 4.6/(SWt^2) +5.4/(SWt^3)
  #y <- log(WS/SW) # - 1/(2*WS) + 1/(2*SW)
  #y <- (log(WS/SW) + log((WS+1)/(SW+1)))/2
  ypred <- B*x - M - log(GC) + log(1 - GC)
  removeNA <- which(WS>0 & SW>0)
  #removeNA <- which(WS>1 & SW>1)
  w <- w[removeNA]
  y <- y[removeNA]
  ypred <- ypred[removeNA]
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
#'
#' @returns The gradient function
#'
#' @noRd
gr_sum_of_squares_BM_err <- function(par,WS,SW,GC) {
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
  # Same expression without 1-e1-e2 that simplifies in ratios
  #WSt2 <- ((1 - e2)*WS - e2*rev(SW))
  #SWt2 <- ((1 - e1)*SW - e1*rev(WS))
  w <- 1/(t_variance(WS) + t_variance(SW))
  x <- c(1:n)/(n+1)
  y <- t_sfs(WS) - t_sfs(SW)
  ypred <- B*x - M - log(GC) + log(1 - GC)
  # Derivative of y as a function of error rates
  # Approximated version, derivative of log(WSt2/SWt2) without additional terms
  dy1 <- d_expected_log_ratio(WS,SW,e1,e2)$d1
  dy2 <- d_expected_log_ratio(WS,SW,e1,e2)$d2
  removeNA <- which(WS!=0 & SW!=0)
  w <- w[removeNA]
  x <- x[removeNA]
  y <- y[removeNA]
  ypred <- ypred[removeNA]
  dy1 <- dy1[removeNA]
  dy2 <- dy2[removeNA]
  grB <- sum(w*(-2*x*(y - ypred))/sum(w))
  grM <- sum(w*(2*(y - ypred))/sum(w))
  gre1 <- sum(w*(2*(y - ypred)*dy1)/sum(w))
  gre2 <- sum(w*(2*(y - ypred)*dy2)/sum(w))
  return(c(grB,grM,gre1,gre2))
}
