# Estimation of selection and GC-biased gene conversion intensity from site frequency spectrum data by a least-square approach
# Sylvain Glemin
# sylvain.glemin@univ-rennes.fr


# Series of SS functions for the different null model




########################################## #
# NULL MODELS ##############################
########################################## #

# This model corresponds to a fit of the intercept only
# It gives the mean and the total sum of square of the dataset
# The totla sum of squares is used for R2 calculations of other models


#' @title Sum of squares of the null model
#'
#' @description Function that returns the weighted total sum of squares for the model with only the intercept
#'
#' @param WS the WS observed SFS
#' @param SW the SW observed SFS
#' @param cor a Boolean to add a correction to the transformed SFS (default = TRUE)
#' @param snpthresh the value below which the SNP category is removed (default = 1)
#'
#' @returns The weighted sum of squares
#'
#' @export
#'
#' @examples
#' sfsWS <- c(100,50,30,15,10)
#' sfsSW <- c(200,80, 30, 10, 5)
#' sum_of_squares_NULL(sfsWS,sfsSW)
sum_of_squares_NULL <- function(WS,SW,cor=COR,snpthresh=SNPTHRESH) {
  #Error checking
  if(!is.numeric(c(WS,SW)) || length(which(c(WS,SW)<0))>0 ) {
    abort("The two SFSs must have positive numeric values")
  }
  if(length(WS)!=length(SW)) {
    abort("The two SFSs, WS and SW, must have the same length")
  }
  #Main
  n <- length(WS)
  # Suppressing values below the threshold
  removeNA <- which(WS>=snpthresh & SW>=snpthresh)
  WS <- WS[removeNA]
  SW <- SW[removeNA]
  w <- 1/(t_variance(WS,cor) + t_variance(SW,cor))
  y <- t_sfs(WS,cor) - t_sfs(SW,cor)
  y_mean <- mean(y)
  return(
    list(
      "mean"=y_mean,
      "SStot"=sum(w*(y - y_mean)^2)/sum(w)
    )
  )
}


# This model corresponds to a fit of the intercept fixed to the predicted effect of GC alone:
# intercept = log(1-GC) -log(GC)
# This is used to test whether the mutation bias is significant or not


#' @title Sum of squares of the GC null model
#'
#' @description Function that returns the weighted total sum of squares for the model with only the intercept
#'
#' @param WS the WS observed SFS
#' @param SW the SW observed SFS
#' @param GC GC content
#' @param cor a Boolean to add a correction to the transformed SFS (default = TRUE)
#' @param snpthresh the value below which the SNP category is removed (default = 1)
#'
#' @returns The weighted sum of squares
#'
#' @export
#'
#' @examples
#' sfsWS <- c(100,50,30,15,10)
#' sfsSW <- c(200,80, 30, 10, 5)
#' sum_of_squares_NULL_GC(sfsWS,sfsSW,0.5)
sum_of_squares_NULL_GC <- function(WS,SW,GC,cor=COR,snpthresh=SNPTHRESH) {
  #Error checking
  if(!is.numeric(c(WS,SW)) || length(which(c(WS,SW)<0))>0 ) {
    abort("The two SFSs must have positive numeric values")
  }
  if(length(WS)!=length(SW)) {
    abort("The two SFSs, WS and SW, must have the same length")
  }
  if(GC<=0 | GC>=1) {
    abort("GC content must be strictly between 0 and 1")
  }
  #Main
  n <- length(WS)
  # Suppressing values below the threshold
  removeNA <- which(WS>=snpthresh & SW>=snpthresh)
  WS <- WS[removeNA]
  SW <- SW[removeNA]
  w <- 1/(t_variance(WS,cor) + t_variance(SW,cor))
  y <- t_sfs(WS,cor) - t_sfs(SW,cor)
  y_pred <- log(1-GC) - log(GC)
  SStot <- sum_of_squares_NULL(SW,WS)$SStot
  SSres <- sum(w*(y - y_pred)^2)/sum(w)
  R2 <- 1 - SSres/SStot
  AIC <- AICls(length(which(WS!=0 & SW!=0)),0,SSres)
  return(
    list("SStot"=SStot,
         "SSres"=SSres,
         "R2"=R2,
         "AIC"=AIC
    )
  )
}


######################################### #
# NULL MODEL WITH ERROR ###################
######################################### #

# This model corresponds to a fit of the intercept fixed to the predicted effect of GC alone:
# intercept = log(1-GC) -log(GC)
# This is used to test whether the mutation bias is significant for model with errors


#' @title Sum of squares of the GC null model with errors
#'
#' @description Function that return the weighted sum of squares between the function of the SFSs and error rates and the linear predictor.
#' The theory predict that the expectation of log(Tws(j)/Tsw(j)) = -log(mutbias) + log(GC) - log(1 - GC)
#' where Tws and Tsw are the "True" SFSs and j/n the frequency class j.
#' Tws and Tsw can be expressed as a function of the Observed SFSs (Ows and Osw) and error rates (e1 and e2):
#' Tws = Ows(1-e2) - rev(Osw)e2 and Tsw = Osw(1-e1) - rev(Ows)*e1
#' The function return the weighted least-square as afunction of: sum(w(j) * (log(Tws(j)/Tsw(j)) - B * j/n + log(mut_bias) )^2)
#' where Tws and Tsw are expressed as a function of Ows, Osw, e1 and e2
#' w(j) = Ows(j)*Osw(j) / (Ows(j)+Osw(j)) is the weight used in the least square
#' Note that this function must be optimized for e1 and e2 whereas the other null model directly give the sum of squares
#'
#'
#'
#' @param par a vector with the for parameters of the model.
#' par(1) = e1
#' par(2) = e2
#' @param WS the WS observed SFS
#' @param SW the SW observed SFS
#' @param GC GC content
#' @param cor a Boolean to add a correction to the transformed SFS (default = TRUE)
#' @param snpthresh the value below which the SNP category is removed (default = 1)
#'
#' @returns The weighted sum of squares
#'
#' @export
#'
#' @examples
#' sfsWS <- c(100,50,30,15,10)
#' sfsSW <- c(200,80,30,10,5)
#' param <- c(0.02,0.01)
#' sum_of_squares_NULL_GC_err(param,sfsWS,sfsSW,0.5)
sum_of_squares_NULL_GC_err <- function(par,WS,SW,GC,cor=COR,snpthresh=SNPTHRESH) {
  #Error checking
  if(!is.numeric(c(WS,SW)) || length(which(c(WS,SW)<0))>0 ) {
    abort("The two SFSs must have positive numeric values")
  }
  if(length(WS)!=length(SW)) {
    abort("The two SFSs, WS and SW, must have the same length")
  }
  if(length(par)!=2) {
    abort("A two values vector must be given as par")
  }
  if(GC<=0 | GC>=1) {
    abort("GC content must be strictly between 0 and 1")
  }
  #Main
  e1 <- par[1]
  e2 <- par[2]
  n <- length(WS)
  # True SFS as a function of observed one.
  WSt <- ((1 - e2)*WS - e2*rev(SW))/(1 - e1 - e2)
  SWt <- ((1 - e1)*SW - e1*rev(WS))/(1 - e1 - e2)
  # Suppressing values below the threshold
  removeNA <- which(WSt>=snpthresh & SWt>=snpthresh & WS>=snpthresh & SW>=snpthresh)
  WSt <- WSt[removeNA]
  SWt <- SWt[removeNA]
  # Variables for the regression
  w <- 1/(t_variance(WSt,cor) + t_variance(SWt,cor))
  x <- c(1:n)/(n+1)
  x <- x[removeNA]
  y <- t_sfs(WSt,cor) - t_sfs(SWt,cor)
  ypred <- rep(log(1 - GC) - log(GC),n)
  ypred <- ypred[removeNA]
  SS <- sum(w*(y-ypred)^2)/sum(w)
  if(is.na(SS)) {
    return(sum_of_squares_NULL(SW,WS)$SStot)
  } else {
    return(SS)
  }
}


