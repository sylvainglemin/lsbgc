# Estimation of selection and GC-biased gene conversion intensity from site frequency spectrum data by a least-square approach
# Sylvain Glemin
# sylvain.glemin@univ-rennes.fr




# Functions to minimize the sum of squares of the different models


#' @title Sum of squares minimization of model M with error
#'
#' @description Function that searches for the three parameters that minimize the sum_of_squares function
#'
#' @param WS the WS observed SFS
#' @param SW the SW observed SFS
#' @param GC GC content
#' @param cor a Boolean to add a correction to the transformed SFS (default = TRUE)
#' @param snpthresh the value below which the SNP category is removed (default = 1)
#' @param Mmin minimum for the range of M, default value = -5
#' @param Mmax maximum for the range of M, default value = 5
#' @param e1max maximum for the range of WS error rate, default value = 0.49
#' @param e2max maximum for the range of SW error rate, default value = 0.49
#' @param e1init initial starting value for optimization for WS error rates, default = 0.01
#' @param e2init initial starting value for optimization for SW error rates, default = 0.01
#' @param Maxit maximum number of iterations (option for optim), see manual, default value = 100
#' @param Factr level of control the convergence (option for optim, see manual), default value = 10^7
#' @param Lmm number of updates in the method (option for optim, see manual)
#' @param Verbose from 0 (default) to 5: level of outputs during optimization (option for optim, see manual)
#' @param Usegr to use (default) or not the analytical gradient (option for optim, see manual)
#' It's better to use the analytical gradient function but the option can be turn off for testing
#'
#' @returns A list of two lists:
#' - param: Optimized parameters
#' - criteria: Model fit criteria
#'
#' @importFrom stats optim
#'
#' @export
#'
#' @examples
#' sfsWS <- c(1000,500,300,200,80,50)
#' sfsSW <- c(2000,800,400,150,50,10)
#' LS <- least_square_M_err(sfsWS,sfsSW,0.5)
#' LS$param$mutbias # mutation bias
#' LS$criteria$AIC # model AIC
least_square_M_err <- function(WS,SW,GC,cor=COR,snpthresh=SNPTHRESH,
                           Mmin=MMIN,Mmax=MMAX,e1max=EMAX,e2max=EMAX,e1init=EINIT,e2init=EINIT,
                           Maxit=MAXIT,Factr=FACTR,Lmm=LMM,Verbose=VERBOSE,Usegr=USEGR) {
  # Determination of initial values for optimization: use of the simple regression
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
  Minit <- mean(log(WS/SW),na.rm=T) + log(1 - GC) - log(GC)
  init <- c(Minit,e1init,e2init)
  # Boundaries for optimization
  inf <- c(Mmin,0,0)
  sup <- c(Mmax,e1max,e2max)
  if(Usegr) gradient <- gr_sum_of_squares_M_err else gradient <- NULL
  SCALE <- abs(init)
  minSSE <- optim(
    par = init,
    fn = sum_of_squares_M_err,
    gr = gradient,
    WS = WS, SW = SW, GC = GC, cor = cor, snpthresh = snpthresh,
    lower = inf,upper = sup,
    method = "L-BFGS-B",
    control=list(parscale=SCALE,maxit=Maxit,factr=Factr,lmm=Lmm,trace=Verbose))
  SStot <- sum_of_squares_NULL(SW,WS)$SStot
  SSres <- minSSE$value
  R2 <- 1 - SSres/SStot
  AIC <- AICls(length(which(WS!=0 & SW!=0)),length(init),SSres)
  return( list(
    "param"=list("mutbias"=exp(minSSE$par[1]),
                 "e1"=minSSE$par[2],
                 "e2"=minSSE$par[3]),
    "criteria"=list("SStot"=SStot,
                    "SSres"=SSres,
                    "R2"=R2,
                    "AIC"=AIC) )
  )
}



#' @title Sum of squares minimization of model B with error
#'
#' @description Function that searches for the four parameters that minimize the sum_of_squares function
#'
#' @param WS the WS observed SFS
#' @param SW the SW observed SFS
#' @param GC GC content
#' @param cor a Boolean to add a correction to the transformed SFS (default = TRUE)
#' @param snpthresh the value below which the SNP category is removed (default = 1)
#' @param Bmin minimum for the range of B, default value = -100
#' @param Bmax maximum for the range of B, default value = 100
#' @param e1max maximum for the range of WS error rate, default value = 0.49
#' @param e2max maximum for the range of SW error rate, default value = 0.49
#' @param e1init initial starting value for optimization for WS error rates, default = 0.01
#' @param e2init initial starting value for optimization for SW error rates, default = 0.01
#' @param Maxit maximum number of iterations (option for optim), see manual, default value = 100
#' @param Factr level of control the convergence (option for optim, see manual), default value = 10^7
#' @param Lmm number of updates in the method (option for optim, see manual)
#' @param Verbose from 0 (default) to 5: level of outputs during optimization (option for optim, see manual)
#' @param Usegr to use (default) or not the analytical gradient (option for optim, see manual)
#' It's better to use the analytical gradient function but the option can be turn off for testing
#'
#' @returns A list of two lists:
#' - param: Optimized parameters
#' - criteria: Model fit criteria
#'
#' @importFrom stats lm
#' @importFrom stats optim
#'
#' @export
#'
#' @examples
#' sfsWS <- c(1000,500,300,200,80,50)
#' sfsSW <- c(2000,800,400,150,50,10)
#' LS <- least_square_B_err(sfsWS,sfsSW,0.5)
#' LS$param$B # population-scaled gBGC
#' LS$criteria$R2 # R2 of the model
least_square_B_err <- function(WS,SW,GC,cor=COR,snpthresh=SNPTHRESH,
                           Bmin=BMIN,Bmax=BMAX,e1max=EMAX,e2max=EMAX,e1init=EINIT,e2init=EINIT,
                           Maxit=MAXIT,Factr=FACTR,Lmm=LMM,Verbose=VERBOSE,Usegr=USEGR) {
  # Determination of initial values for optimization: use of the simple regression
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
  x <- c(1:n)/(n+1)
  NONZERO <- which(WS!=0 & SW!=0)
  # Determination of initial values for optimization: use of the simple regression
  reginit <- lm(log(WS[NONZERO]/SW[NONZERO]) ~ x[NONZERO] - 1)
  Binit <- reginit$coef
  init <- c(Binit,e1init,e2init)
  # Boundaries for optimization
  inf <- c(Bmin,0,0)
  sup <- c(Bmax,e1max,e2max)
  if(Usegr) gradient <- gr_sum_of_squares_B_err else gradient <- NULL
  SCALE <- abs(init)
  minSSE <- optim(
    par = init,
    fn = sum_of_squares_B_err,
    gr = gradient,
    WS = WS, SW = SW, GC = GC, cor = cor, snpthresh = SNPTHRESH,
    lower = inf,upper = sup,
    method = "L-BFGS-B",
    control=list(parscale=SCALE,maxit=Maxit,factr=Factr,lmm=Lmm,trace=Verbose))
  SStot <- sum_of_squares_NULL(SW,WS)$SStot
  SSres <- minSSE$value
  R2 <- 1 - SSres/SStot
  AIC <- AICls(length(which(WS!=0 & SW!=0)),length(init),SSres)

  return( list(
    "param"=list("B"=minSSE$par[1],
                 "e1"=minSSE$par[2],
                 "e2"=minSSE$par[3]),
    "criteria"=list("SStot"=SStot,
                    "SSres"=SSres,
                    "R2"=R2,
                    "AIC"=AIC) )
  )
}


#' @title Sum of squares minimization of model BM with error
#'
#' @description Function that searches for the four parameters that minimize the sum_of_squares function
#'
#' @param WS the WS observed SFS
#' @param SW the SW observed SFS
#' @param GC GC content
#' @param cor a Boolean to add a correction to the transformed SFS (default = TRUE)
#' @param snpthresh the value below which the SNP category is removed (default = 1)
#' @param Bmin minimum for the range of B, default value = -100
#' @param Bmax maximum for the range of B, default value = 100
#' @param Mmin minimum for the range of M, default value = -5
#' @param Mmax maximum or the range of M, default value = 5
#' @param e1max maximum for the range of WS error rate, default value = 0.49
#' @param e2max maximum for the range of SW error rate, default value = 0.49
#' @param e1init initial starting value for optimization for WS error rates, default = 0.01
#' @param e2init initial starting value for optimization for SW error rates, default = 0.01
#' @param Maxit maximum number of iterations (option for optim), see manual, default value = 100
#' @param Factr level of control the convergence (option for optim, see manual), default value = 10^7
#' @param Lmm number of updates in the method (option for optim, see manual)
#' @param Verbose from 0 (default) to 5: level of outputs during optimization (option for optim, see manual)
#' @param Usegr to use (default) or not the analytical gradient (option for optim, see manual)
#' It's better to use the analytical gradient function but the option can be turn off for testing
#'
#' @returns A list of two lists:
#' - param: Optimized parameters
#' - criteria: Model fit criteria
#'
#' @importFrom stats lm
#' @importFrom stats optim
#'
#' @export
#'
#' @examples
#' sfsWS <- c(1000,500,300,200,80,50)
#' sfsSW <- c(2000,800,400,150,50,10)
#' LS <- least_square_BM_err(sfsWS,sfsSW,0.5)
#' LS$B # population-scaled gBGC
#' LS$mutbias # mutation bias
least_square_BM_err <- function(WS,SW,GC,cor=COR,snpthresh=SNPTHRESH,
                            Bmin=BMIN,Bmax=BMAX,Mmin=MMIN,Mmax=MMAX,e1max=EMAX,e2max=EMAX,e1init=EINIT,e2init=EINIT,
                            Maxit=MAXIT,Factr=FACTR,Lmm=LMM,Verbose=VERBOSE,Usegr=USEGR) {
  # Determination of initial values for optimization: use of the simple regression
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
  x <- c(1:n)/(n+1)
  NONZERO <- which(WS!=0 & SW!=0)
  # Determination of initial values for optimization: use of the simple regression
  reginit <- lm(log(WS[NONZERO]/SW[NONZERO]) ~ x[NONZERO])
  Minit <- -reginit$coef[1] + log(1 - GC) - log(GC)
  Binit <- reginit$coef[2]
  # To avoid negative values in SFS after transformation the error rates must be bounded as follows:
  # This solution is too strict empirically
  #e1max <- max(ZERO,min(SW/(SW + rev(WS)),na.rm=T) * ONE)
  #e2max <- max(ZERO,min(WS/(WS + rev(SW)),na.rm=T) * ONE)
  init <- c(Binit,Minit,e1init,e2init)
  # Boundaries for optimization
  inf <- c(Bmin,Mmin,0,0)
  sup <- c(Bmax,Mmax,e1max,e2max)
  if(Usegr) gradient <- gr_sum_of_squares_BM_err else gradient <- NULL
  SCALE <- abs(init)
  minSSE <- optim(
    par = init,
    fn = sum_of_squares_BM_err,
    gr = gradient,
    WS = WS, SW = SW, GC = GC, cor = cor, snpthresh = snpthresh,
    lower = inf,upper = sup,
    method = "L-BFGS-B",
    control=list(parscale=SCALE,maxit=Maxit,factr=Factr,lmm=Lmm,trace=Verbose)
  )
  SStot <- sum_of_squares_NULL(SW,WS)$SStot
  SSres <- minSSE$value
  R2 <- 1 - SSres/SStot
  AIC <- AICls(length(which(WS!=0 & SW!=0)),length(init),SSres)
  # Note that the number of observations is n-1
  return( list(
    "param"=list("B"=minSSE$par[1],
                 "mutbias"=exp(minSSE$par[2]),
                 "e1"=minSSE$par[3],
                 "e2"=minSSE$par[4]),
    "criteria"=list("SStot"=SStot,
                    "SSres"=SSres,
                    "R2"=R2,
                    "AIC"=AIC) )
  )
}


#' @title Sum of squares minimization of model Hotpsot 1 with error
#'
#' @description Function that searches for the six parameters that minimize the sum_of_squares function
#'
#' @param WS the WS observed SFS
#' @param SW the SW observed SFS
#' @param GC GC content
#' @param cor a Boolean to add a correction to the transformed SFS (default = TRUE)
#' @param snpthresh the value below which the SNP category is removed (default = 1)
#' @param Bmin minimum for the range of B0, default value = -100
#' @param Bmax maximum for the range of B0, default value = 100
#' @param Mmin minimum for the range of M, default value = -5
#' @param Mmax maximum or the range of M, default value = 5
#' @param e1max maximum for the range of WS error rate, default value = 0.49
#' @param e2max maximum for the range of SW error rate, default value = 0.49
#' @param e1init initial starting value for optimization for WS error rates, default = 0.01
#' @param e2init initial starting value for optimization for SW error rates, default = 0.01
#' @param Maxit maximum number of iterations (option for optim), see manual, default value = 100
#' @param Factr level of control the convergence (option for optim, see manual), default value = 10^7
#' @param Lmm number of updates in the method (option for optim, see manual)
#' @param Verbose from 0 (default) to 5: level of outputs during optimization (option for optim, see manual)
#' @param Usegr to use (default) or not the analytical gradient (option for optim, see manual)
#' It's better to use the analytical gradient function but the option can be turn off for testing
#'
#' @returns A list of two lists:
#' - param: Optimized parameters
#' - criteria: Model fit criteria
#'
#' @export
#'
#' @importFrom stats lm
#' @importFrom stats optim
#'
#' @examples
#' sfsWS <- c(1000,500,300,200,80,50)
#' sfsSW <- c(2000,800,400,150,50,10)
#' LS <- least_square_hotspot1_err(sfsWS,sfsSW,0.5)
#' LS$B # population-scaled hotspot gBGC
#' LS$f # proportion of hotspots
least_square_hotspot1_err <- function(WS,SW,GC,cor=COR,snpthresh=SNPTHRESH,
                                  Bmin=BMIN,Bmax=BMAX,Mmin=MMIN,Mmax=MMAX,e1max=EMAX,e2max=EMAX,e1init=EINIT,e2init=EINIT,
                                  Maxit=MAXIT,Factr=FACTR,Lmm=LMM,Verbose=VERBOSE,Usegr=USEGR) {
  # Determination of initial values for optimization: use of the simple regression
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
  x <- c(1:n)/(n+1)
  NONZERO <- which(WS!=0 & SW!=0)
  # Arbitray starting value fo f. Note that we assume 0 ≤ f ≤ 1/2
  finit <- 0.1
  # Determination of initial values for optimization: use of the simple regression
  reginit <- lm(log(WS[NONZERO]/SW[NONZERO]) ~ x[NONZERO])
  Minit <- -reginit$coef[1] + log(1 - GC) - log(GC)
  Binit <- reginit$coef[2]/finit
  init <- c(Binit,finit,Minit,e1init,e2init)
  # Boundaries for optimization
  inf <- c(Bmin,0,Mmin,0,0)
  sup <- c(Bmax,1/2,Mmax,e1max,e2max)
  if(Usegr) gradient <- gr_sum_of_squares_hotspot1_err else gradient <- NULL
  SCALE <- abs(init)
  minSSE <- optim(
    par = init,
    fn = sum_of_squares_hotspot1_err,
    gr = gradient,
    WS = WS, SW = SW, GC = GC, cor = cor, snpthresh = snpthresh,
    lower = inf,upper = sup,
    method = "L-BFGS-B",
    control=list(parscale=SCALE,maxit=Maxit,factr=Factr,lmm=Lmm,trace=Verbose))
  SStot <- sum_of_squares_NULL(SW,WS)$SStot
  SSres <- minSSE$value
  R2 <- 1 - SSres/SStot
  AIC <- AICls(length(which(WS!=0 & SW!=0)),length(init),SSres)
  # Note that the number of observations is n-1
  return( list(
    "param"=list("B"=minSSE$par[1],
                 "f"=minSSE$par[2],
                 "mutbias"=exp(minSSE$par[3]),
                 "e1"=minSSE$par[4],
                 "e2"=minSSE$par[5]),
    "criteria"=list("SStot"=SStot,
                    "SSres"=SSres,
                    "R2"=R2,
                    "AIC"=AIC) )
  )
}



#' @title Sum of squares minimization of model Hotpsot 2 with error
#'
#' @description Function that searches for the six parameters that minimize the sum_of_squares function
#'
#' @param WS the WS observed SFS
#' @param SW the SW observed SFS
#' @param GC GC content
#' @param cor a Boolean to add a correction to the transformed SFS (default = TRUE)
#' @param snpthresh the value below which the SNP category is removed (default = 1)
#' @param B0min minimum for the range of B0, default value = -100
#' @param B0max maximum for the range of B0, default value = 100
#' @param B1min minimum for the range of B1, default value = -100
#' @param B1max maximum for the range of B1, default value = 100
#' @param Mmin minimum for the range of M, default value = -5
#' @param Mmax maximum or the range of M, default value = 5
#' @param e1max maximum for the range of WS error rate, default value = 0.49
#' @param e2max maximum for the range of SW error rate, default value = 0.49
#' @param e1init initial starting value for optimization for WS error rates, default = 0.01
#' @param e2init initial starting value for optimization for SW error rates, default = 0.01
#' @param Maxit maximum number of iterations (option for optim), see manual, default value = 100
#' @param Factr level of control the convergence (option for optim, see manual), default value = 10^7
#' @param Lmm number of updates in the method (option for optim, see manual)
#' @param Verbose from 0 (default) to 5: level of outputs during optimization (option for optim, see manual)
#' @param Usegr to use (default) or not the analytical gradient (option for optim, see manual)
#' It's better to use the analytical gradient function but the option can be turn off for testing
#'
#' @returns A list of two lists:
#' - param: Optimized parameters
#' - criteria: Model fit criteria
#'
#' @export
#'
#' @importFrom stats lm
#' @importFrom stats optim
#'
#' @examples
#' sfsWS <- c(1000,500,300,200,80,50)
#' sfsSW <- c(2000,800,400,150,50,10)
#' LS <- least_square_hotspot2_err(sfsWS,sfsSW,0.5)
#' LS$B0 # population-scaled background gBGC
#' LS$B1 # population-scaled hotspot gBGC
#' LS$f # proportion of hotspots
least_square_hotspot2_err <- function(WS,SW,GC,cor=COR,snpthresh=SNPTHRESH,
                                  B0min=BMIN,B0max=BMAX,B1min=BMIN,B1max=BMAX,Mmin=MMIN,Mmax=MMAX,e1max=EMAX,e2max=EMAX,e1init=EINIT,e2init=EINIT,
                                  Maxit=MAXIT,Factr=FACTR,Lmm=LMM,Verbose=VERBOSE,Usegr=USEGR) {
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
  x <- c(1:n)/(n+1)
  NONZERO <- which(WS!=0 & SW!=0)
  # Arbitrary starting value fo f. Note that we assume 0 ≤ f ≤ 1/2
  finit <- 0.1
  # Determination of initial values for optimization: use of the simple regression
  reginit <- lm(log(WS[NONZERO]/SW[NONZERO]) ~ x[NONZERO])
  Minit <- -reginit$coef[1] + log(1 - GC) - log(GC)
  # Starting values such that B1 = 5*B0
  B0init <- reginit$coef[2]/(1 + 4*finit)
  B1init <- 5*reginit$coef[2]/(1 + 4*finit)
  init <- c(B0init,B1init,finit,Minit,e1init,e2init)
  # Boundaries for optimization
  inf <- c(B0min,B1min,0,Mmin,0,0)
  sup <- c(B0max,B1max,1/2,Mmax,e1max,e2max)
  if(Usegr) gradient <- gr_sum_of_squares_hotspot2_err else gradient <- NULL
  SCALE <- abs(init)
  minSSE <- optim(
    par = init,
    fn = sum_of_squares_hotspot2_err,
    gr = gradient,
    WS = WS, SW = SW, GC = GC, cor = cor, snpthresh = snpthresh,
    lower = inf,upper = sup,
    method = "L-BFGS-B",
    control=list(parscale=SCALE,maxit=Maxit,factr=Factr,lmm=Lmm,trace=Verbose))
  SStot <- sum_of_squares_NULL(SW,WS)$SStot
  SSres <- minSSE$value
  R2 <- 1 - SSres/SStot
  AIC <- AICls(length(which(WS!=0 & SW!=0)),length(init),SSres)
  # Note that the number of observations is n-1
  return( list(
    "param"=list("B0"=minSSE$par[1],
                 "B1"=minSSE$par[2],
                 "f"=minSSE$par[3],
                 "mutbias"=exp(minSSE$par[4]),
                 "e1"=minSSE$par[5],
                 "e2"=minSSE$par[6]),
    "criteria"=list("SStot"=SStot,
                    "SSres"=SSres,
                    "R2"=R2,
                    "AIC"=AIC) )
  )
}



#' @title Sum of squares minimization of model Hotpsot 2bis with error
#'
#' @description Function that searches for the five parameters that minimize the sum_of_squares function
#' This model is the same as model 2 except that f is fixed by the user
#'
#' @param WS the WS observed SFS
#' @param SW the SW observed SFS
#' @param GC GC content
#' @param f proportion of hotspots (fixed by the user: 0 ≤ f ≤ 1/2)
#' @param cor a Boolean to add a correction to the transformed SFS (default = TRUE)
#' @param snpthresh the value below which the SNP category is removed (default = 1)
#' @param B0min minimum for the range of B0, default value = -100
#' @param B0max maximum for the range of B0, default value = 100
#' @param B1min minimum for the range of B1, default value = -100
#' @param B1max maximum for the range of B1, default value = 100
#' @param Mmin minimum for the range of M, default value = -5
#' @param Mmax maximum or the range of M, default value = 5
#' @param e1max maximum for the range of WS error rate, default value = 0.49
#' @param e2max maximum for the range of SW error rate, default value = 0.49
#' @param e1init initial starting value for optimization for WS error rates, default = 0.01
#' @param e2init initial starting value for optimization for SW error rates, default = 0.01
#' @param Maxit maximum number of iterations (option for optim), see manual, default value = 100
#' @param Factr level of control the convergence (option for optim, see manual), default value = 10^7
#' @param Lmm number of updates in the method (option for optim, see manual)
#' @param Verbose from 0 (default) to 5: level of outputs during optimization (option for optim, see manual)
#' @param Usegr to use (default) or not the analytical gradient (option for optim, see manual)
#' It's better to use the analytical gradient function but the option can be turn off for testing
#'
#' @returns A list of two lists:
#' - param: Optimized parameters
#' - criteria: Model fit criteria
#'
#' @importFrom stats lm
#' @importFrom stats optim
#'
#' @export
#'
#' @examples
#' sfsWS <- c(1000,500,300,200,80,50)
#' sfsSW <- c(2000,800,400,150,50,10)
#' LS <- least_square_hotspot2bis_err(sfsWS,sfsSW,0.5,0.2)
#' LS$B0 # population-scaled background gBGC
#' LS$B1 # population-scaled hotspot gBGC
least_square_hotspot2bis_err <- function(WS,SW,GC,f,cor=COR,snpthresh=SNPTHRESH,
                                     B0min=BMIN,B0max=BMAX,B1min=BMIN,B1max=BMAX,Mmin=MMIN,Mmax=MMAX,e1max=EMAX,e2max=EMAX,e1init=EINIT,e2init=EINIT,
                                     Maxit=MAXIT,Factr=FACTR,Lmm=LMM,Verbose=VERBOSE,Usegr=USEGR) {
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
  if(f<0 | f>1/2) {
    abort("f must be between 0 and 1/2")
  }
  #Main
  n <- length(WS)
  x <- c(1:n)/(n+1)
  NONZERO <- which(WS!=0 & SW!=0)
  # Determination of initial values for optimization: use of the simple regression
  reginit <- lm(log(WS[NONZERO]/SW[NONZERO]) ~ x[NONZERO])
  Minit <- -reginit$coef[1] + log(1 - GC) - log(GC)
  # Starting values such that B1 = 5*B0
  B0init <- reginit$coef[2]/(1 + 4*f)
  B1init <- 5*reginit$coef[2]/(1 + 4*f)
  init <- c(B0init,B1init,Minit,e1init,e2init)
  # Boundaries for optimization
  inf <- c(B0min,B1min,Mmin,0,0)
  sup <- c(B0max,B1max,Mmax,e1max,e2max)
  if(Usegr) gradient <- gr_sum_of_squares_hotspot2bis_err else gradient <- NULL
  SCALE <- abs(init)
  minSSE <- optim(
    par = init,
    fn = sum_of_squares_hotspot2bis_err,
    gr = gradient,
    WS = WS, SW = SW, GC = GC, f = f, cor = cor, snpthresh = snpthresh,
    lower = inf,upper = sup,
    method = "L-BFGS-B",
    control=list(parscale=SCALE,maxit=Maxit,factr=Factr,lmm=Lmm,trace=Verbose))
  SStot <- sum_of_squares_NULL(SW,WS)$SStot
  SSres <- minSSE$value
  R2 <- 1 - SSres/SStot
  AIC <- AICls(length(which(WS!=0 & SW!=0)),length(init),SSres)
  # Note that the number of observations is n-1
  return( list(
    "param"=list("B0"=minSSE$par[1],
                 "B1"=minSSE$par[2],
                 "mutbias"=exp(minSSE$par[3]),
                 "e1"=minSSE$par[4],
                 "e2"=minSSE$par[5]),
    "criteria"=list("SStot"=SStot,
                    "SSres"=SSres,
                    "R2"=R2,
                    "AIC"=AIC) )
  )
}
