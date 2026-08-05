# Estimation of selection and GC-biased gene conversion intensity from site frequency spectrum data by a least-square approach
# Sylvain Glemin
# sylvain.glemin@univ-rennes.fr

# Functions to transform initial data


#' @title Log transformation of the SFS plus additional terms to remove bias
#' (see article for mathematical justification)
#'
#' @param sfs a site frequency spectrum (a numerical vector)
#' @param cor a Boolean to add a correction to the transformed SFS (default = FALSE)
#'
#' @returns the transformed sfs
#' @export
#'
#' @examples
#' # A WS and SW SFS
#' sfs_ws <- c(100,50,30,20,3)
#' sfs_sw <- c(300,140,80,50,7)
#' # Simple log ratio
#' l_ratio <- log(sfs_ws) - log(sfs_sw)
#' # Transformed ratio with bias correction
#' t_ratio <- t_sfs(sfs_ws) - t_sfs(sfs_sw)
#' # Comparison
#' plot(l_ratio)
#' points(t_ratio,col="red")
t_sfs <- function(sfs,cor=COR) {
  #Error checking
  if(!is.numeric(sfs)){
    abort("sfs must be a numeric vector")
  }
  # Main
  if(COR) {
    return( log(sfs) + 1/(2*sfs) - 1/(2*sfs^2) - 0.04*exp(-sfs) - 1.19*exp(-sfs^2))
    #return( log(sfs) + 1/(2*sfs) - 1/(2*sfs^2))
    # return( log(sfs) + 1/(2*sfs) - 1/(2*sfs^2) - 1/(3*sfs^3) - 1/(4*sfs^4) )
    #return(log(sfs) -0.05/sfs + 4.6/(sfs^2) -5.4/(sfs^3))
  }
  else {
    return(log(sfs) + 1/(2*sfs))
  }

}

#' @title Expected variance of the transformed sfs
#' Assuming that each class follows a Poisson dsitribution
#'
#' @param sfs a site frequency spectrum (a numerical vector)
#' @param cor a Boolean to add a correction to the transformed SFS (default = FALSE)
#'
#' @returns the expected variance of the sfs
#' @export
#'
#' @examples
#' #' A WS and SW SFS
#' sfs_ws <- c(100,50,30,20,3)
#' sfs_sw <- c(300,140,80,50,7)
#' # Transformed ratio with bias correction
#' t_ratio <- t_sfs(sfs_ws) - t_sfs(sfs_sw)
#' # The weight for each class is given by the inverse of the variance
#' weight <- 1/(t_variance(sfs_ws) + t_variance(sfs_sw))
t_variance <- function(sfs,cor=COR) {
  #Error checking
  if(!is.numeric(sfs)){
    abort("sfs must be a numeric vector")
  }
  # Main
  if(cor) {
    #a0 <- (1 - 2*sfs)^2/(4*sfs^3)
    #a1 <- -(3 - 7*sfs + 4*sfs^3)/(4*sfs^3)
    #a2 <- (3 -2*sfs - 4*sfs^2)/(4*sfs^3)
    #a3 <- -(1 + sfs)/(4*sfs^3)
    #return( a0 + a1*exp(-sfs) + a2*exp(-2*sfs) + a3*exp(-3*sfs) )
    return(1/sfs - 1/(sfs^2) + 1/(4*sfs^3))
    #return(1/sfs - 1/(sfs^2))
    # return(1/sfs)
    # var <- 1/sfs + 0.1/(sfs)^2
    # corvar <- c(0, 0.54, 0.43, 0.24, 0.11, 0.04,rep(0,max(sfs))) # Manual correcting factor
    # # Interpolation for non-integer values
    # corvarextrapol <- function(x) {
    #   xplus <- ceiling(x)
    #   xminus <- max(1,floor(x))
    #   f <- xplus - x
    #   return(f*corvar[xminus]+(1-f)*corvar[xplus])
    # }
    # var <- var+corvarextrapol(sfs)
    return(var)
  }
  else {
    return(1/sfs)
  }

}


