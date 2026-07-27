# Estimation of selection and GC-biased gene conversion intensity from site frequency spectrum data by a least-square approach
# Sylvain Glemin
# sylvain.glemin@univ-rennes.fr

# Functions to transform initial data


#' @title Log transformation of the SFS plus additional terms to remove bias
#' (see article for mathematical justification)
#'
#' @param sfs a site frequency spectrum (a numerical vector)
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
t_sfs <- function(sfs) {
  #Error checking
  if(!is.numeric(sfs)){
    abort("sfs must be a numeric vector")
  }
  # Main
  #log(sfs) -0.0654/sfs + 4.716/(sfs^2) -5.526/(sfs^3)
  log(sfs) -0.05/sfs + 4.6/(sfs^2) -5.4/(sfs^3)

}

#' @title Expected variance of the transformed sfs
#' Assuming that each class follows a Poisson dsitribution
#'
#' @param sfs a site frequency spectrum (a numerical vector)
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
t_variance <- function(sfs) {
  #Error checking
  if(!is.numeric(sfs)){
    abort("sfs must be a numeric vector")
  }
  # Main
  var <- 1/sfs + 0.1/(sfs)^2
  corvar <- c(0, 0.54, 0.43, 0.24, 0.11, 0.04,0) # Manual correcting factor computed for integers between 2 and 6
  # Interpolation for non-integer values
  corvarextrapol <- function(x) {
    xplus <- ceiling(x)
    xminus <- floor(x)
    f <- xplus - x
    return(f*corvar[xminus]+(1-f)*corvar[xplus])
  }
  var <- ifelse(sfs>6,var,
                ifelse(sfs<2,1.1, # 1 and 2 have roughly the same value: plateau between them
                       var+corvarextrapol(sfs) # linear interpolation of the correcting factor
                       )
                )
  return(var)
}


