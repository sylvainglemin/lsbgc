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
  #log(sfs) + 0.444/sfs - 0.107/(sfs^2) - 0.927/(sfs^3)
  log(sfs) + 0.987/sfs - 1.557/(sfs^2)
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
  return(0.5/sfs - 0.69/(sfs)^2)
  # This is the full version of the variance
  # However exponentials lead to numerical issues
  # exp_sfs <- exp(sfs)
  # exp_neg_sfs <- exp(-sfs)
  # denominator <- 4 * (-1 + exp_neg_sfs)^2 * (-1 + exp_sfs) * sfs^8
  # term_exp <- exp(-10 * sfs)
  # term_cubed <- (-1 + exp_sfs)^3
  # term_linear <- (-1 + exp_sfs - sfs)
  # inner_A1 <- 33.48 + 0.66 * sfs
  # inner_A2 <- -33.48 + (-1.32 + 0.88 * sfs) * sfs
  # inner_A3 <- 11.16 + sfs * (0.66 + sfs * (-0.88 + 1 * sfs))
  # term_A <- -11.16 + exp_sfs * (inner_A1 + exp_sfs * (inner_A2 + exp_sfs * inner_A3))
  # inner_B1 <- 8.37 + 0.22 * sfs
  # inner_B2 <- -8.37 + (-0.44 + 0.44 * sfs) * sfs
  # inner_B3 <- 2.79 + sfs * (0.22 + sfs * (-0.44 + 1 * sfs))
  # term_B <- -2.79 + exp_sfs * (inner_B1 + exp_sfs * (inner_B2 + exp_sfs * inner_B3))
  # combined_terms <- (-1) * (-1 + exp_sfs -  sfs) * (term_A)^2 +
  #   2 * exp_sfs * sfs * (term_B)^2
  # t_var <- (1 / denominator) * term_exp * term_cubed * term_linear * combined_terms
  # return(t_var)
}


