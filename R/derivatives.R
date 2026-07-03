# Estimation of selection and GC-biased gene conversion intensity from site frequency spectrum data by a least-square approach
# Sylvain Glemin
# sylvain.glemin@univ-rennes.fr


# Derivatives of several sub-functions used in the computation of gradients


#' @title derivatives of the expectations of the log of the ratio the WS and SW SFS
#'
#' @param WS W->S SFS
#' @param SW S->W SFS
#' @param e1 W->S polarization error
#' @param e2 S->W polarization error
#'
#' @returns a list with the derivatives of the two ratios
#'
#' @noRd
d_expected_log_ratio <- function(WS,SW,e1,e2){
  d1 <- (SW + rev(WS))/((1 - e1)*SW - e1*rev(WS)) - # Log term
    1/((1 - e2)*WS - e2*rev(SW)) + # First ratio term
    ((1 - e2)*rev(WS) - e2*SW)/(((1 - e1)*SW - e1*rev(WS))^2) # Second ratio term
  d2 <- (rev(SW) + WS)/(e2*rev(SW) - (1 - e2)*WS) + # Log term
    ((1 - e1)*rev(SW) - e1*WS)/(((1 - e2)*WS - e2*rev(SW))^2) - # First ratio term
    1/((e1*rev(WS) - (1 - e1)*SW))
  return(list("d1"=d1,"d2"=d2))
}

#' @title ratio of expectation of W->S and S->W SFS
#'
#' @description
#' This function is an intermediate function for some derivative functions
#'
#' @param B population scale BGC coefficent = 4Neb
#' @param x allele frequency of the S allele
#'
#' @returns ratio of expectation of W->S and S->W SFS
#'
#' @noRd
ratio_B <- function(B,x){
  if(B==0) return(1 - x)
  return(
    (exp(B) - exp(B*x))/((exp(B) - 1))
  )
}
ratio_B <- Vectorize(ratio_B)

#' @title derivative of the ratio of expectation of W->S and S->W SFS for constant B
#'
#' @param B population scale BGC coefficent = 4Neb
#' @param x allele frequency of the S allele
#'
#' @returns derivative of the ratio of expectation of W->S and S->W SFS
#'
#' @noRd
d_ratio_B <- function(B,x){
  if(B==0) return(
    x*(1 - x)/2
  )
  return(
    (-exp(B) + exp(B*x)*(exp(B)* (1 - x) + x))/(-1 + exp(B))^2
  )
}
d_ratio_B <- Vectorize(d_ratio_B)

#' @title derivative of the ratio of expectation of W->S and S->W SFS for hotspot model 1
#'
#' @param B population scale BGC coefficent in hotspots = 4Neb
#' @param f fraction of hotpots
#' @param x allele frequency of the S allele
#'
#' @returns derivative of the ratio of expectation of W->S and S->W SFS
#'
#' @noRd
d_hotspot1 <- function(B,f,x) {
  dB <- -f*(
    d_ratio_B(B,x)/((1-f)*(1-x) + f*ratio_B(B,x)) +
      d_ratio_B(-B,x)/((1-f)*(1-x) + f*ratio_B(-B,x))
  )
  df <- ((1-x)*ratio_B(-B,x)-(1-x)*ratio_B(B,x))/
    ( ((1-f)*(1-x) + f*ratio_B(B,x))*
        ((1-f)*(1-x) + f*ratio_B(-B,x)) )
  return(list("dB"=dB,"df"=df))
}

#' @title derivative of the ratio of expectation of W->S and S->W SFS for hotspot model 1
#'
#' @param B0 population scale BGC coefficent in the background = 4Neb0
#' @param B1 population scale BGC coefficent in hotspots = 4Neb1
#' @param f fraction of hotpots
#' @param x allele frequency of the S allele
#'
#' @returns derivative of the ratio of expectation of W->S and S->W SFS
#'
#' @noRd
d_hotspot2 <- function(B0,B1,f,x) {
  dB0 <- -(1-f)*(
    d_ratio_B(B0,x)/((1-f)*ratio_B(B0,x) + f*ratio_B(B1,x)) +
      d_ratio_B(-B0,x)/((1-f)*ratio_B(-B0,x) + f*ratio_B(-B1,x))
  )
  dB1 <- -f*(
    d_ratio_B(B1,x)/((1-f)*ratio_B(B0,x) + f*ratio_B(B1,x)) +
      d_ratio_B(-B1,x)/((1-f)*ratio_B(-B0,x) + f*ratio_B(-B1,x))
  )
  df <- (ratio_B(B0,x)*ratio_B(-B1,x)-ratio_B(-B0,x)*ratio_B(B1,x))/
    ( ((1-f)*ratio_B(B0,x) + f*ratio_B(B1,x))*
        ((1-f)*ratio_B(-B0,x) + f*ratio_B(-B1,x)) )
  return(list("dB0"=dB0,"dB1"=dB1,"df"=df))
}
