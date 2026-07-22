# Set of function to simulate sfs

#' @title Expected SFS constant gBGC
#'
#' @description
#' Generate the expected WS and SW SFS.
#' The default is a neutral SFS without bias so that the two generated SFS are equal.
#'
#' @param n sample size
#' @param theta population scale mutation rate for WS mutations = 4NeuWS
#' @param B population scale gBGC coefficient = 4Neb
#' @param mut_bias S to W mutation bias so uSW = mut_bias*uWS
#' @param GC GC content
#' @param eWS polarization error for WS mutations (default = 0)
#' @param eSW polarization error for SW mutations (default = 0)
#' @param vect_r vector of noise parameters (default = rep(1,(n-1)))
#'
#' @returns the expected WS and SW SFS
#'
#' @export
#'
#' @examples
#' sfs <- expected_sfs_constant(20,100,1,2,0.5,0,0)
#' barplot(sfs$WS)
#' barplot(sfs$SW)
expected_sfs_constant <- function(n,theta,B=0,mut_bias=1,GC=0.5,eWS=0,eSW=0,vect_r=rep(1,(n-1))){
  #Error checking
  if(!is.numeric(n) || floor(n)!=n || n < 2){
    abort("n must be a positive integer higher than 1")
  }
  if(!is.numeric(theta) ||  theta < 0){
    abort("theta must be a positive numerical value")
  }
  if(!is.numeric(B) ){
    abort("B must be a numerical value")
  }
  if(!is.numeric(mut_bias) ||  mut_bias < 0){
    abort("mut_bias must be a positive numerical value")
  }
  if(!is.numeric(GC) ||  GC <= 0 || GC >= 1){
    abort("GC must be a numerical value strictly higher than 0 and lower than 1")
  }
  if(!is.numeric(eWS) ||  eWS < 0 || eWS > 1){
    abort("eWS must be a numerical value higher than 0 and lower than 1")
  }
  if(!is.numeric(eSW) ||  eSW < 0 || eSW > 1){
    abort("eSW must be a numerical value higher than 0 and lower than 1")
  }
  if (!is.numeric(vect_r) || length(vect_r) != (n - 1) || any(vect_r < 0)) {
    abort("vect_r must be a vector of positive numerical values of length n-1")
  }
  #Main
  thetaWS <- theta*(1-GC)
  thetaSW <- theta*mut_bias*GC
  sfsWS <- vect_r*ens_constant_err(thetaWS,thetaSW,B,eWS,eSW,n,j = 1:(n-1))
  sfsSW <- vect_r*ens_constant_err(thetaSW,thetaWS,-B,eSW,eWS,n,j = 1:(n-1))
  return(list("WS"=sfsWS,"SW"=sfsSW))
}



#' @title Expected SFS hotspot gBGC
#'
#' @description
#' Generate the expected WS and SW SFS for the hotspot model
#' If B0 = B1 = 0 it corresponds to the neutral model.
#' Yet, if you want to generate a neutral sfs it is simpler to use the "expected_sfs_constant" function
#'
#' @param n sample size
#' @param theta population scale mutation rate for WS mutations = 4NeuWS
#' @param B0 population scale gBGC coefficient for the background = 4Neb0
#' @param B1 population scale gBGC coefficient for hotspots = 4Neb1
#' @param f fraction of hotspots
#' @param mut_bias S to W mutation bias so uSW = mut_bias*uWS
#' @param GC GC content
#' @param eWS polarization error for WS mutations (default = 0)
#' @param eSW polarization error for SW mutations (default = 0)
#' @param vect_r vector of noise parameters (default = rep(1,(n-1)))
#'
#' @returns expected WS and SW SFS
#'
#' @export
#'
#' @examples
#' sfs <- expected_sfs_hotspot(20,100,3,0.1,0.1,2,0.5,0,0)
#' barplot(sfs$WS)
#' barplot(sfs$SW)
expected_sfs_hotspot <- function(n,theta,B0,B1,f,mut_bias,GC,eWS=0,eSW=0,vect_r=rep(1,(n-1))){
  #Error checking
  if(!is.numeric(n) || floor(n)!=n || n < 2){
    abort("n must be a positive integer higher than 1")
  }
  if(!is.numeric(theta) ||  theta < 0){
    abort("theta must be a positive numerical value")
  }
  if(!is.numeric(B0) ){
    abort("B0 must be a numerical value")
  }
  if(!is.numeric(B1) ){
    abort("B1 must be a numerical value")
  }
  if(!is.numeric(f) ||  f < 0 || f > 1){
    abort("f must be a numerical value higher than 0 and lower than 1")
  }
  if(!is.numeric(mut_bias) ||  mut_bias < 0){
    abort("mut_bias must be a positive numerical value")
  }
  if(!is.numeric(GC) ||  GC <= 0 || GC >=1){
    abort("GC must be a numerical value strictly higher than 0 and lower than 1")
  }
  if(!is.numeric(eWS) ||  eWS < 0 || eWS > 1){
    abort("eWS must be a numerical value higher than 0 and lower than 1")
  }
  if(!is.numeric(eSW) ||  eSW < 0 || eSW > 1){
    abort("eSW must be a numerical value higher than 0 and lower than 1")
  }
  if (!is.numeric(vect_r) || length(vect_r) != (n - 1) || any(vect_r < 0)) {
    abort("vect_r must be a vector of positive numerical values of length n-1")
  }
  #Main
  thetaWS <- theta*(1-GC)
  thetaSW <- theta*mut_bias*GC
  sfsWS <- vect_r*ens_hotspot2_err(thetaWS,thetaSW,B0,B1,f,eWS,eSW,n,j = 1:(n-1))
  sfsSW <- vect_r*ens_hotspot2_err(thetaSW,thetaWS,-B0,-B1,f,eSW,eWS,n,j = 1:(n-1))
  return(list("WS"=sfsWS,"SW"=sfsSW))
}


#' @title Simulated SFS constant gBGC
#'
#' @description
#' Generate a random WS and SW SFS.
#' The default is a neutral SFS without bias so that the two generated SFS are equal.
#'
#' @param n sample size
#' @param theta population scale mutation rate for WS mutations = 4NeuWS
#' @param B population scale gBGC coefficient = 4Neb
#' @param mut_bias S to W mutation bias so uSW = mut_bias*uWS
#' @param GC GC content
#' @param eWS polarization error for WS mutations
#' @param eSW polarization error for SW mutations
#' @param vect_r vector of noise parameters (default = rep(1,(n-1)))
#' @param fix_snp choice between three possibilities:
#' "free" (default): in each frequency class a Poisson number of SNPs is drawn so that the total number is also a random variable.
#' "fixed_sfs": for each SFS (WS and SW) the number of sampled SNPs is fixed to the expected number of the model.
#' Then a multinomial sampling is drawn for each SFS
#' "fixed_tot": the total number of sampled SNPs is fixed (WS + SW)
#' Then a multinomial sampling is drawn fro the combined SFS
#'
#' @returns the expected SFS
#'
#' @importFrom stats rpois rmultinom
#'
#' @export
#'
#' @examples
#' sim_sfs <- simulated_sfs_constant(20,100,1,2,0.5,0,0)
#' barplot(sim_sfs$WS)
#' barplot(sim_sfs$SW)
simulated_sfs_constant <- function(n,theta,B=0,mut_bias=1,GC=0.5,eWS=0,eSW=0,vect_r=rep(1,(n-1)),fix_snp="free"){
  #Error checking
  if(!is.numeric(n) || floor(n)!=n || n < 2){
    abort("n must be a positive integer higher than 1")
  }
  if(!is.numeric(theta) ||  theta < 0){
    abort("theta must be a positive numerical value")
  }
  if(!is.numeric(B) ){
    abort("B must be a numerical value")
  }
  if(!is.numeric(mut_bias) ||  mut_bias < 0){
    abort("mut_bias must be a positive numerical value")
  }
  if(!is.numeric(GC) ||  GC <= 0 || GC >=1){
    abort("GC must be a numerical value strictly higher than 0 and lower than 1")
  }
  if(!is.numeric(eWS) ||  eWS < 0 || eWS > 1){
    abort("eWS must be a numerical value higher than 0 and lower than 1")
  }
  if(!is.numeric(eSW) ||  eSW < 0 || eSW > 1){
    abort("eSW must be a numerical value higher than 0 and lower than 1")
  }
  if (!is.numeric(vect_r) || length(vect_r) != (n - 1) || any(vect_r < 0)) {
    abort("vect_r must be a vector of positive numerical values of length n-1")
  }
  #Main
  exp_sfs <- expected_sfs_constant(n, theta,B,mut_bias,GC,eWS,eSW,vect_r)
  if(fix_snp=="free") {
    sfsWS <- sapply(exp_sfs$WS,function(x) rpois(1,x))
    sfsSW <- sapply(exp_sfs$SW,function(x) rpois(1,x))
    return(list("WS"=sfsWS,"SW"=sfsSW))
  }
  if(fix_snp=="fixed_sfs") {
    nWS <- sum(exp_sfs$WS)
    probaWS <- exp_sfs$WS/nWS
    sfsWS <- rmultinom(1,nWS,probaWS)
    nSW <- sum(exp_sfs$SW)
    probaSW <- exp_sfs$SW/nSW
    sfsSW <- rmultinom(1,nSW,probaSW)
    return(list("WS"=sfsWS,"SW"=sfsSW))
  }
  if(fix_snp=="fixed_tot") {
    ntot <- sum(exp_sfs$WS+exp_sfs$SW)
    probatot <- c(exp_sfs$WS,exp_sfs$SW)/ntot
    sfstot <- rmultinom(1,ntot,probatot)
    sfsWS <- sfstot[1:(n-1)]
    sfsSW <- sfstot[n:(n-2)]
    return(list("WS"=sfsWS,"SW"=sfsSW))
  } else {
    abort("The chosen option is not available: must be free, fixed_sfs or fixed_tot")
  }
}


#' @title Simulated SFS hotspot model
#'
#' @description
#' Generate a random WS and SW SFS.
#' If B0 = B1 = 0 it corresponds to the neutral model.
#' Yet, if you want to generate a neutral sfs it is simpler to use the "simulated_sfs_constant" function
#'
#' @param n sample size
#' @param theta population scale mutation rate for WS mutations = 4NeuWS
#' @param B0 population scale gBGC coefficient for the background = 4Neb0
#' @param B1 population scale gBGC coefficient for hotspots = 4Neb1
#' @param f fraction of hotspots
#' @param mut_bias S to W mutation bias so uSW = mut_bias*uWS
#' @param GC GC content
#' @param eWS polarization error for WS mutations
#' @param eSW polarization error for SW mutations
#' @param vect_r vector of noise parameters (default = rep(1,(n-1)))
#' @param fix_snp choice between three possibilities:
#' "free" (default): in each frequency class a Poisson number of SNPs is drawn so that the total number is also a random variable.
#' "fixed_sfs": for each SFS (WS and SW) the number of sampled SNPs is fixed to the expected number of the model.
#' Then a multinomial sampling is drawn for each SFS
#' "fixed_tot": the total number of sampled SNPs is fixed (WS + SW)
#' Then a multinomial sampling is drawn fro the combined SFS
#'
#' @returns the expected SFS
#'
#' @importFrom stats rpois rmultinom
#'
#' @export
#'
#' @examples
#' sim_sfs <- simulated_sfs_constant(20,100,1,2,0.5,0,0)
#' barplot(sim_sfs$WS)
#' barplot(sim_sfs$SW)
simulated_sfs_hotspot <- function(n,theta,B0,B1,f,mut_bias,GC,eWS=0,eSW=0,vect_r=rep(1,(n-1)),fix_snp="free"){
  #Error checking
  if(!is.numeric(n) || floor(n)!=n || n < 2){
    abort("n must be a positive integer higher than 1")
  }
  if(!is.numeric(theta) ||  theta < 0){
    abort("theta must be a positive numerical value")
  }
  if(!is.numeric(B0) ){
    abort("B0 must be a numerical value")
  }
  if(!is.numeric(B1) ){
    abort("B1 must be a numerical value")
  }
  if(!is.numeric(f) ||  f < 0 || f > 1){
    abort("f must be a numerical value higher than 0 and lower than 1")
  }
  if(!is.numeric(mut_bias) ||  mut_bias < 0){
    abort("mut_bias must be a positive numerical value")
  }
  if(!is.numeric(GC) ||  GC <= 0 || GC >=1){
    abort("GC must be a numerical value strictly higher than 0 and lower than 1")
  }
  if(!is.numeric(eWS) ||  eWS < 0 || eWS > 1){
    abort("eWS must be a numerical value higher than 0 and lower than 1")
  }
  if(!is.numeric(eSW) ||  eSW < 0 || eSW > 1){
    abort("eSW must be a numerical value higher than 0 and lower than 1")
  }
  if (!is.numeric(vect_r) || length(vect_r) != (n - 1) || any(vect_r < 0)) {
    abort("vect_r must be a vector of positive numerical values of length n-1")
  }
  #Main
  exp_sfs <- expected_sfs_hotspot(n,theta,B0,B1,f,mut_bias,GC,eWS,eSW,vect_r)
  if(fix_snp=="free") {
    sfsWS <- sapply(exp_sfs$WS,function(x) rpois(1,x))
    sfsSW <- sapply(exp_sfs$SW,function(x) rpois(1,x))
    return(list("WS"=sfsWS,"SW"=sfsSW))
  }
  if(fix_snp=="fixed_sfs") {
    nWS <- sum(exp_sfs$WS)
    probaWS <- exp_sfs$WS/nWS
    sfsWS <- rmultinom(1,nWS,probaWS)
    nSW <- sum(exp_sfs$SW)
    probaSW <- exp_sfs$SW/nSW
    sfsSW <- rmultinom(1,nSW,probaSW)
    return(list("WS"=sfsWS,"SW"=sfsSW))
  }
  if(fix_snp=="fixed_sfs") {
    ntot <- sum(exp_sfs$WS+exp_sfs$SW)
    probatot <- c(exp_sfs$WS,exp_sfs$SW)/ntot
    sfstot <- rmultinom(1,ntot,probatot)
    sfsWS <- sfstot[1:(n-1)]
    sfsSW <- sfstot[n:(n-2)]
    return(list("WS"=sfsWS,"SW"=sfsSW))
  } else {
    abort("The chosen option is not available: must be free, fixed_sfs or fixed_tot")
  }
}
