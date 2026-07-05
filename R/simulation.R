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
#'
#' @returns the expected WS and SW SFS
#'
#' @export
#'
#' @examples
#' sfs <- expected_sfs_constant(20,100,1,2,0.5,0,0)
#' barplot(sfs$WS)
#' barplot(sfs$SW)
expected_sfs_constant <- function(n,theta,B=0,mut_bias=1,GC=0.5,eWS=0,eSW=0){
  thetaWS <- theta*(1-GC)
  thetaSW <- theta*mut_bias*GC
  freq <- c(1:(n-1))
  sfsWS <- sapply(freq, function(j) ens_constant_err(thetaWS,thetaSW,B,eWS,eSW,n,j))
  sfsSW <- sapply(freq, function(j) ens_constant_err(thetaSW,thetaWS,-B,eSW,eWS,n,j))
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
#'
#' @returns expected WS and SW SFS
#'
#' @export
#'
#' @examples
#' sfs <- expected_sfs_hotspot(20,100,3,0.1,0.1,2,0.5,0,0)
#' barplot(sfs$WS)
#' barplot(sfs$SW)
expected_sfs_hotspot <- function(n,theta,B0,B1,f,mut_bias,GC,eWS=0,eSW=0){
  thetaWS <- theta*(1-GC)
  thetaSW <- theta*mut_bias*GC
  freq <- c(1:(n-1))
  sfsWS <- sapply(freq, function(j) ens_hotspot2_err(thetaWS,thetaSW,B0,B1,f,eWS,eSW,n,j))
  sfsSW <- sapply(freq, function(j) ens_hotspot2_err(thetaSW,thetaWS,-B0,-B1,f,eSW,eWS,n,j))
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
#' @param fix_snp choice between three possibilities:
#' "free" (default): in each frequency class a Poisson number of SNPs is drawn so that the total number is also a random variable.
#' "fixed_sfs": for each SFS (WS and SW) the number of sampled SNPs is fixed to the expected number of the model.
#' Then a multinomial sampling is drawn for each SFS
#' "fixed_tot": the total number of sampled SNPs is fixed (WS + SW)
#' Then a multinomial sampling is drawn fro the combined SFS
#'
#' @returns the expected SFS
#'
#' @export
#'
#' @examples
#' sim_sfs <- simulated_sfs_constant(20,100,1,2,0.5,0,0)
#' barplot(sim_sfs$WS)
#' barplot(sim_sfs$SW)
simulated_sfs_constant <- function(n,theta,B=0,mut_bias=1,GC=0.5,eWS=0,eSW=0,fix_snp="free"){
  exp_sfs <- expected_sfs_constant(n, theta,B,mut_bias,GC,eWS,eSW)
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
    print("ERROR: the chosen option is not available")
    return(NA)
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
#' @param fix_snp choice between three possibilities:
#' "free" (default): in each frequency class a Poisson number of SNPs is drawn so that the total number is also a random variable.
#' "fixed_sfs": for each SFS (WS and SW) the number of sampled SNPs is fixed to the expected number of the model.
#' Then a multinomial sampling is drawn for each SFS
#' "fixed_tot": the total number of sampled SNPs is fixed (WS + SW)
#' Then a multinomial sampling is drawn fro the combined SFS
#'
#' @returns the expected SFS
#'
#' @export
#'
#' @examples
#' sim_sfs <- simulated_sfs_constant(20,100,1,2,0.5,0,0)
#' barplot(sim_sfs$WS)
#' barplot(sim_sfs$SW)
simulated_sfs_hotspot <- function(n,theta,B0,B1,f,mut_bias,GC,eWS=0,eSW=0,fix_snp="free"){
  exp_sfs <- expected_sfs_hotspot(n,theta,B0,B1,f,mut_bias,GC,eWS,eSW)
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
    print("ERROR: the chosen option is not available")
    return(NA)
  }
}
