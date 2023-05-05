#' @title Estimation of latent variables in the continuous case
#'
#' @description This function computes the estimation of a latent variables for each cluster using the conditional a posteriori median.
#'
#' @param u      vector of values in (0,1)
#' @param family  copula family: "gaussian" , "t" , "clayton" ,  "joe", "frank" , "fgm", gumbel", "plackett", "galambos", "huesler-reiss"
#' @param rot     rotation: 0 (default), 90, 180 (survival), or 270.
#' @param thC0k   vector of copula parameters
#' @param dfC     degrees of freedome for the Student copula (default is NULL)
#' @param nq      number of nodes and weighted for Gaussian quadrature of the product of conditional copulas; default is 31.
#' @return \item{condmed}{Conditional a posteriori median.}
#'
#' @references Krupskii, Nasri & Remillard (2023). On factor copula-based mixed regression models
#' @author Pavel Krupskii, Bouchra R. Nasri and Bruno N. Remillard
#' @import  statmod, matrixStats
#' @examples
#' u = c(0.5228155, 0.3064417, 0.2789849, 0.5176489, 0.3587144)
#' thC0k=rep(17.54873,5)
#' MAP.continuous(u,"clayton",rot=90,thC0k,nq=35)
#' @export

MAP.continuous = function(u,family,rot,thC0k,dfC=NULL,nq=35){
  d=length(u)
  gl=statmod::gauss.quad.prob(nq)
  nl=gl$nodes
  wl=gl$weights
  uu=rep(u,nq)
  nn=rep(nl,each=d)

  param=rep(thC0k,nq)
  param=cbind(param,dfC);
  tem=dcop(uu,nn,family,rot,param)

  pdf=matrix(tem,ncol=nq)
  den=matrixStats::colProds(pdf)
  normC=sum(wl*den)
  cdff=function(x,thx){
    nnx=x*nn
    temx=dcop(uu,nnx,family,rot,param)
    pdfx=matrix(temx,ncol=nq)
    denx=matrixStats::colProds(pdfx)
    cdfx=x*sum(wl*denx)/normC
    cdfx
  }
  condmed=invfunc(0.5,cdff,0,lb=1e-10,ub=1-1e-10,tol=1e-8)
  return(condmed)
}
