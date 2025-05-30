#' @title Estimation of latent variable in the dicrete case
#'
#' @description This function computes the estimation of a latent variables foe=r each cluster using the conditional a posteriori median.
#'
#' @param uu      vector of values in (0,1)
#' @param vv      vector of values in (0,1)
#' @param family  copula family "gaussian" , "t" , "clayton" ,  "joe", "frank" , "fgm", gumbel", "plackett", "galambos", "huesler-reiss"
#' @param rot     rotation: 0 (default), 90, 180 (survival), or 270.
#' @param thC0k   vector of copula parameters
#' @param dfC     degrees of freedom for the Student copula (default is NULL)
#' @param nq      number of nodes and weighted for Gaussian quadrature of the product of conditional copulas; default is 31.
#' @return \item{condmed}{Conditional a posteriori median.}
#'
#' @references Krupskii, Nasri & Remillard (2023). On factor copula-based mixed regression models
#' @author Pavel Krupskii, Bouchra R. Nasri and Bruno N. Remillard
#' @import  statmod, matrixStats
#' @examples
#' uu = c(0.5228155, 0.3064417, 0.2789849, 0.5176489, 0.3587144)
#' vv = c(0.7816627, 0.6688788, 0.6351364, 0.7774917, 0.7264787)
#' thC0k=rep(17.54873,5)
#' MAP.discrete(vv,uu,"clayton",rot=90,thC0k,nq=35)
#' @export



MAP.discrete = function(vv,uu,family,rot,thC0k,dfC=NULL,nq=35){
  d=length(vv)
  gl=statmod::gauss.quad.prob(nq)
  nl=gl$nodes
  wl=gl$weights
  diff = (abs(vv-uu) < 1e-5)
  uu[diff] = vv[diff] - 1e-5
  uu=pmax(uu,0)#new
  vvv=rep(vv,nq)
  uuu=rep(uu,nq)
  nn=rep(nl,each=d)
  param = rep(thC0k,nq)


  tem =pcond(vvv,nn,family,rot,param,dfC)
  tem1=pcond(uuu,nn,family,rot,param,dfC)


  pdf=matrix(tem-tem1,ncol=nq)
  pdf0=as.vector(pdf)
  pdf0=pmax(1e-100,pdf0)
  pdf=matrix(pdf0,ncol=nq)
  #den=matrixStats::colProds(adj*pdf)

  lden=colSums(log(pdf))
  lmax = max(lden)
  wprd = wl*exp(lden - lmax)
  normC = sum(wprd)

  #normC=sum(wl*den)
  cdff=function(x,thx){
    nnx=x*nn

    temx =pcond(vvv,nnx,family,rot,param,dfC)
    temx1=pcond(uuu,nnx,family,rot,param,dfC)

    pdfx=matrix(temx-temx1,ncol=nq)
    pdfx0=as.vector(pdfx)
    pdfx0=pmax(1e-100,pdfx0)
    pdfx=matrix(pdfx0,ncol=nq)
    #denx=matrixStats::colProds(adj*pdfx)

    lpdfx=colSums(log(pdfx))
    wprd = wl*exp(lpdfx - lmax)
    cdfx = x*sum(wprd)/normC

    #denx=matrixStats::colProds(pdfx)
    #cdfx=x*sum(wl*denx)/normC
    return(cdfx)
  }
  condmed=invfunc(0.5,cdff,0,lb=1e-10,ub=1-1e-10,tol=1e-8)
  return(condmed)
}
