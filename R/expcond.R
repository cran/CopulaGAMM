#' @title Conditional expectation
#'
#' @description This function computes the conditional expectation for a given copula family and a given margin variables for a clustered data model. The clusters ar3e independent but the observations with clusters are dependent, according to a one-factor copula model.
#'
#' @param w       value of the conditioning random variable
#' @param family  copula model: "gaussian" , "t" , "clayton" ,"joe", "frank" , "gumbel",  "plackett"
#' @param rot     rotation: 0 (default), 90, 180 (survival), or 270
#' @param cpar    copula parameter
#' @param margin  marginal distribution function
#' @param dfC     degrees of freedom for the Student copula (default is NULL)
#' @param subs    number of subdivisions for the integrals (default=1000)
#'
#' @author Pavel Krupskii and Bruno N. Remillard
#' @return \item{mest}{Conditional expectations}
#'
#' @examples
#' margin = function(x){ppois(x,10)}
#' expcond(0.4,'clayton',cpar=2,margin=margin)
#'
#' @export

expcond = function(w,family,rot=0,cpar,margin,dfC=NULL,subs=1000)
{
  nn=length(w)
  mest=numeric(nn)

  for(i in 1:nn){
    v=w[i]
    func1 = function(y){1-pcond(margin(y),v,family,rot,cpar,dfC)}
    func  = function(y){  pcond(margin(y),v,family,rot,cpar,dfC)}
    mest[i]=  integrate(func1,0,Inf,subdivisions=subs)$value- integrate(func,-Inf,0,subdivisions=subs)$value
  }
  return(mest)
}



