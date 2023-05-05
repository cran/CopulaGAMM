#' @title Inverse conditional expectation for a vector of probabilities
#'
#' @description This function computes the inverse conditional expecatation for a given copula family and a given margin variables for a clustered data model. The clusters ar3e independent but the observations with clusters are dependent, according to a one-factor copula model.
#'
#' @param u       conditional expectation
#' @param family  copula model: "gaussian" , "t" , "clayton"   "joe", "frank" , "gumbel",  "plackett"
#' @param cpar    copula parameter
#' @param rot     rotation: 0 (default), 90, 180 (survival), or 270
#' @param margin  marginal distribution function of the response
#' @param subs    number of subdivisions for the integrals (default=1000)
#' @param eps     precision required
#'
#' @return \item{minv}{Inverse conditional expectation}
#'
#' @author        Pavel Krupskii and Bruno N. Remillard
#' @export
#'

expcondinv = function(u,family,cpar,rot=0,margin,subs=1000,eps=1e-4)
{
n = length(u)
minv = numeric(n)
for(i in 1:n){minv[i] = expcondinv1(u[i],family,cpar,rot,margin,subs,eps)}
return(minv)
}

#' @title Inverse conditional expectation for a single value
#'
#' @param u       conditional expectation
#' @param family  copula model: "gaussian" , "t" , "clayton"   "joe", "frank" , "gumbel",  "plackett"
#' @param cpar    copula parameter
#' @param rot     rotation: 0 (default), 90, 180 (survival), or 270
#' @param margin  marginal distribution function of the response
#' @param subs    number of subdivisions for the integrals (default=1000)
#' @param eps     precision required
#'
#' @return \item{minv}{Inverse conditional expectation}
#'
expcondinv1 = function(u,family,cpar,rot=0,margin,subs=1000,eps=1e-4)
{

  a = 0.0001;
  b = 1-a;
  c1 = (a+b)/2
  f1 = expcond(a,family,cpar,rot,margin,subs)
  f2 = expcond(b,family,cpar,rot,margin,subs)

  f00=min(f1,f2)
  f01 = max(f1,f2)
  if( u< f00 || u> f01){ c1 = NA}
  else{
  f3 = expcond(c1,family,cpar,rot,margin,subs)
  s = f2-f1
  delta = abs(f3-u)
  if(s>0)
  {
    while(delta>eps){
      if(f3>u)
      {b=c1}else
      {
        a=c1}

      c1= (a+b)/2
      f1 = expcond(a,family,cpar,rot,margin,subs)
      f2 = expcond(b,family,cpar,rot,margin,subs)
      f3 = expcond(c1,family,cpar,rot,margin,subs)
      delta = abs(f3-u)
    }
  }else{
    while(delta>eps){
      if(f3>u)
      {a=c1}else
      {b=c1}

      c1= (a+b)/2
      f1 = expcond(a,family,cpar,rot,margin,subs)
      f2 = expcond(b,family,cpar,rot,margin,subs)
      f3 = expcond(c1,family,cpar,rot,margin,subs)
      delta = abs(f3-u)
    }

  }
  }
  return(c1)
}
