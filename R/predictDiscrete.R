#' @title Conditional expectation for a copula-based estimation of mixed regression models for discrete response
#' @description Compute the conditional expectation of a copula-based  2-level hierarchical model for disctrete response.
#'
#' @param object  Object of class ``EstDiscrete`` generated by EstDiscrete.
#' @param newdata List of variables for be predicted (``clu`` for clusters, ``xc`` for the copula covariates, and ``xm`` for the margins covariates). The covariates can be NULL.
#' @param m   Number of points for the numerical integration (default is 100).
#'
#' @return \item{mest}{Conditional expectations (conditional probabilities for the multinomial case}
#' @references Krupskii, Nasri & Remillard (2023). On factor copula-based mixed regression models
#' @author Pavel Krupskii and Bruno N. Remillard, January 20, 2023
#' @import statmod, matrixStats
#'
#' @examples
#' data(out.poisson)
#' newdata = list(clu=c(1:50),xc=rep(0.2,50),xm=rep(0.5,50))
#' pred= predictDiscrete(out.poisson,newdata,m=100)
#' @export

predictDiscrete=function(object,newdata,m=100)
{

  family  = object$family
  dfC     = object$dfC
  V       = object$V
  par     = object$coefficients
  cluster = object$cluster
  disc  = object$disc
  rot     = object$rot

  clu  = newdata$clu
  nx   = length(clu)
  xc   = newdata$xc
  xm   = newdata$xm

  if(is.null(xc))
  {Matxc = matrix(1,nrow=nx,ncol=1)}else
  {Matxc = cbind(1,xc)}

  if(is.null(xm))
  {Matxm = matrix(1,nrow=nx,ncol=1)}else
  {Matxm = cbind(1,xm)}

  k1 = ncol(Matxc)
  k2 = ncol(Matxm)
  nx = nrow(Matxc)




  if(length(clu)!=nx){warning("cluster sizes do not match the total number of variables"); return(NULL) }

  # if(model=="bernoulli"){model=="binomial"}
  # switch(model,
  #        "binomial" = { disc= 1  },
  #
  #        "poisson" = {  disc = 2 },
  #
  #        "nbinom" = {   disc = 3 },
  #
  #        "geometric" = { disc = 4  },
  #
  #        "multinomial" = { disc = 5  }
  #
  # )


  thC = par$copula

  L1 = dim(par$margin)[1]

  thF = Matxm %*% t(par$margin)


  mest=matrix(0,nrow=nx,ncol=L1)
  for(i in 1:nx)
  {
    k = which(cluster ==clu[i])
    Matxck = Matxc[i,]
    thCk = sum(thC*Matxck)
    thFk = thF[i,]

    thC0=linkCop(thCk,family)$cpar

    switch(disc,
           {  p = 1/(1+exp(-thFk))
           mest[i] = 1 - pcond(1-p,V[k],family,rot,thC0,dfC)
           },
           { p = exp(thFk)
           mest[i] = sum( 1 - pcond(ppois(0:m,p),V[k],family,rot,thC0,dfC) )
           },
           {size = par[k1+k2+1];  p = 1/(1+exp(-thFk))
           mest[i] = sum( 1 - pcond(pnbinom(0:m,size,p),family,rot,thC0,dfC) )
           },
           { p = 1/(1+exp(-thFk))
           mest[i] = sum( 1 - pcond(pgeom(0:m,p),V[k],family,rot,thC0,dfC) )
           } ,
           {
             p = exp(thFk);
             p0 = 1/(1+sum(p))
             p = p0*p

             cump = p0+cumsum(p)
             cump[L1]=1
             cump0 = c(p0,cump[1:(L1-1)])
             mest[i,] = pcond(cump,V[k],family,rot,thC0,dfC)-pcond(cump0,V[k],family,rot,thC0,dfC)
           }
    )
  }
  return(mest)

}
