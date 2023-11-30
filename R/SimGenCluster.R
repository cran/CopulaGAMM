#' @title Simulation of clustered data
#'
#' @description Generate a random sample of observations from a copula-based mixed regression model.
#'
#' @param parC    vector of copula parameters; k1  is the number of covariates + constant for the copula
#' @param parM    vector of margin parameters; k2  is the number of covariates + constant for the margins
#' @param clu     vector of clusters (can be a factor)
#' @param xc      matrix (N x k1) of covariates for the copula, not including the constant (can be NULL)
#' @param xm      matrix (N x k2) of covariates for the margins, not including the constant (can be NULL)
#' @param family  copula family: "gaussian" , "t" , "clayton" ,  "joe", "frank" , "gumbel", "plackett"
#' @param rot     rotation: 0 (default), 90, 180 (survival), or 270
#' @param dfC     degrees of freedom for the Student copula (default is NULL)
#' @param model   marginal distribution: "binomial" (bernoulli), "poisson", "nbinom" (mean is the parameter),"nbinom1" (p is the parameter), "geometric", "multinomial", exponential", "weibull", "normal" (gaussian),"t", "laplace"
#' @param dfM     degrees of freedom for the Student margins (default is NULL)
#' @param offset  offset for the margins (default is NULL)
#'
#' @return \item{y}{Simulated response}
#'
#' @author  Bruno N. Remillard
#' @return \item{y}{Simulated values}
#'
#' @examples
#' K=50 #number of clusters
#' n=5  #size of each cluster
#' N=n*K
#' set.seed(1)
#' clu=rep(c(1:K),each=n)
#' parC = 0 # yields tau = 0.5 for Clayton
#' parM= c(1,-1,4)
#' xm = runif(N)
#' y=SimGenCluster(parC,parM,xm,family="clayton",rot=90,clu=clu,model="gaussian")
#' @export
#'
SimGenCluster<-function(parC, parM, clu, xc=NULL, xm=NULL,
                        family, rot=0, dfC=NULL,
                        model, dfM=NULL, offset=NULL)
{

  if(model=="multinomial"){
    y = SimMultinomial(parC,parM,clu,xc,xm,family,rot,dfC,offset)
  }else{
   N=length(clu)
  K = max(clu)

  if(model=="gaussian"){model="normal"}
  if(model=="bernoulli"){model="binomial"}

  if(is.null(offset)){offset1=0}
  if(is.null(xc))
  {Matxc = matrix(1,nrow=N,ncol=1)}else
  {Matxc = cbind(1,xc)}
  if(is.null(xm))
  {Matxm = matrix(1,nrow=N,ncol=1)}else
  {Matxm = cbind(1,xm)}

  k1 = ncol(Matxc)
  k2 = ncol(Matxm)


  thC = matrixStats::colSums2(parC*t(Matxc))
  cpar = linkCop(thC,family)$cpar
  size=NULL

  thF = matrixStats::colSums2(parM[1:k2]*t(Matxm))+offset1

  if(length(parM)>k2)
    {
       size = parM[(k2+1)]
    }




  V=rep(0,N)
  for(k in 1:K){
    ind=(clu==k)
    v = runif(1)
    V[ind]=v
  }

  w = runif(N)
  U = rep(0,N)
  y = rep(0,N)

  if(family=="t") cpari = cbind(cpar,dfC)

  U=qcond(w,V,family=family,cpar,rot)

  #for(i in 1:N){
  #  cpari = cpar[i]
  #  if(family=="t") cpari = c(cpari,dfC)
  #  U[i]=qcond(w[i],V[i],family=family,cpari,rot)
  #  }


  switch(model,
         "binomial" = {
           p = 1/(1+exp(-thF)) #P(Y=1)
           y= as.numeric(U>1-p)
         },

         "poisson" = {
           mu  = exp(thF)
           y = qpois(U,mu)
         },

         "nbinom" = {
           mu  = exp(thF)
           out = qnbinom(U,size=size,mu=mu)
           },

         "nbinom1" = {
           p = 1/(1+exp(-thF))
           out =qnbinom(U,size=size,prob=p)
           },

         "geometric" = {
           p = 1/(1+exp(-thF))
           out = qgeom(U,p)
           },


         "exponential" = {
           rate = exp(thF)
           out = -rate*log(1-U)
         },

         "weibull" = {
           shape  = exp(thF)
           out = qweibull(U,shape,scale=size)
         },

         "normal" = {
           y = thF+size*qnorm(U)
         },

         "t" = {
           y = thF+size*qt(U,dfM)
         },

         "laplace" = {
           y = thF+size*qlap(U)
         }


  )
  }
  return(y)

}

