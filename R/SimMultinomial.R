#' @title Simulation of multinomial clustered data
#'
#' @description Generate a random sample of multinomial observations from a copula-based mixed regression model.
#'
#' @param parC    copula parameters
#' @param parM    matrix of dimension (L-1)x k2 of margin parameters; L is the number of levels and k2 is the number of covariates+constant for the margins
#' @param clu     vector of clusters (can be a factor)
#' @param xc      matrix of covariates for the copula, not including the constant (can be NULL)
#' @param xm      matrix  of covariates for the margins, not including the constant (can be NULL)
#' @param family  copula family: "gaussian" (normal), "t" , "clayton" ,  "joe", "frank" , "gumbel", "plackett"
#' @param rot    rotation: 0 (default), 90, 180 (survival), or 270
#' @param dfC     degrees of freedom for student copula (default is NULL)
#' @param offset  offset for the margins (default is NULL)
#' @return \item{out}{List of simulated factors (y) and cluster factors (V)}
#'
#' @author  Bruno N. Remillard
#'
#' @examples
#' K=50 #number of clusters
#' n=5  #size of each cluster
#' N=n*K
#' set.seed(1)
#' clu=rep(c(1:K),each=n)
#' parC = 2
#' parM=matrix(c(1,-1,0.5,2),byrow=TRUE,ncol=2)
#' xm = runif(N)
#' y=SimMultinomial(parC,parM,clu,xm=xm,family="clayton",rot=90)$y
#' @export
#'
SimMultinomial<-function(parC,parM,clu,xc=NULL,xm=NULL,
                         family, rot=0,  dfC=NULL, offset=NULL)
{
  N=length(clu)
  K = max(clu)
  L= 1+dim(parM)[1]
  val=c(1:L)

  if(is.null(offset)){offset1=0}
  if(is.null(xc))
  {Matxc = matrix(1,nrow=N,ncol=1)}else
  {Matxc = cbind(1,xc)}
  if(is.null(xm))
  {Matxm = matrix(1,nrow=N,ncol=1)}else
  {Matxm = cbind(1,xm)}

  k1 = ncol(Matxc)
  k2 = ncol(Matxm)



  thC = colSums(parC*t(Matxc))

  thC0=linkCop(thC,family)$cpar

  ind0 = c(1:k2)
  par0=as.vector(t(parM))

  thF = matrix(0,ncol= L-1,nrow=N)
  for(j in 1:(L-1))
  {
    ind1 = (j-1)*k2+ind0
    thF[,j] = colSums(par0[ind1]*t(Matxm))+offset1
  }

  th1 = cbind(0,thF)
  p = exp(th1)
  p = p/rowSums(p)
  if(L==2){p[,1]=1-p[,2]}else{  p[,1]=1-rowSums(p[,2:L])}
  P = matrixStats::rowCumsums(p)
  P[,L]=1
  V=rep(0,N)
  for(k in 1:K){
    ind=(clu==k)
    v = runif(1)
    V[ind]=v
  }
  w = runif(N)
  U = rep(0,N)
  y = rep(0,N)
  for(i in 1:N){
    U[i]=qcond(w[i],V[i],family=family,thC0[i],rot)
    j=min(which(P[i,]>U[i]))
    y[i]=val[j]
  }
  y=as.factor(y)
  out=list(y=y,V=V)
  out
}

