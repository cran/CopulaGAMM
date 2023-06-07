#' @title Margins  cdf/pdf and their derivatives
#' @description This function computes the cdf, pdf, and associated derivatives
#' @param z     vector of responses
#' @param th    linear combination of covariates (can be negative)
#' @param model model for margin: "binomial" (bernoulli), "poisson", "nbinom" (mean is the parameter),"nbinom1" (p is the parameter), "geometric", "multinomial", "exponential", "weibull", "normal","t", "laplace"
#' @param x     covariates for the multinomial margin (default is NULL)
#' @param dfM   degrees of freedom for the Student margin (default is NULL)
#' @author Pavel Krupskii and Bruno N. Remillard, January 20, 2022
#' @return \item{out}{Matrix of conditional cdf, derivative with respect to parameter, pdf, }
#' @examples
#' out = margins(0,2.5,"binomial")
#' @export

margins = function(z,th,model,x=NULL,dfM=NULL){
  if(model=="bernoulli"){model=="binomial"}
  switch(model,
         "binomial" = {

           out= berncpdf(z,th)
         },

         "poisson" = {
           out = poiscpdf(z,th)
         },

         "nbinom" = {
           out = nbinomcpdf(z,th)

         },

         "nbinom1" = {
           out = nbinom1cpdf(z,th)

         },

         "geometric" = {
           out = geomcpdf(z,th)

         },

         "multinomial" = {

           out = multinomcpdf(z,th,x)
         },

         "exponential" = {
           out = expcpdf(z,th)
         },

         "weibull" = {
           out = weibcpdf(z,th)
         },

         "normal" = {
           out = normcpdf(z,th)
         },

         "t" = {
           out = tcpdf(z,th,dfM)
         },

         "laplace" = {
           out = lapcpdf(z,th)
         }

  )
  return(out)

}


#' @title Bernoulli with p = 1/(1+exp(-th)) cdf/pdf and derivatives
#' @description This function computes the cdf, pdf, and associated derivatives
#' @param z vector of responses
#' @param th linear combination of covariates (can be negative)
#'
#' @author Pavel Krupskii and Bruno N. Remillard, January 20, 2022
#' @return \item{out}{Matrix of conditional cdf and pdf with derivative with respect to parameters}
#' @examples
#' out = berncpdf(0,2.5)
#'
#' @export
#'
berncpdf=function(z,th){
  n=length(z)
  p = 1/(1+exp(-th))
  ind1 = (z==1)
  ind0 = (z==0)
  pdf  = vector(mode='numeric',n)
  cdf  = vector(mode='numeric',n)
  pdf1 = vector(mode='numeric',n)
  cdf1 = vector(mode='numeric',n)
  pdf[ind1]  = p
  pdf[ind0]  = 1-p
  cdf[ind1]  = 1
  cdf[ind0]  = 1-p
  pdf1[ind1] =  p*(1-p)
  pdf1[ind0] = -p*(1-p)
  cdf1[ind0] = pdf1[ind0]

  cbind(cdf,pdf,cdf1,0,pdf1,0)
}


#' @title Poisson cdf/pdf and ders
#' @description This function computes the cdf, pdf, and associated derivatives
#' @param z vector of responses
#' @param th values of lambda >0
#'
#' @author Pavel Krupskii and Bruno N. Remillard, January 20, 2022
#' @return \item{out}{Matrix of conditional cdf, derivative with respect to parameter, pdf, }
#' @examples
#' out = poiscpdf(0,2.5)
#' @export
#'
poiscpdf=function(z,th){
  cdf=ppois(z,th)
  pdf= dpois(z,th)
  cdf1=-pdf
  pdf1= dpois(z-1,th)-dpois(z,th)
  cbind(cdf,pdf,cdf1,0,pdf1,0)
}

#' @title Negative binomial cdf/pdf and ders
#' @description This function computes the cdf, pdf, and associated derivatives
#' @param z vector of responses
#' @param th th[,1] is size > 0 and th[,2] is p, with 0<p<1; size  does not have to be integer
#'
#' @author Pavel Krupskii and Bruno N. Remillard, January 20, 2022
#' @return \item{out}{Matrix of conditional cdf, derivative with respect to parameter, pdf, }
#' @examples
#' out = nbinomcpdf(0,c(1,0.5))
#' @export
#'
nbinomcpdf=function(z,th){
  n = length(z)
  if(is.vector(th)) th  = matrix(th, ncol=2, nrow=n, byrow=TRUE)
  size=th[,1]
  p = th[,2]

  i0 = z < 0
  z = abs(z)

  cdf=pnbinom(z,size=size,prob=p)
  pdf=dnbinom(z,size=size,prob=p)

  pdf2 = pdf*(size/p - z/(1-p))
  pdf1 = pdf*( log(p) +digamma(z+size)-digamma(size))

  tem1 = rep(0,n)
  tem2 = rep(0,n)
  for(i in 1:n){
    pdfi = dnbinom(0:z[i],size[i],prob=p[i])
    tem1[i] = sum(pdfi[-1]*(1:z[i]))
    tem2[i] = sum(pdfi*digamma(0:z[i] + size[i]))
  }

  cdf2=cdf*size/p - tem1/(1-p)
  cdf1=(log(p)-digamma(size))*cdf + tem2

  out=cbind(cdf,pdf,cdf1,cdf2,pdf1,pdf2)
  out[i0,] = 0
  out
}
#' @title Negative binomial cdf/pdf and ders
#' @description This function computes the cdf, pdf, and associated derivatives
#' @param z vector of responses
#' @param th th[,1] is size > 0 and th[,2] is mean > 0; size  does not have to be integer
#'
#' @author Pavel Krupskii and Bruno N. Remillard, January 20, 2022
#' @return \item{out}{Matrix of conditional cdf, derivative with respect to parameter, pdf, }
#' @examples
#' out = nbinom1cpdf(0,c(1,0.5))
#' @export
#'

nbinom1cpdf=function(z,th){
  n = length(z)
  if(is.vector(th)) th  = matrix(th, ncol=2, nrow=n, byrow=TRUE)
  size=th[,1]
  mu = th[,2]

  p = size/(size+mu)

  i0 = z < 0
  z = abs(z)

  cdf=pnbinom(z,size=size,mu=mu)
  pdf=dnbinom(z,size=size,mu=mu)


  pdf1=pdf*(digamma(z+size)-digamma(size)+(mu-z)/(size+mu) + log(p))
  pdf2=pdf*p*(-1+z/mu)

  tem1 = rep(0,n)
  tem2 = rep(0,n)
  for(i in 1:n){
    pdfi = dnbinom(0:z[i],size[i],mu=mu[i])
    tem1[i] = sum(pdfi[-1]*(1:z[i]))
    tem2[i] = sum(pdfi*digamma(0:z[i] + size[i]))
  }

  cdf2=-p*cdf + (p/mu)*tem1
  cdf1=tem2+(1-p+log(p)-digamma(size))*cdf - tem1/(mu+size)

  out=cbind(cdf,pdf,cdf1,cdf2,pdf1,pdf2)
  out[i0,] = 0
  out
}

#' @title Geometric with p = 1/(1+exp(-th)) cdf/pdf and ders
#' @description This function computes the cdf, pdf, and associated derivatives
#' @param z vector of responses
#' @param th linear combination of covariates (can be negative)
#'
#' @author Pavel Krupskii and Bruno N. Remillard, January 20, 2022
#' @return \item{out}{Matrix of conditional cdf, derivative with respect to parameter, pdf, }
#' @examples
#' out = geomcpdf(0,-3)
#' @export
#'
geomcpdf=function(z,th){

  p = 1/(1+exp(th))

  cdf  = pgeom(z,p)
  pdf  = dgeom(z,p)
  pdf1 = pdf*(p*z-1+p)# there was an error here pdf*(2*p-1)
  cdf1 = -pdf*(1+z)*(1-p)


  cbind(cdf,pdf,cdf1,0,pdf1,0)
}

#' @title Multinomial with p = 1/(1+exp(-th)) cdf/pdf and ders
#' @description This function computes the cdf, pdf, and associated derivatives
#' @param z vector of responses taking values in {1,...,nL}: as.number(z) if z is a factor!
#' @param th th is a n x (L-1) matrix of parameters, i.e., mpar = a=[a_{1,1},...a_{1,k2},a_{2,1},...a_{2,k2},... a_{L-1,1}... a_{L-1,k2}], and first level is the baseline.
#' @param x matrix of covariates (including the constant)
#' @author Pavel Krupskii and Bruno N. Remillard, January 20, 2022
#' @return \item{out}{Matrix of conditional cdf, derivative with respect to parameter, pdf, }
#' @examples
#' x=matrix(c(1,1,-1,-1,0,2),nrow=2)
#' z = c(1,3)
#' th = matrix(c(1,2,3,4,5,6),nrow=2)
#' out = multinomcpdf(z,th,x = x)
#' @export
#'
multinomcpdf=function(z,th,x){
 if(is.vector(th)){th=matrix(th,ncol=1)}
  if(is.vector(x)){x=matrix(x,ncol=1)}
   th1 = cbind(0,th)
  p = exp(th1)
  L = dim(p)[2]
  p = p/rowSums(p)
  if(L==2){p[,1]=1-p[,2]}else{p[,1]=1-rowSums(p[,2:L])}

  P = matrixStats::rowCumsums(p)
  P[,L]=1
  n=length(z)

  k2 = dim(x)[2] #k2
  ind1 = c(1:n)
  ind0 = (z>0);
  pdf  = rep(0,n);

  cdf  = rep(0,n);

  Gpdf = matrix(0,c(1:n*k2*(L-1)),nrow=n);
  pdf1=Gpdf
  cdf1=Gpdf
  Gcdf =Gpdf
  A0= Gpdf
  G1= Gpdf
  G2=G1
  p0 = p[ind0,]
  z0 = z[ind0];
  x0 = x[ind0,];
  if(is.vector(x0)){x0=matrix(x0,ncol=1)}


  ind = cbind(ind1[ind0],z0)
  ind00 = cbind(ind1[ind0],z0)
  f0 = p[ind00]
  F0 = P[ind00]
  pdf[ind0]  = f0
  cdf[ind0]  = F0

  p00=p0[,-1]
  if(L==2){
    p00=as.matrix(p00,ncol=1)
  }


  A = wprod(p00,x0)
  Gcdf[ind0,]= F0*A
  Gpdf[ind0,]= f0*A
  A0[ind0,] = A

  indk2 = c(1:k2)
  for(i in 1:n){
    if(z[i]>1)
    {
      m = (z[i]-2)*k2+indk2 ;
      G1[i,m] = A0[i,m];
    }
  }

  pdf1 = G1-Gpdf;

  for(i in 1:n){
    if(z[i]>1)
    {
      m = c(1: ((z[i]-1)*k2));
      G2[i,m] = A0[i,m];
    }
  }

  cdf1 = G2-Gcdf;
  # ind=(z==L)
  # cdf1[ind,]=0

  cbind(cdf,pdf,cdf1,0,pdf1,0)
}



wprod<-function(M,x)
{
  # M: n x L; L>1
  # x: n x k

  L = dim(M)[2];
  A = M[,1]*x;
  if(L>1)
    {
    for(j in 2:L)
      {
       A = cbind(A,M[,j]*x);  #(x * M[,1]| ... | x * M[,L])
    }
  }
  return(A) # n x (Lxk) matrix

}

#' @title Exponential cdf/pdf and ders
#' @description This function computes the cdf, pdf, and associated derivatives
#' @param z vector of responses
#' @param th th is rate > 0
#'
#' @author Pavel Krupskii and Bruno N. Remillard, January 20, 2022
#' @return \item{out}{Matrix of conditional cdf, derivative with respect to parameter, pdf, }
#' @examples
#' out = expcpdf(2,3)
#' @export
#'
expcpdf=function(z,th){
  if(is.vector(th)){ th = matrix(th,ncol=1)}
  th1=th[,1]
  ex=exp(-z*th1)
  cdf=1-ex
  pdf=th1*ex
  cdf1=z*ex
  pdf1=(1-z*th1)*ex
  cbind(cdf,pdf,cdf1,0,pdf1,0)
}

#' @title Weibul cdf/pdf and ders
#' @description This function computes the cdf, pdf, and associated derivatives
#' @param z vector of responses
#' @param th th[,1] is rate>0, th[,2] is shape  > 0;
#'
#' @author Pavel Krupskii and Bruno N. Remillard, January 20, 2022
#' @return \item{out}{Matrix of conditional cdf, derivative with respect to parameter, pdf, }
#' @examples
#' out = weibcpdf(2,c(2,3))
#' @export
#'
weibcpdf=function(z,th){
  if(is.vector(th)){ la=th[1]; al=th[2] }
  if(is.matrix(th)){ la=th[,1]; al=th[,2] }
  #la=th[1]
  #al=th[2]
  za1= (z*la)^(al-1)
  za = za1*z*la
  lx = log(la*z)
  ex=exp(-za)
  cdf=1-ex
  pdf=al*la*za1*ex
  cdf1=(al/la)*za*ex
  cdf2=za*lx*ex
  pdf1=al*al*za1*(1-za)*ex
  pdf2=za1*(la-al*la*za*lx+al*la*lx)*ex
  cbind(cdf,pdf,cdf1,cdf2,pdf1,pdf2)
}

#' @title normal cdf/pdf and ders
#' @description This function computes the cdf, pdf, and associated derivatives
#' @param z vector of responses
#' @param th th[,1] is mean, th[,2] is standard deviation  > 0;
#'
#' @author Pavel Krupskii and Bruno N. Remillard, January 20, 2022
#' @return \item{out}{Matrix of conditional cdf, derivative with respect to parameter, pdf, }
#' @examples
#' out = normcpdf(2,c(-3,4))
#' @export
#'
normcpdf=function(z,th){
  if(is.vector(th)){ m=th[1]; s=th[2] }
  if(is.matrix(th)){ m=th[,1]; s=th[,2] }
  stz = (z-m)/s
  stz2= stz*stz
  cdf=pnorm(z,mean=m,sd=s)
  pdf=exp(-.5*stz2)/sqrt(2*pi)/s
  cdf1=-pdf; cdf2=-stz*pdf
  pdf1=pdf*stz/s; pdf2=-pdf/s+stz2*pdf/s
  cbind(cdf,pdf,cdf1,cdf2,pdf1,pdf2)
}

#' @title Student  cdf/pdf and ders
#' @description This function computes the cdf, pdf, and associated derivatives
#' @param z vector of responses
#' @param th th[,1] is mean, th[,2] is standard deviation  > 0
#' @param df degrees of freedom
#'
#' @author Pavel Krupskii and Bruno N. Remillard, January 20, 2022
#' @return \item{out}{Matrix of conditional cdf, derivative with respect to parameter, pdf, }
#' @examples
#' out = tcpdf(2,c(-3,4),25)
#' @export
#'
tcpdf=function(z,th,df){
  if(is.vector(th)){ m=th[1]; s=th[2] }
  if(is.matrix(th)){ m=th[,1]; s=th[,2] }
  df=df[1]
  stz = (z-m)/s
  stz2= stz*stz/df
  cdf=pt(stz,df=df)
  C=gamma((df+1)/2)/sqrt(df*pi)/gamma(df/2)
  pdf=C*(1+stz2)^(-(df+1)/2)/s
  cdf1=-pdf
  cdf2=-pdf*stz
  C0=(1+1/df)/(1+stz2)
  pdf1=C0*pdf*stz/s
  pdf2=-pdf/s+stz*pdf1
  cbind(cdf,pdf,cdf1,cdf2,pdf1,pdf2)
}


#' @title Laplace cdf/pdf and ders
#' @description This function computes the cdf, pdf, and associated derivatives
#' @param z vector of responses
#' @param th th[,1] is mean, th[,2] is standard deviation  > 0
#'
#' @author Pavel Krupskii and Bruno N. Remillard, January 20, 2022
#' @return \item{out}{Matrix of conditional cdf, derivative with respect to parameter, pdf, }
#' @examples
#' out = lapcpdf(2,c(-3,4))
#' @export
#'
lapcpdf=function(z,th){
  if(is.vector(th)){ m=th[1]; s=th[2] }
  if(is.matrix(th)){ m=th[,1]; s=th[,2] }
  stz = (z-m)/s

  cdf=plap(stz)
  pdf=dlap(stz)/s

  cdf1=-pdf;
  cdf2= stz*cdf1;

  pdf1= sign(stz)*pdf/s;
  pdf2=-pdf/s+stz*pdf1;
  out=cbind(cdf,pdf,cdf1,cdf2,pdf1,pdf2)
  return(out)
}

#Inverse Laplace cdf
qlap = function(u){
  ind0=(u<=0.5)
  ind1 = (u>0.5)
  n = length(u)
  q=rep(0,n)
  q[ind0] = log(2*u[ind0])
  q[ind1] = -log(2*(1-u[ind1]))
  return(q)
}

# Laplace cdf
plap = function(x){
  ind0=(x<=0)
  ind1 = (x>0)
  n = length(x)
  f=rep(0,n)
  f[ind0] = 0.5*exp(x[ind0])
  f[ind1] = 1-0.5*exp(-x[ind1])
  return(f)
}

# Laplace pdf
dlap = function(x){
  return(0.5*exp(-abs(x)))
}
