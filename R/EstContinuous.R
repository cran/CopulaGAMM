#' @title Copula-based estimation of mixed regression models for continuous response
#'
#' @description This function computes the estimation of a copula-based  2-level hierarchical model.
#'
#' @param y       n x 1 vector of response variable (assumed continuous).
#' @param model   function for margins: "gaussian" (normal), "t" (Student with known df=dfM), "laplace" , "exponential", "weibull".
#' @param family  copula family: "gaussian" (normal), "t" , "clayton" ,  "frank" , "fgm", "gumbel".
#' @param rot     rotation: 0 (default), 90, 180 (survival), or 270
#' @param clu     variable of size n defining the clusters; can be a factor
#' @param xc      covariates of size n for the estimation of the copula, in addition to the constant; default is NULL.
#' @param xm      covariates of size n for the estimation of the mean of the margin, in addition to the constant; default is NULL.
#' @param start   starting point for the estimation; default (NULL) are the ones associated with a Gaussian-copula model defined by lme.
#' @param LB      lower bound for the parameters.
#' @param UB      upper bound for the parameters.
#' @param nq      number of nodes and weighted for Gaussian quadrature of the product of conditional copulas; default is 25.
#' @param dfC     degrees of freedom for a Student margin; default is 5.
#' @param dfM     degrees of freedom for a Student margin; default is NULL for non-t distribution.
#' @param prediction  logical variable for prediction of latent variables V; default is TRUE.
#' @author Pavel Krupskii, Bouchra R. Nasri and Bruno N. Remillard
#' @return \item{coefficients}{List of estimated parameters: copula, margin, size}
#' @return \item{sd}{Standard deviations of the estimated parameters}
#' @return \item{tstat}{T statistics for the estimated parameters}
#' @return \item{pval}{P-values of the t statistics for the estimated parameters}
#' @return \item{gradient}{Gradient of the log-likelihood}
#' @return \item{loglik}{Log-likelihood}
#' @return \item{aic}{AIC coefficient}
#' @return \item{bic}{BIC coefficient}
#' @return \item{cov}{Covariance matrix of the estimations}
#' @return \item{grd}{Gradients by clusters}
#' @return \item{clu}{Cluster values}
#' @return \item{Matxc}{Matrix of covariates defining the copula parameters, including a constant}
#' @return \item{Matxm}{Matrix of covariates defining the margin parameters, including a constant}
#' @return \item{V}{Estimated value of the latent variable by clusters (if prediction=TRUE)}
#' @return \item{cluster}{Unique values of clusters}
#' @return \item{family}{Copula family}
#' @return \item{tau}{Kendall's tau by observation}
#' @return \item{thC0}{Estimated parameters of the copula by observation}
#' @return \item{thF}{Estimated parameters of the margins by observation}
#' @return \item{pcond}{Conditional copula cdf}
#' @return \item{fcpdf}{Margin functions (cdf and pdf)}
#' @return \item{dfM}{Degrees of freedom for Student margin (default is NULL)}
#' @return \item{dfC}{Degrees of freedom for the Student copula (default is NULL)}
#' @references Krupskii, Nasri & Remillard (2023). On factor copula-based mixed regression models
#'
#' @import  statmod, matrixStats
#' @examples
#' data(normal) #simulated data with normal margins
#' start=c(0,0,0,1); LB=c(rep(-10,3),0.001);UB=c(rep(10,3),10)
#' y=normal$y; clu=normal$clu;xm=normal$xm
#' out=EstContinuous(y,model="gaussian",family="clayton",rot=90,clu=clu,xm=xm,start=start,LB=LB,UB=UB)
#' @export

EstContinuous =function(y,model,family,rot=0,clu,xc=NULL,xm=NULL,start=NULL, LB=NULL, UB=NULL, nq=31,dfC=NULL,dfM=NULL,prediction=TRUE){


  if(model=="gaussian"){model="normal"}
  d = length(y)
  if(length(clu)!=d){warning("clusters do not match the total number of variables"); return(NULL) }
  if(is.null(xc))
  {Matxc = matrix(1,nrow=d,ncol=1)}else
  {Matxc = cbind(1,xc)}

  if(is.null(xm))
  {Matxm = matrix(1,nrow=d,ncol=1)}else
  {Matxm = cbind(1,xm)}

  k1 = ncol(Matxc)
  k2 = ncol(Matxm)
  if(is.factor(clu))
  {cluster=levels(clu)}else
  {cluster=unique(clu)}  # all possible cluster values

  nclu = length(cluster) # number of clusters

  ## Automated values if initial values or upper/lower bounds are NULL
  if(is.null(start) || is.null(LB) || is.null(UB)){
    mmle = EstContinuousMnormal(y,xm=xm)
    mstart = mmle$coefficients
    cstart = c(-log(1.5),rep(0,k1-1))
    start = c(cstart,mstart)
    LB = start - 1000
    UB = start + 1000
    LB[k1+k2+1] = 1e-10
    }

  likf=function(par){
    gl=statmod::gauss.quad.prob(nq)
    nl=gl$nodes
    wl=gl$weights
    if(min(par - LB) < 0 || max(par - UB) > 0) return(1e100)

    thC = colSums(par[1:k1]*t(Matxc))  #phi(beta,x_ki)
    thF = cbind(colSums(par[k1+(1:k2)]*t(Matxm)), par[k1+k2+1]) #h(alpha,x_ki)

    cpdf= margins(y,thF,model,dfM)

    grd = matrix(0,nrow=nclu,ncol=(k1+k2+1))  # needed to estimate the covariance matrix
    # grdsum = rep(0,k1+k2+1)

    out = 0
    fk  = rep(0,nclu)

    for(k in 1:nclu){
      ind = (clu==cluster[k])
      n_k = sum(ind) # number of elements in the cluster


      Matxck = Matxc[ind,]
      if(is.vector(Matxck)){Matxck=matrix(Matxck,ncol=1)}
      nclk = clu[ind]
      Matxmk = Matxm[ind,]
      if(is.vector(Matxmk)){Matxmk=matrix(Matxmk,ncol=1)}
      cpdfk=cpdf[ind,]
      u   = cpdfk[,1]  #G(y_ki)
      u1  = cpdfk[,3]  #cdf1
      u2  = cpdfk[,4]  #cdf2
      mpdfk  = cpdfk[,2]  #pdf
      mpdf1k = cpdfk[,5]  #pdf1
      mpdf2k = cpdfk[,6]  #pdf2
      uu  = rep(u,nq)
      nn  = rep(nl,each=n_k)
      thCk = thC[ind]
      out0=linkCop(thCk,family)
      thC0 = out0$cpar
      thC1 = rep(thC0,nq)
      thCd = out0$hder*Matxck
      tem = coplik(uu,nn,family,rot,thC1,dfC,TRUE)

      tem3 = tem[,3] #c_ki
      tem4 = tem[,4] #\dot c_ki
      tem5 = tem[,5] #\partial_u \dot c_ki
      tem3[tem3<1e-100] = 1e-100

      mtem3 = matrix(tem3,ncol=nq)
      ltem3 = colSums(log(mtem3))
      lmax3 = max(ltem3)
      wprd = wl*exp(ltem3 - lmax3)
      intf = sum(wprd)


      sl = sum(log(mpdfk))+log(intf) #log(f_k)
      out = out - sl - lmax3

      M=matrix(tem4/tem3,nrow=nq,byrow=TRUE)  #\dot c_ki/c_ki
      M1 = matrix(tem5/tem3,nrow=nq,byrow=TRUE)
      grd[k,1:k1]      =   - colSums(thCd*colSums(wprd*M))/intf
      grd[k,k1+(1:k2)] =   - colSums( wprd*M1%*%(u1*Matxmk) )/intf - colSums(Matxmk*mpdf1k/mpdfk)
      grd[k,k1+k2+1]   =  - sum( u2*t(wprd*matrix(tem5/tem3,nrow=nq,byrow=TRUE)) )/intf-sum(mpdf2k/mpdfk)

    }

    grdsum = colSums(grd)
    attr(out, "gradient") = grdsum
    attr(out, "grd") = grd
    out
  }
  mle = nlm(likf,p=start,check.analyticals=F,print.level=0,iterlim=200)

  V=NULL



  par = mle$estimate


  thC = colSums(par[1:k1]*t(Matxc))  #phi(beta,x_ki)

  out0=likf(par)
  grd=attributes(out0)$grd
  thC0 = linkCop(thC,family)$cpar

  thF = cbind(colSums(par[k1+(1:k2)]*t(Matxm)), par[k1+k2+1])
  if(prediction)
  {
    ## prediction of V using the posterior median
    V=rep(0,nclu)

    cpdf= margins(y,thF,model,dfM)



    for(k in 1:nclu){

      ind = (clu==cluster[k])
      n_k = sum(ind) # number of elements in the cluster
      Matxck = Matxc[ind,]
      thC0k=thC0[ind]
      u   = cpdf[ind,1]


      V[k] = MAP.continuous(u,family,rot,thC0k,dfC,nq)

    }
  }

  covar = cov(grd)
  st.dev = sqrt(diag(covar)/nclu)
  tstat = mle$estimate/st.dev
  pval = 2*pnorm(abs(tstat),lower.tail = FALSE)
  parC=mle$estimate[1:k1]
  parM=mle$estimate[(k1+1):(k1+k2)]
  size=NULL
  k=length(mle$estimate)
  if(k>(k1+k2)){size=mle$estimate[k1+k2+1]}
  k = length(size)+length(parC)+length(parM)
  LL = mle$minimum #-log-likelihood
  AIC = 2*k+2*LL
  BIC = k*log(d)+2*LL
  coef=list(copula=parC,margin=parM,size=size)
  out=list(coefficients=coef, sd = st.dev, tstat=tstat, pval=pval,
           gradient=mle$gradient,loglik=-LL, aic=AIC,bic=BIC,cov=covar,
           grd=grd,clu=clu,Matxc=Matxc,Matxm=Matxm,V=V,cluster=cluster,
           family=family,thC0=thC0,thF=thF, rot=rot,model=model,
           dfC=dfC,dfM=dfM,disc=NULL)
 # class(out) <-  "EstContinuous"
  out

}



#MARGINAL LIKELIHOOD IN THE GENERAL CASE
EstContinuousM =function(y,model,clu,xm=NULL,start, dfM=NULL){

  if(model=="gaussian"){model="normal"}
  d = length(y)
  if(length(clu)!=d){warning("clusters do not match the total number of variables"); return(NULL) }


  if(is.null(xm))
  {Matxm = matrix(1,nrow=d,ncol=1)}else
  {Matxm = cbind(1,xm)}

  k2 = ncol(Matxm)
  if(is.factor(clu))
  {cluster=levels(clu)}else
  {cluster=unique(clu)}  # all possible cluster values

  nclu = length(cluster) # number of clusters


  likf=function(par){

    #LB = -100
    #UB =  100
    #if(min(par - LB) < 0 || max(par - UB) > 0) return(1e100)

    thF = cbind(matrixStats::colSums2(par[(1:k2)]*t(Matxm)), par[k2+1]) #h(alpha,x_ki)

    cpdf= margins(y,thF,model,dfM)
    mpdf  = cpdf[,2]  #pdf
    mpdf1 = cpdf[,5]  #pdf1
    mpdf2 = cpdf[,6]  #pdf2


    grd = matrix(0,nrow=nclu,ncol=(k2+1))  # needed to estimate the covariance matrix

    out = 0
    fk  = rep(0,nclu)

    for(k in 1:nclu){
      ind = (clu==cluster[k])
      n_k = sum(ind) # number of elements in the cluster



      nclk = clu[ind]
      Matxmk = Matxm[ind,]
      if(is.vector(Matxmk)){Matxmk=matrix(Matxmk,ncol=1)}
      cpdfk=cpdf[ind,]
      u   = cpdfk[,1]  #G(y_ki)
      u1  = cpdfk[,3]  #cdf1
      u2  = cpdfk[,4]  #cdf2
      mpdfk  = cpdfk[,2]  #pdf
      mpdf1k = cpdfk[,5]  #pdf1
      mpdf2k = cpdfk[,6]  #pdf2


      out = out-sum(log(mpdfk))
      grd[k,(1:k2)] =  - matrixStats::colSums2(Matxmk*mpdf1k/mpdfk)
      grd[k,k2+1]   =  - sum(mpdf2k/mpdfk)

    }

    grdsum = colSums(grd)
    attr(out, "gradient") = grdsum
    attr(out, "grd") = grd
    out
  }
  mle = nlm(likf,p=start,check.analyticals=F,print.level=0,iterlim=500)

  par = mle$estimate
  k = length(par)
  LL = mle$minimum #-log-likelihood
  AIC = 2*k+2*LL
  BIC = k*log(d)+2*LL

  out0=likf(par)

  grd=attributes(out0)$grd

  thF = cbind(matrixStats::colSums2(par[(1:k2)]*t(Matxm)), par[k2+1])

  u = margins(y,thF,model,dfM=dfM)[,1]

  out=list(coefficients=mle$estimate, loglik=-LL, aic=AIC,bic=BIC,
           clu=clu,Matxm=Matxm,cluster=cluster,
           model=model,thF=thF, u=u, dfM=dfM,disc=NULL)
  class(out) <-  "EstContinuous"
  out

}


#MARGINAL LIKELIHOOD WITH NORMAL MARGINALS
EstContinuousMnormal =function(y,xm=NULL){

  d = length(y)

  if(is.null(xm))
  {Matxm = matrix(1,nrow=d,ncol=1)}else
  {Matxm = cbind(1,xm)}


  pr = solve(t(Matxm)%*%Matxm,t(Matxm)%*%y)
  thF = Matxm%*%pr
  k1 = length(pr)
  err = y-thF
  sig2= mean(err*err)
  sig = sqrt(sig2)
  LL = -sum(dnorm(err,0,sig,log=TRUE))


  par = c(pr,sig)
  k = length(par)
  AIC = 2*k+2*LL
  BIC = k*log(d)+2*LL


  u = pnorm(err,mean=0,sd=sig)


  out=list(coefficients=par, loglik=-LL, aic=AIC,bic=BIC,Matxm=Matxm,
           #clu=clu,Matxm=Matxm,cluster=cluster,
           model="gaussian",thF=thF, u=u, dfM=NULL,disc=NULL)
  class(out) <-  "EstContinuous"
  out

}




