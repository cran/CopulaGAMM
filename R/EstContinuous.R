#' @title Copula-based estimation of mixed regression models for continuous response
#'
#' @description This function computes the estimation of a copula-based  2-level hierarchical model.
#'
#' @param y       n x 1 vector of response variable (assumed continuous).
#' @param model   function for margins: "gaussian" (normal), "t" (Student with known df=dfM), laplace" , "exponential", "weibull".
#' @param family  copula family: "gaussian" , "t" , "clayton" ,  "frank" , "fgm", "gumbel".
#' @param rot     rotation: 0 (default), 90, 180 (survival), or 270
#' @param clu     variable of size n defining the clusters; can be a factor
#' @param xc      covariates of size n for the estimation of the copula, in addition to the constant; default is NULL.
#' @param xm      covariates of size n for the estimation of the mean of the margin, in addition to the constant; default is NULL.
#' @param start   starting point for the estimation; could be the ones associated with a Gaussian-copula model defined by lmer.
#' @param LB      lower bound for the parameters.
#' @param UB      upper bound for the parameters.
#' @param nq      number of nodes and weighted for Gaussian quadrature of the product of conditional copulas; default is 25.
#' @param dfM      degrees of freedom for a Student margin; default is 0 for non-t distribution,
#' @param dfC     degrees of freedom for a Student margin; default is 5.
#' @param prediction  logical variable for prediction of latent variables V; default is TRUE.
#' @author Pavel Krupskii, Bouchra R. Nasri and Bruno N. Remillard
#' @return \item{coefficients}{Estimated parameters}
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

EstContinuous =function(y,model,family,rot=0,clu,xc=NULL,xm=NULL,start, LB, UB, nq=31,dfM=NULL,dfC=NULL,prediction=TRUE){

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
  likf=function(par){
    gl=statmod::gauss.quad.prob(nq)
    nl=gl$nodes
    wl=gl$weights
    if(min(par - LB) < 0 || max(par - UB) > 0) return(1e5)

    thC = colSums(par[1:k1]*t(Matxc))  #phi(beta,x_ki)
    thF = cbind(colSums(par[k1+(1:k2)]*t(Matxm)), par[k1+k2+1]) #h(alpha,x_ki)

    cpdf= margins(y,thF,model,dfM)
    # mpdf  = cpdf[,2]  #pdf
    # mpdf1 = cpdf[,5]  #pdf1
    # mpdf2 = cpdf[,6]  #pdf2


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
      wprd = wl*matrixStats::colProds(matrix(tem3,ncol=nq))
      intf = sum(wprd)  +1e-100#c_k
      sl = sum(log(mpdfk))+log(intf) #log(f_k)

      fk[k] = exp(sl) #f_k
      out  = out - sl
      M=matrix(tem4/tem3,nrow=nq,byrow=TRUE)  #\dot c_ki/c_ki
      M1 = matrix(tem5/tem3,nrow=nq,byrow=TRUE)
      grd[k,1:k1] =   - colSums(thCd*colSums(wprd*M))/intf
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
  k = length(par)
  LL = mle$minimum #-log-likelihood
  AIC = 2*k+2*LL
  BIC = k*log(d)+2*LL
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
  out=list(coefficients=mle$estimate, sd = st.dev, tstat=tstat, pval=pval,
           gradient=mle$gradient,loglik=-LL, aic=AIC,bic=BIC,cov=covar,
           grd=grd,clu=clu,Matxc=Matxc,Matxm=Matxm,V=V,cluster=cluster,
           family=family,thC0=thC0,thF=thF, rot=rot,model=model,
           dfC=dfC,dfM=dfM,disc=NULL)
  class(out) <-  "EstContinuous"
  out

}





summary.EstContinous <- function(object, ...) {

  Coef = object$coefficients
  d=  length(Coef)
  source = c("copula","(Intercept)")
  if(d>2)
  {
    for(k in 3:d)
    {
      source=c(source,paste("b",k-2,sep=""))
    }
  }
  if(is.null(object$disc)){source[d]="sigma"}else{if(object$disc==3){source[d]="size"}}
  cat("\nCoefficients:\n")
   out = data.frame(
    Source = source,
    Estimate = Coef,
    Gradient = object$gradient,
    t.value = object$tstat,
    pvalue  =object$pval
  )
  print(out)
  cat("\nAIC:", format(object$aic), "\n")
  cat("\nBIC:", format(object$bic), "\n")

}
