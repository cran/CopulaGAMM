#' @title Copula-based estimation of mixed regression models for discrete response
#'
#' @description This function computes the estimation of a copula-based  2-level hierarchical model.
#'
#' @param y       n x 1 vector of response variable (assumed continuous).
#' @param model   margins: "binomial" or "bernoulli","poisson", "nbinom" (Negative Binomial), "geometric", "multinomial".
#' @param family  copula family: "gaussian" , "t" , "clayton" ,  "frank" , "fgm", gumbel".
#' @param rot     rotation: 0 (default), 90, 180 (survival), or 270
#' @param clu     variable of size n defining the clusters; can be a factor
#' @param xc      covariates of size n for the estimation of the copula, in addition to the constant; default is NULL.
#' @param xm      covariates of size n for the estimation of the mean of the margin, in addition to the constant; default is NULL.
#' @param start   starting point for the estimation; could be the ones associated with a Gaussian-copula model defined by lmer.
#' @param LB      lower bound for the parameters.
#' @param UB      upper bound for the parameters.
#' @param nq      number of nodes and weighted for Gaussian quadrature of the product of conditional copulas; default is 25.
#' @param dfC     degrees of freedom for a Student margin; default is 0.
#' @param offset  offset (default is NULL)
#' @param prediction  logical variable for prediction of latent variables V (default is TRUE).
#'
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
#' @return \item{cluster}{Unique clusters}
#' @return \item{family}{Copula family}
#' @return \item{thC0}{Estimated parameters of the copula by observation}
#' @return \item{thF}{Estimated parameters of the margins by observation}
#' @return \item{rot}{rotation}
#' @return \item{dfC}{Degrees of freedom for the Student copula}
#' @return \item{model}{Name of the margins}
#' @return \item{disc}{Discrete margin number}
#'
#' @references Krupskii, Nasri & Remillard (2023). On factor copula-based mixed regression models
#' @author Pavel Krupskii, Bouchra R. Nasri and Bruno N. Remillard
#' @import  statmod, matrixStats
#' @examples
#' data(sim.poisson) #simulated data with Poisson margins
#' start=c(2,8,3,-1); LB =    c(-3,  3,  -7,  -6);UB=c( 7, 13,   13,   4)
#' y=sim.poisson$y; clu=sim.poisson$clu;
#' xc=sim.poisson$xc; xm=sim.poisson$xm
#' model = "poisson"; family="frank"
#' out.poisson=EstDiscrete(y,model,family,rot=0,clu,xc,xm,start,LB,UB,nq=31,prediction=TRUE)
#' @export


EstDiscrete=  function(y,model,family, rot = 0, clu,
                          xc=NULL,xm=NULL,start, LB, UB, nq=25,
                          dfC=NULL,offset=NULL, prediction=TRUE)
{
  
  if(model=="bernoulli"){model=="binomial"}
  switch(model,
         "binomial" = { disc= 1  },
         
         "poisson" = {  disc = 2 },
         
         "nbinom" = {   disc = 3 },
         
         "geometric" = { disc = 4  },
         
         "multinomial" = { disc = 5  }
         
  )
  
  
  d = length(y)
  L=2
  z=y
  if(is.character(z)){z=as.factor(z)}
  if(is.factor(z))
  { z=as.numeric(z)
  L = max(z)
  }
  if(is.null(offset)==F)
  { noff=length(offset)
  if(noff!=d){warning("length of offset does not match the total number of variables"); return(NULL) }
  }
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
  offset1=offset
  ind0 = c(1:k2)
  if(is.null(offset1)){ offset1=0}

  thF = matrix(0,ncol= L-1,nrow=d)

  Lk2=(L-1)*k2
  a2 = 2+Lk2            #3:a2   cdf1
  b2 = 1+a2+Lk2         #a2+2:b2  pdf1

  #############  Likelihood function
  likf=function(par){

    gl=statmod::gauss.quad.prob(nq)
    nl=gl$nodes
    wl=gl$weights
    if(min(par - LB) < 0 || max(par - UB) > 0) return(1e100)
    thC =colSums(par[1:k1]*t(Matxc))
    for(j in 1:(L-1))
    {
      ind1 = (j-1)*k2+ind0
      thF[,j] = offset1+colSums(par[k1+ind1]*t(Matxm))
    }
    if(L==2){thF=as.numeric(thF)}
    switch(disc,
           {
             p = 1/(1+exp(-thF)); #P(Y=1)
             p1 = p*(1-p)*Matxm;
             u = 1-p; #P(Y=0)
             u1 = -p1;
           }, #disc = 1
           {
             p  = exp(thF);
             a  = dpois(z,p);
             p1 = p*Matxm;
             u  = ppois(z-1,p);
             v  = u + a;
             u1 = -p1*dpois(z-1,p);
             v1 = -p1*a;
           }, #disc = 2
           {
             size = par[k1+k2+1];
             p  = 1/(1+exp(-thF));
             p1 = p*(1-p)*Matxm;
             thh = cbind(size,p);
             nbcpdf0=nbinomcpdf(z-1,thh);
             nbcpdf1=nbinomcpdf(z,thh);
             u  = nbcpdf0[,1];
             u2 = nbcpdf0[,3];
             u1 = p1*nbcpdf0[,4];
             v  = nbcpdf1[,1];
             v2 = nbcpdf1[,3];
             v1 = p1*nbcpdf1[,4];
           }, #disc = 3
           {
             p = 1/(1+exp(-thF));
             a = dgeom(z,p);
             u = pgeom(z-1,p);
             v = u +a ;
             u1 = z*a*Matxm;
             v1 = a*(1+z)*(1-p)*Matxm;
           }, #disc = 4
           {
             multinom0=multinomcpdf(z-1,thF,Matxm);
             multinom1=multinomcpdf(z  ,thF,Matxm);
             u  = multinom0[,1];
             u1 = multinom0[,3:a2]; #cdf1
             u2 = multinom0[,(a2+2):b2]; #pdf1
             v  = multinom1[,1];
             v1 = multinom1[,3:a2]; #cdf1
             v2 = multinom1[,(a2+2):b2]; #pdf1
           } #disc = 5
    )

    grd = matrix(0,nrow=nclu,ncol=(k1+ Lk2))
    if(disc==3) grd = cbind(grd,0) #  parameter
    # grdsum = rep(0,k1+k2+1)

    out = 0
    fk  = rep(0,nclu)

    for(k in 1:nclu){
      ind = (clu==cluster[k])
      n_k = sum(ind) # number of elements in the cluster


      Matxck = Matxc[ind,]

      Matxmk = Matxm[ind,]


      uu  = rep(u[ind],nq)
      zz  = rep(z[ind],nq)
      uu1 = u1[ind,]
      if(disc>=2){ vv = rep(v[ind],nq); vv1=v1[ind,]; }
      if(disc==3){ uu2=u2[ind]; vv2=v2[ind]; }
      nn  = rep(nl,each=n_k)
      thCk = thC[ind]
      out0=linkCop(thCk,family)
      thC0 = out0$cpar
      thCd = out0$hder*Matxck
      tem = coplik(uu,nn,family,rot,thC0,dfC,TRUE)
      if(disc>=2){ tem0 = coplik(vv,nn,family,rot,thC0,dfC,FALSE)}


      if(disc==1){
        tem1 = zz + (1-2*zz)*tem[,1]
        tem2 = (1-2*zz)*tem[,2]
        tem3 = (1-2*zz)*tem[,3]
        tem1[tem1 < 1e-20] = 1e-20
      }
      if(disc>=2){
        tem1 = tem0[,1] - tem[,1];
        tem2 = tem0[,2] - tem[,2]
        tem3u = tem[,3];
        tem3v = tem0[,3]
        tem1[tem1 < 1e-20] = 1e-20
      }
      
      
      mat1 = matrix(tem1,ncol=nq)
      lmat = colSums(log(mat1))
      lmax = max(lmat)
      wprd = wl*exp(lmat - lmax)

      intf = sum(wprd)
      sl = log(intf) #log(f_k)

      fk[k] = intf #f_k

      ##out  = out - sl + n_k*log(adj)
      out  = out - sl - lmax
      #adj = exp(-lmax/n_k)

      M=matrix(tem2/tem1,nrow=nq,byrow=TRUE)

      grd[k,1:k1] =   - sum(thCd*colSums(wprd*M))/intf

      
      if(disc==1) grd[k,k1+(1:k2)] =  - colSums( wprd*matrix(tem3/tem1,nrow=nq,byrow=TRUE)%*%uu1 )/intf;
      if(disc>=2){
        Mv = matrix(tem3v/tem1,nrow=nq,byrow=TRUE)
        Mu = matrix(tem3u/tem1,nrow=nq,byrow=TRUE)
        grd[k,k1+(1:Lk2)] =  colSums( wprd*Mu%*%uu1 -wprd*Mv%*%vv1)/intf }
      if(disc==3){
        grd[k,k1+k2+1] =   sum(wprd*Mu%*%uu2- wprd*Mv%*%vv2 )/intf;  }
    }

    grdsum = colSums(grd)
    attr(out, "gradient") = grdsum
    attr(out, "grd") = grd
    ##attr(out, "adj") = adj
    out
  }
  mle = nlm(likf,p=start,check.analyticals=FALSE,print.level=0,iterlim=400)

  V=NULL
  thC0=NULL
  par = mle$estimate
  parC = par[1:k1]
  parM = matrix(par[(k1+1):(k1+(L-1)*k2)],ncol=k2,byrow=TRUE)
  size=NULL
  if(disc==3){size = par[k1+k2+1]}
  k = length(par)
  LL = mle$minimum #-log-likelihood
  AIC = 2*k+2*LL
  BIC = k*log(d)+2*LL
  for(j in 1:(L-1))
  {
    ind1 = (j-1)*k2+ind0
    thF[,j] = offset1+colSums(par[k1+ind1]*t(Matxm))
  }


  thC  = colSums(par[1:k1]*t(Matxc))  #phi(beta,x_ki)
  thC0 = linkCop(thC,family)$cpar


  out0=likf(par)
  grd=attributes(out0)$grd
  adj=attributes(out0)$adj


  switch(disc,
         { p  = 1/(1+exp(-thF));
         v  = z +(1-z)*(1-p);
         u = z*(1-p); }, #disc=1

         { p = exp(thF);
         v = ppois(z,p);
         u = ppois(z-1,p); },#disc=2

         { size = par[k1+k2+1];
         p  = 1/(1+exp(-thF));
         u = pnbinom(z-1, size, p);
         v = u+ dnbinom(z, size, p); }, #disc=3

         {  p = 1/(1+exp(-thF));
         u = pgeom(z-1,p);
         v = u+dgeom(z,p);}, #disc=4

         {
           multinom0=multinomcpdf(z-1,thF,Matxm);
           multinom1=multinomcpdf(z  ,thF,Matxm);
           u   = multinom0[,1];
           v   = multinom1[,1];
         } #disc = 5

  )

  if(prediction)
  {
    V=rep(0,nclu)





    for(k in 1:nclu){
      ind = (clu==cluster[k])
      n_k = sum(ind) # number of elements in the cluster
      Matxck = Matxc[ind,]
      nclk = clu[ind]
      thC0k=thC0[ind]

      vv = v[ind]
      uu = u[ind]




      V[k] = MAP.discrete(vv,uu,family,rot,thC0k,dfC,nq)

    }
  }

  covar  = cov(grd)
  st.dev = sqrt(diag(covar)/nclu)
  tstat  = mle$estimate/st.dev
  pval   = 2*pnorm(abs(tstat),lower.tail = FALSE)
  coef=list(copula=parC,margin=parM,size=size)
  out=list(coefficients=coef, sd = st.dev, tstat=tstat, pval=pval,
           gradient=mle$gradient,loglik=-LL, aic=AIC,bic=BIC,cov=covar,
           grd=grd,clu=clu,Matxc=Matxc,Matxm=Matxm, cluster=cluster, V=V,
           family=family,thC0=thC0,thF=thF, dfC=dfC, rot=rot,model=model,disc=disc)

 # class(out) <- "EstDiscrete"
  out


}




