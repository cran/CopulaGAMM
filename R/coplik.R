#' @title Copula cdf/pdf and ders
#' @description Derivatives C(u|v), C'_dl(u|v), c(u,v), c'_dl(u,v), c'_u(u,v) for the linking copula
#' @param u vector of values in (0,1)
#' @param v conditioning variable in (0,1)
#' @param family  copula family: "gaussian", "t", "clayton",  "frank", "fgm", "gumbel", "joe", "plackett".
#' @param rot    rotation: 0 (default), 90, 180 (survival), or 270
#' @param cpar copula parameter
#' @param dfC degrees of freedom for the Student copula (default is NULL)
#' @param du logical value (default = FALSE) for the derivative of the copula density with respect to u
#'
#' @author Pavel Krupskii and Bruno N. Remillard, January 20, 2022
#' @return \item{out}{Matrix of conditional cdf, pdf, and derivatives with respect to parameters}
#' @examples
#' out = coplik(0.3,0.5,"clayton",cpar=2,du=TRUE)
#' @export
#'
coplik = function(u,v,family,rot=0,cpar,dfC=NULL,du = FALSE){
  switch(family,
         "gaussian" = {
           out= fnorders(u,v,cpar,du)
         },

         "t" = {
          out = ftders(u,v,cpar,dfC,du)
         },

         "clayton" = {
           if(rot==0){
             out = fmtcjders(u,v,cpar,du)
           }else if(rot==90){
              out = fmtcjders(1-u,v,cpar,du)
              out[,1]=1-out[,1]
              out[,c(2,5)] = -out[,c(2,5)]

           }else if(rot==180){
               out = fmtcjders(1-u,1-v,cpar,du)
               out[,1]=1-out[,1]
               out[,c(2,5)] = -out[,c(2,5)]

             }else
             {
               out = fmtcjders(u,1-v,cpar,du)

             }

         },
         "plackett" = {
           if(rot==0){
             out = fpladers(u,v,cpar,du)
           }else if(rot==90){
             out = fpladers(1-u,v,cpar,du)
             out[,1]=1-out[,1]
             out[,c(2,5)] = -out[,c(2,5)]

           }else if(rot==180){
             out = fpladers(1-u,1-v,cpar,du)
             out[,1]=1-out[,1]
             out[,c(2,5)] = -out[,c(2,5)]

           }else
           {
             out = fpladers(u,1-v,cpar,du)

           }

         },

         "joe" = {
           if(rot==0){
             out = fjoeders(u,v,cpar,du)
           }else if(rot==90){
             out = fjoeders(1-u,v,cpar,du)
             out[,1]=1-out[,1]
             out[,c(2,5)] = -out[,c(2,5)]

           }else if(rot==180){
             out = fjoeders(1-u,1-v,cpar,du)
             out[,1]=1-out[,1]
             out[,c(2,5)] = -out[,c(2,5)]

           }else
           {
             out = fjoeders(u,1-v,cpar,du)

           }

         },


         "frank" = {

           out = ffrkders(u,v,cpar,du)
         },

         "gumbel" = {
           if(rot==0){
             out = fgumders(u,v,cpar,du)
           }else if(rot==90){
             out = fgumders(1-u,v,cpar,du)
             out[,1]=1-out[,1]
             out[,c(2,5)] = -out[,c(2,5)]

           }else if(rot==180){
             out = fgumders(1-u,1-v,cpar,du)
             out[,1]=1-out[,1]
             out[,c(2,5)] = -out[,c(2,5)]

           }else
           {
             out = fgumders(u,1-v,cpar,du)

           }
         },


         "fgm" = {
           if(rot==0){
             out = ffgmders(u,v,cpar,du)
           }else if(rot==90){
             out = ffgmders(1-u,v,cpar,du)
             out[,1]=1-out[,1]
             out[,c(2,5)] = -out[,c(2,5)]

           }else if(rot==180){
             out = ffgmders(1-u,1-v,cpar,du)
             out[,1]=1-out[,1]
             out[,c(2,5)] = -out[,c(2,5)]

           }else
           {
             out = ffgmders(u,1-v,cpar,du)

           }
         }


           )
  return(out)

}

#' @title Plackett copula cdf/pdf and ders
#' @description Derivatives C(u|v), C'_dl(u|v), c(u,v), c'_dl(u,v), c'_u(u,v) for the linking copula
#' @param u vector of values in (0,1)
#' @param v conditioning variable in (0,1)
#' @param cpar copula parameter > 0
#' @param du logical value (default = FALSE) for the derivative of the copula density with respect to u
#'
#' @author Pavel Krupskii and Bruno N. Remillard, January 20, 2022
#' @return \item{out}{Matrix of conditional cdf, pdf, and derivatives with respect to parameter}
#' @examples
#' out = fpladers(0.3,0.5,2,TRUE)
#' @export
#'
fpladers=function(u,v,cpar,du = FALSE){
  eta = cpar-1
  uv = u*v
  uv2 = u+v-2*uv
  tem1 = sqrt((1+eta*(u+v))^2-4*eta*cpar*uv)
  tem3 = tem1^3
  tem5 = tem1^5
  temd = 2*(1+eta*(u+v))*(u+v)-4*uv*(eta+cpar)
  num1 = eta*v+1-(eta+2)*u
  num2 = 1+eta*uv2
  ccond = .5 - .5*num1/tem1
  ccond1 = -.5*(v-u)/tem1 + .25*num1*temd/tem3
  pdf = cpar*num2/tem3
  pdf1 = (1+(2*cpar-1)*uv2)/tem3 - 1.5*cpar*num2*temd/tem5
  pdfu = 0
  if(du==T){
    pdfu = cpar*eta*(1-2*v)/tem3-1.5*cpar*num2*(2*eta*(1+eta*(u+v))-4*eta*cpar*v)/tem5
  }
  out = cbind(ccond,pdf,ccond1,pdf1,pdfu)
  out
}


#' @title Joe copula cdf/pdf and ders
#' @description Derivatives C(u|v), C'_dl(u|v), c(u,v), c'_dl(u,v), c'_u(u,v) for the linking copula
#' @param u vector of values in (0,1)
#' @param v conditioning variable in (0,1)
#' @param cpar copula parameter > 1
#' @param du logical value (default = FALSE) for the derivative of the copula density with respect to u
#'
#' @author Pavel Krupskii and Bruno N. Remillard, January 20, 2022
#' @return \item{out}{Matrix of conditional cdf, pdf, and derivatives with respect to parameter}
#' @examples
#' out = fjoeders(0.3,0.5,2,TRUE)
#' @export
#'
fjoeders=function(u,v,cpar,du = FALSE){
  ud = (1-u)^cpar
  vd = (1-v)^cpar
  lu = log(1-u)
  lv = log(1-v)
  ctem1 = 1+ud/vd-ud
  ctem2 = ud+vd-ud*vd
  num1 = ud*lu+vd*lv-ud*vd*(lu+lv)
  ccond = ctem1^(-1+1/cpar)*(1-ud)
  pdf = ctem2^(-2+1/cpar)*ud*vd*(cpar-1+ctem2)/(1-u)/(1-v)
  ccond1 = -log(ctem1)/cpar^2 - ud*lu/(1-ud)+(-1+1/cpar)*((lu-lv)*ud/vd-ud*lu)/ctem1
  ccond1 = ccond1*ccond
  pdf1 = -log(ctem2)/cpar^2 + (-2+1/cpar)*num1/ctem2 + (1+num1)/(cpar-1+ctem2) + lu + lv
  pdf1 = pdf1*pdf
  pdfu=0
  if(du==T){
    tem2 = cpar*ud*(vd-1)/(1-u)
    pdfu = (-2+1/cpar)*tem2/ctem2+(1-cpar)/(1-u)+tem2/(cpar-1+ctem2)
    pdfu = pdf*pdfu
  }
  out = cbind(ccond,pdf,ccond1,pdf1,pdfu)
  out[u==0,3] = 0
  out[v==1,] = 0
  out[u==1,2:5] = 0
  out
}

#' @title Clayton copula cdf/pdf and ders
#' @description Derivatives C(u|v), C'_dl(u|v), c(u,v), c'_dl(u,v), c'_u(u,v) for the linking copula
#' @param u vector of values in (0,1)
#' @param v conditioning variable in (0,1)
#' @param cpar copula parameter > 0
#' @param du logical value (default = FALSE) for the derivative of the copula density with respect to u
#'
#' @author Pavel Krupskii and Bruno N. Remillard, January 20, 2022
#' @return \item{out}{Matrix of conditional cdf, pdf, and derivatives with respect to parameter}
#' @examples
#' out = fmtcjders(0.3,0.5,2,TRUE)
#' @export
#'
fmtcjders=function(u,v,cpar,du = FALSE){
    ud=u^(-cpar);
    vd=v^(-cpar);
    cpar2=1/cpar/cpar;
    lu = log(u);
    lv = log(v);
    ud1=ud/u;
    vd1=vd/v;
    uvd=ud+vd-1;
    luvd=log(uvd);
    cdf=(ud+vd-1)^(-1/cpar);
    ccond=vd1*cdf/uvd;
    pdf=(1+cpar)*ud1*ccond/uvd;
    tem1=-(lu*ud+lv*vd)/uvd;
    lccond1=-lv+cpar2*luvd-(1/cpar+1)*tem1;
    ccond1=ccond*lccond1;
    lpdf1=1/(cpar+1)-lu-lv+cpar2*luvd-(1/cpar+2)*tem1;
    pdf1=pdf*lpdf1;

    pdfu=0;
    if(du==TRUE) pdfu = pdf*(-cpar-1+(2*cpar+1)*ud/(ud+vd-1))/u

    out = cbind(ccond,ccond1,pdf,pdf1,pdfu)
    out[u==0,1:5] = 0
    out[v==0,2:5] = 0
    out[v==0,1] = 1
    out
}


#' @title Gumbel copula cdf/pdf and ders
#' @description Derivatives C(u|v), C'_dl(u|v), c(u,v), c'_dl(u,v), c'_u(u,v) for the linking copula
#' @param u vector of values in (0,1)
#' @param v conditioning variable in (0,1)
#' @param cpar copula parameter > 1
#' @param du logical value (default = FALSE) for the derivative of the copula density with respect to u
#'
#' @author Pavel Krupskii and Bruno N. Remillard, January 20, 2022
#' @return \item{out}{Matrix of conditional cdf, pdf, and derivatives with respect to parameter}
#' @examples
#' out = fgumders(0.3,0.5,2,TRUE)
#' @export
#'
fgumders=function(u,v,cpar,du = FALSE){
  ul=-log(u);
  uul=u*ul;
  vl=-log(v);
  ull=log(ul);
  vll=log(vl);
  ut=ul^cpar;
  vt=vl^cpar;
  uvt=ut+vt;
  luvt=log(uvt);
  uvd=uvt^(1/cpar);
  cdf=exp(-uvd);
  ul1=ut/uul;
  vl1=vt/v/vl;
  uvdt=uvd/uvt;
  l1v=vl1*uvdt;
  l1u=ul1*uvdt;
  l2uv=(cpar-1)*l1u*vl1/uvt;
  ccond=cdf*l1v;

  l1uv=l1u*l1v;
  l1cdf=cdf*l1uv;
  l2cdf=cdf*l2uv;
  pdf=l1cdf+l2cdf;

  pdfu=0;
  if(du==TRUE){
    cf1=-l1u/(uvd+cpar-1);
    cf2=-(1-2*cpar)*ul1/uvt;
    cf3=-(cpar-1)/uul;
    pdfu=pdf*(l1u+cf1+cf2+cf3-1/u);
  }

  tem1=(ut*ull+vt*vll)/uvt;
  tem2=-luvt/cpar/cpar;
  uvd1=uvd*(tem2+tem1/cpar);
  tem12=tem2+(1/cpar-1)*tem1
  ccond1=ccond*(-uvd1+vll+tem12);
  pdf1=-pdf*uvd1+l1cdf*(ull+vll+2*tem12)+l2cdf*(1/(cpar-1)+ull+vll+tem12-tem1);
  out=cbind(ccond, ccond1, pdf, pdf1, pdfu)
  out[u==0,] = 0
  out[u==1,2:5] = 0
  out[u==1,1] = 1
  out[v==0,2:5] = 0
  out[v==0,1] = 1
  out[v==1,] = 0
  out
}



#' @title Frank copula cdf/pdf and ders
#' @description Derivatives  C(u|v), C'_dl(u|v), c(u,v), c'_dl(u,v), c'_u(u,v) for the linking copula
#' @param u vector of values in (0,1)
#' @param v conditioning variable in (0,1)
#' @param cpar copula parameter
#' @param du logical value (default = FALSE) for the derivative of the copula density with respect to u
#'
#' @author Pavel Krupskii and Bruno N. Remillard, January 20, 2022
#' @return \item{out}{Matrix of conditional cdf, pdf, and derivatives with respect to parameter}
#' @examples
#' out = ffrkders(0.3,0.5,2,TRUE)
#' @export
#'
ffrkders=function(u,v,cpar,du = FALSE){
  ind0 = (cpar==0)
  ind1 = (cpar!=0)
  n=length(u)
  ccond=rep(0,n)
  pdf = ccond
  ccond1=pdf
  pdf1=pdf
  pdfu=pdf
  if(sum(ind0)>0)
  { ccond[ind0]=u[ind0]
  pdf[ind0] = 1
  pdfu[ind0]=0
  ccond1[ind0]= u[ind0] *(1-u[ind0])*(1-v[ind0])
  pdf1[ind0] = u[ind0]*v[ind0] + (1-u[ind0])*(1-v[ind0])
  }else{

    u = u[ind1]
    v = v[ind1]
    cpar = cpar[ind1]
    de=exp(-cpar)
    ue=exp(-cpar*u)
    ve=exp(-cpar*v)
    uve=ue*ve
    uue=u*ue
    vve=v*ve
    uvee=uve*(u+v)
    den=(1-de)-(1-ue)*(1-ve)
    den1 = de-uue-vve+uvee
    den2 = den*den
    ccond[ind1]=ve*(1-ue)/den
    ccond1[ind1]=(-vve+uvee)/den - ccond[ind1]*den1/den
    pdf[ind1] = cpar*(1-de)*uve/den2
    num1= (1-de)*uve + de*cpar*uve - cpar*(1-de)*uvee
    pdf1[ind1]= num1/den2-2*pdf[ind1]*den1/den

    pdfu[ind1] = 0;
    if(du==TRUE) pdfu[ind1] = cpar*pdf[ind1]*(2*ue*(1-ve)/den-1)

  }
  cbind(ccond, ccond1, pdf, pdf1, pdfu)
}

#' @title Farlie-Gumbel-Morgenstern copula cdf/pdf and ders
#' @description Derivatives  C(u|v), C'_dl(u|v), c(u,v), c'_dl(u,v), c'_u(u,v) for the linking copula
#' @param u vector of values in (0,1)
#' @param v conditioning variable in (0,1)
#' @param cpar copula parameter in [-1,1]
#' @param du logical value (default = FALSE) for the derivative of the copula density with respect to u
#'
#' @author Pavel Krupskii and Bruno N. Remillard, January 20, 2022
#' @return \item{out}{Matrix of conditional cdf, pdf, and derivatives with respect to parameter}
#' @examples
#' out = ffgmders(0.3,0.5,2,TRUE)
#' @export
#'
ffgmders=function(u,v,cpar,du = FALSE){
  ccond=  u*(1+cpar*(1-u)*(1-2*v))
  ccond1= u*(1-u)*(1-2*v)
  pdf = 1+cpar*(1-2*u)*(1-2*v)
  pdf1 = (1-2*u)*(1-2*v)

  pdfu=0;
  if(du==TRUE) pdfu = -2*cpar*(1-2*v)

  cbind(ccond, ccond1, pdf, pdf1, pdfu)
}

#' @title Farlie-Gumbel-Morgenstern copula cdf/pdf and ders
#' @description Derivatives  C(u|v), C'_dl(u|v), c(u,v), c'_dl(u,v), c'_u(u,v) for the linking copula
#' @param u vector of values in (0,1)
#' @param v conditioning variable in (0,1)
#' @param cpar copula parameter in (-1,1)
#' @param du logical value (default = FALSE) for the derivative of the copula density with respect to u
#'
#' @author Pavel Krupskii and Bruno N. Remillard, January 20, 2022
#' @return \item{out}{Matrix of conditional cdf, pdf, and derivatives with respect to parameter}
#' @examples
#' out = fnorders(0.3,0.5,0.6,TRUE)
#' @export
#'

fnorders=function(u,v,cpar,du = FALSE){
  rho = cpar
  Cpi = sqrt(2*pi)
  qu = qnorm(u)
  qv = qnorm(v)
  qu2= qu*qu
  qv2= qv*qv
  quv= qu*qv
  rho2 = rho*rho
  srho2 = sqrt(1-rho2)
  srho3 = srho2*(1-rho2)
  zuv = (qu-rho*qv)/srho2
  ccond=pnorm(zuv)
  euv = exp(-.5*zuv*zuv)/Cpi
  ccond1=euv*(rho*qu-qv)/srho3
  pdf = exp(.5*(-rho2*(qu2+qv2)+2*rho*quv)/(1-rho2))/srho2
  pdf1 = pdf*((quv*(1+rho2)-rho*(qu2+qv2))/(1-rho2)/(1-rho2) + rho/(1-rho2))


  pdfu= 0;
  if(du){ pdfu = pdf*Cpi*rho*(qv-rho*qu)*exp(.5*qu2)/(1-rho2)}


  out=cbind(ccond, ccond1, pdf, pdf1, pdfu)
  out[u==0,2:5] = 0
  out[u==1,2:5] = 1
  out
}

#' @title Student copula cdf/pdf and ders
#' @description Derivatives  C(u|v), C'_dl(u|v), c(u,v), c'_dl(u,v), c'_u(u,v) for the linking copula
#' @param u vector of values in (0,1)
#' @param v conditioning variable in (0,1)
#' @param cpar copula parameter in (-1,1)
#' @param dfC degrees of freedom
#' @param du logical value (default = FALSE) for the derivative of the copula density with respect to u
#'
#' @author Pavel Krupskii and Bruno N. Remillard, January 20, 2022
#' @return \item{out}{Matrix of conditional cdf, pdf, and derivatives with respect to parameter}
#' @examples
#' out = ftdersP(0.3,0.5,2,25,TRUE)
#' @export
#'
ftdersP=function(u,v,cpar,dfC,du = FALSE){
  rho = cpar
  tu = qt(u,dfC)
  tv = qt(v,dfC)
  tu2 = tu*tu
  tv2 = tv*tv
  rho2 = (1-rho*rho)
  den = rho2*(dfC+tv2)/(dfC+1)
  sden = sqrt(den)
  tuv = (tu - rho*tv)/sden
  ccond = pt(tuv,dfC+1)
  ccond1 = dt(tuv,dfC+1)*(rho*tu - tv)/rho2/sden
  cf = gamma(dfC/2+1)*gamma(dfC/2)/(gamma(dfC/2+1/2))^2/sqrt(rho2)
  tem0 = (tu2+tv2-2*rho*tu*tv)/dfC/rho2
  tem1 = (2*rho*(tu2+tv2)-2*(2-rho2)*tu*tv)/dfC/rho2/rho2
  pdf  = cf*(1+tem0)^(-dfC/2-1)*((1+tu2/dfC)*(1+tv2/dfC))^(dfC/2+1/2)
  pdf1 = pdf*(rho/rho2-(dfC/2+1)*tem1/(tem0+1))
  pdfu= 0;
  if(du==TRUE){
    dtu = dt(tu,dfC)
    tem1u = 2*(tu-rho*tv)/dfC/rho2
    pdfu = tu*(dfC+1)/(dfC+tu2) - (dfC/2+1)*tem1u/(1+tem0)
    pdfu = pdfu*pdf/dtu
  }
  out=cbind(ccond, ccond1, pdf, pdf1, pdfu)
  cf1 = sqrt(rho2/(dfC+1))
  out[u==0,2:4] = 0
  out[u==1,2:4] = 0
  out[u==0,5] = Inf
  out[u==1,5] = -Inf
  out[v==0,3:5] = 0
  out[v==1,3:5] = 0
  out[v==0, 1] = 1-pt(-rho/cf1,dfC+1)
  out[v==1, 1] = pt(-rho/cf1,dfC+1)
  out[v==0, 2] = dt(rho/cf1,dfC+1)/rho2/cf1
  out[v==1, 2] = -dt(rho/cf1,dfC+1)/rho2/cf1
  out
}

#' @title Student copula cdf/pdf and ders
#' @description Derivatives  C(u|v), C'_dl(u|v), c(u,v), c'_dl(u,v), c'_u(u,v) for the linking copula
#' @param u vector of values in (0,1)
#' @param v conditioning variable in (0,1)
#' @param cpar copula parameter in (-1,1)
#' @param nu degrees of freedom >0
#' @param du logical value (default = FALSE) for the derivative of the copula density with respect to u
#'
#' @author Pavel Krupskii and Bruno N. Remillard, January 20, 2022
#' @return \item{out}{Matrix of conditional cdf, pdf, and derivatives with respect to parameter}
#' @examples
#' out = ftders(0.3,0.5,2,25)
#' @export
#'
ftders=function(u,v,cpar,nu,du = FALSE){

  x= qt(u,nu);
  y= qt(v,nu);
  rho=cpar
  Omega = 1-rho^2;
  mu = y*rho;
  z = (nu+y^2)/(nu+1);
  sig = sqrt(Omega*z);
  w = (x-mu)/sig;
  ccond = pt(w,nu+1);
  wu = 1/sig/dt(x,nu)
  pdf = dt(w,nu+1)*wu;
  dd = (y - rho*w*sig/Omega)
  w1 = -dd/sig;

  ccond1 = dt(w,nu+1)*w1 ;
  pdf1 = pdf*( rho/Omega -  w*(nu+2)/(nu+1) *w1/(1+w^2/(nu+1))  )

  pdfu= 0;
  if(du) {
    pdfu = pdf*wu*(x*sig*(nu+1)/(nu+x^2) - w*(nu+2)/(nu+1+w^2)) }

  out=cbind(ccond, ccond1, pdf, pdf1, pdfu)
  out[u==0,2:5] = 0
  out[u==1,2:5] = 1
  out
}


#' @title Estimation of the parameter of  a bivariate copula (Clayton, Frank, Gumbel)
#' @description Computes the MLE estimation for a bivariate copula using gradient. The likelihood is likelihood is c(u,v;theta)
#' @param u vector of values in (0,1)
#' @param v vector of values in (0,1)
#' @param fcopders ffrkders, fgumders or fmtcjders
#' @param start starting value for the parameter (default =2)
#' @param LB   lower bound for the parameter (default is 1.01)
#' @param UB   upper bound for the parameter (default is 7)
#' @author Pavel Krupskii
#' @return \item{mle}{List of outputs from nlm function}
#' @examples
#' set.seed(2)
#' v = runif(250)
#' w = runif(250)
#' u = 1/sqrt(1+(w^(-2/3)-1)/v^2) # Clayton copula with parameter 2 (tau=0.5)
#' out = mlecop(u,v,fmtcjders)
#' @export
#'
mlecop=function(u,v,fcopders,start=2, LB=1.01, UB=7){
likf=function(par){
  if(par < LB || par > UB) return(1e5)
  tem = fcopders(u,v,par)
  out = -sum(log(tem[,3]))
  grd = -sum(tem[,4]/tem[,3])
  attr(out, "gradient") = grd
  out
  }
mle = nlm(likf,p=start,check.analyticals=FALSE)
return(mle)
}

#' @title Estimation of the parameter of a bivariate copula (Clayton, Frank, Gumbel) when the first observation is 0 or 1
#' @description Computes the MLE estimation for a bivariate copula using gradient. The likelihood is likelihood is C(1-p|v;theta) if y=0 and 1-C(1-p|v;theta) if y=1
#' @param y vector of binary values 0 or 1
#' @param v vector of values in (0,1)
#' @param fcopders ffrkders, fgumders or fmtcjders
#' @param start starting value for the parameter (default =2)
#' @param LB   lower bound for the parameter (default is 1.01)
#' @param UB   upper bound for the parameter (default is 7)
#' @author Pavel Krupskii
#' @return \item{mle}{List of outputs from nlm function}
#' @examples
#' set.seed(2)
#' v = runif(250)
#' w = runif(250)
#' u = 1/sqrt(1+(w^(-2/3)-1)/v^2) #Clayton with parameter 2
#' y = as.numeric(u>0.6) # if one takes (u<4), one obtains a rotation of the Clayton!
#' out = mlecop.disc(y,v,fmtcjders)
#' @export
#'
mlecop.disc=function(y,v,fcopders,start=2, LB=1.01, UB=7){
  pr = mean(y)
  U=cbind(1-pr,v)
likf=function(par){
  if(par < LB || par > UB) return(1e5)
  tem = fcopders(U[,1],U[,2],par)
  tem1 = matrix(tem[y==1,],ncol=5)
  tem11 = 1-tem1[,1]+1-1e-20
  tem0 = matrix(tem[y==0,],ncol=5)
  tem00 = tem0[,1]+1e-20
  out0 = -sum(log(tem00))
  grd0 = -sum(tem0[,2]/tem00)
  out1 = -sum(log(tem11))
  grd1 = sum(tem1[,2]/tem11)
  out  = out0+out1
  attr(out, "gradient") = grd0+grd1
  out
  }
mle = nlm(likf,p=start,check.analyticals=F,print.level=0)
return(mle)
}

