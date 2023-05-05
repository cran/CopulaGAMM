#' @title Copula pdf
#' @description Evaluates the copula density at given points (u,v)#'
#' @param u vector of values in (0,1)
#' @param v conditioning variable in (0,1)
#' @param family  copula family: "gaussian" ("normal), "t", "clayton",  "frank", "fgm", "galambos", "gumbel", "joe", "huesler-reiss", "plackett".
#' @param cpar copula parameter
#' @param rot    rotation: 0 (default), 90, 180 (survival), or 270
#' @param dfC  degrees of freedom for the Student copula (default is NULL)
#' @return \item{out}{Copula density}
#' @author Pavel Krupskii and Bruno Remillard Mai 1, 2023
#' @return \item{out}{Vector of pdf values}
#' @examples
#' out = dcop(0.3,0.7,"clayton",270,2)
#' @export
#'
dcop = function(u,v,family,rot=0, cpar,dfC=NULL){
  switch(family,
         "gaussian" = {

           out= dbvncop(u,v,cpar)
         },

         "t" = {
           out = dbvtcop(u,v,cpar,dfC)
         },

         "clayton" = {
           if(rot==0){
             out = dmtcj(u,v,cpar)
           }else if(rot==90){
             out = dmtcj(1-u,v,cpar)
             }else if(rot==180){
             out = dmtcj(1-u,1-v,cpar)
            }else
           { out = dmtcj(u,1-v,cpar)}
            },

         "frank" = { out = dfrk(u,v,cpar)
         },

         "gumbel" = {
           if(rot==0){
             out = dgum(u,v,cpar)
           }else if(rot==90){
             out = dgum(1-u,v,cpar)
             }else if(rot==180){
             out = dgum(1-u,1-v,cpar)
             }else
           { out = dgum(u,1-v,cpar)}
           },

           "joe" = {
             if(rot==0){
               out = djoe(u,v,cpar)
             }else if(rot==90){
               out = djoe(1-u,v,cpar)
             }else if(rot==180){
               out = djoe(1-u,1-v,cpar)
             }else
             {
               out = djoe(u,1-v,cpar)}
             },

             "plackett" = {
               if(rot==0){
                 out = dpla(u,v,cpar)
               }else if(rot==90){
                 out = dpla(1-u,v,cpar)
               }else if(rot==180){
                 out = dpla(1-u,1-v,cpar)
               }else
               { out = dpla(u,1-v,cpar)}

               },

               "galambos" = {
                 if(rot==0){
                   out = dgal(u,v,cpar)
                 }else if(rot==90){
                   out = dgal(1-u,v,cpar)
                 }else if(rot==180){
                   out = dgal(1-u,1-v,cpar)
                 }else
                 { out = dgal(u,1-v,cpar)}

                 },
                 "huesler-reiss" = {
                   if(rot==0){
                     out = dhr(u,v,cpar)
                   }else if(rot==90){
                     out = dhr(1-u,v,cpar)
                   }else if(rot==180){
                     out = dhr(1-u,1-v,cpar)
                   }else
                   { out = dhr(u,1-v,cpar)}

                   },

         "fgm" = {
           if(rot==0){
             out = dfgm(u,v,cpar)
           }else if(rot==90){
             out = dfgm(1-u,v,cpar)


           }else if(rot==180){
             out = dfgm(1-u,1-v,cpar)


           }else
           {
             out = dfgm(u,1-v,cpar)

           }
         }


  )
  return(out)

}

#' @title B2 Plackett copula density
#' @description Density at (u,v)
#' @param u vector of values in (0,1)
#' @param v vector of values in (0,1)
#' @param cpar copula parameter > 0
#'
#' @author Pavel Krupskii
#' @return \item{out}{Vector of densities}
#' @examples
#' out = dpla(0.3,0.5,2)
#' @export
#'
dpla=function(u,v,cpar)
{ #if(cpar==1.) return(1.)
  cpar1=cpar-1.;
  tem=1.+cpar1*(u+v); tem1=tem*tem-4.*cpar*cpar1*u*v;
  tem2=sqrt(tem1);
  #pdf=cpar*(1.-u-v+2.*u*v+cpar*(u+v-2.*u*v))/(tem1*tem2);
  tem0=2.*cpar1*u*v; tem3=tem-tem0;
  pdf=cpar*tem3/tem1/tem2;
  #ifelse(cpar==1., 1,pdf)
  pdf
}


#' @title B3 bivariate Frank copula density
#' @description Density at (u,v)
#' @param u vector of values in (0,1)
#' @param v vector of values in (0,1)
#' @param cpar copula parameter, cpar>0 or cpar<0
#'
#' @author Pavel Krupskii
#' @return \item{out}{Vector of densities}
#' @examples
#' out = dfrk(0.3,0.5,2)
#' @export
#'
dfrk=function(u,v,cpar)
{ t1=1.-exp(-cpar);
  tem1=exp(-cpar*u); tem2=exp(-cpar*v);
  pdf=cpar*tem1*tem2*t1;
  tem=t1-(1.-tem1)*(1.-tem2);
  pdf=pdf/(tem*tem);
  pdf
}

#' @title B4 MTCJ copula density,  cpar>0
#' @description Density at (u,v)
#' @param u vector of values in (0,1)
#' @param v vector of values in (0,1)
#' @param cpar copula parameter > 0
#'
#' @author Pavel Krupskii
#' @return \item{out}{Vector of densities}
#' @examples
#' out = dmtcj(0.3,0.5,2)
#' @export
#'
dmtcj=function(u,v,cpar)
{ tem1=u^(-cpar); tem2=v^(-cpar);
  pdf=(tem1+tem2-1)^(-1/cpar-2)*(1+cpar)*tem1*tem2/(u*v)
  pdf
  pdf[v<1e-4]=0
  pdf[u<1e-4]=0
  pdf

}


#' @title B5 Joe copula density
#' @description Density at (u,v)
#' @param u vector of values in (0,1)
#' @param v vector of values in (0,1)
#' @param cpar copula parameter > 1
#'
#' @author Pavel Krupskii
#' @return \item{out}{Vector of densities}
#' @examples
#' out = djoe(0.3,0.5,2)
#' @export
#'
djoe=function(u,v,cpar)
{ f1=1-u; f2=1-v;
  tem1=f1^cpar; tem2=f2^cpar;
  sm=tem1+tem2-tem1*tem2; tem=sm^(1./cpar);
  #cdf=1.-tem;
  pdf=tem*((cpar-1.)*tem1*tem2+tem1*tem1*tem2+tem1*tem2*tem2-tem1*tem1*tem2*tem2)
  pdf=pdf/(sm*sm);
  pdf=pdf/(f1*f2);
  pdf
}

#' @title B6 Gumbel copula density, cpar>1
#' @description Density at (u,v)
#' @param u vector of values in (0,1)
#' @param v vector of values in (0,1)
#' @param cpar copula parameter > 0
#'
#' @author Pavel Krupskii
#' @return \item{out}{Vector of densities}
#' @examples
#' out = dgum(0.3,0.5,2)
#' @export
#'
dgum=function(u,v,cpar)
{ l1= -log(u); l2= -log(v);
  tem1=l1^cpar; tem2=l2^cpar; sm=tem1+tem2; tem=sm^(1./cpar);
  cdf=exp(-tem);
  pdf=cdf*tem*tem1*tem2*(tem+cpar-1.);
  pdf=pdf/(sm*sm*l1*l2*u*v);
  pdf
}

#' @title B7 Galambos copula density, cpar>0
#' @description Density at (u,v)
#' @param u vector of values in (0,1)
#' @param v vector of values in (0,1)
#' @param cpar copula parameter > 0
#'
#' @author Pavel Krupskii
#' @return \item{out}{Vector of densities}
#' @examples
#' out = dgal(0.3,0.5,2)
#' @export
#'
dgal=function(u,v,cpar)
{ cpar1= -cpar;
  l1= -log(u); l2= -log(v);
  tem1=l1^(cpar1-1.); tem2=l2^(cpar1-1.);
  sm=tem1*l1+tem2*l2; tem=sm^(1./cpar1-1.);
  lcdf=-(l1+l2-tem*sm);  # cdf=exp(lcdf);
  deriv12= 1.+ (tem/sm)*tem1*tem2*( 1.+cpar+tem*sm) -tem*(tem1+tem2);
  pdf=lcdf+l1+l2;
  pdf=exp(pdf)*deriv12;
  pdf
}

#' @title B8 Huesler-Reiss copula density, cpar>0
#' @description Density at (u,v)
#' @param u vector of values in (0,1)
#' @param v vector of values in (0,1)
#' @param cpar copula parameter > 0
#'
#' @author Pavel Krupskii
#' @return \item{out}{Vector of densities}
#' @examples
#' out = dhr(0.3,0.5,2)
#' @export
#'
dhr=function(u,v,cpar)
{ x=-log(u); y=-log(v);
  z=x/y; cpar1=1./cpar; lz=log(z);
  tem1=cpar1+.5*cpar*lz; tem2=cpar1-.5*cpar*lz;
  p1=pnorm(tem1); p2=pnorm(tem2);
  cdf=exp(-x*p1-y*p2);
  pdf=cdf*(p1*p2+.5*cpar*dnorm(tem1)/y)/(u*v);
  pdf
}





#' @title Normal copula density
#' @description Density at (u,v)
#' @param u vector of values in (0,1)
#' @param v vector of values in (0,1)
#' @param cpar copula parameter, -1< cpar<1
#'
#' @author Pavel Krupskii
#' @return \item{out}{Vector of densities}
#' @examples
#' out = dbvncop(0.3,0.5,-0.5)
#' @export
#'
dbvncop=function(u,v,cpar)
{ x1=qnorm(u); x2=qnorm(v)
  qf=x1^2+x2^2-2*cpar*x1*x2
  qf=qf/(1-cpar^2)
  con=sqrt(1-cpar^2)*(2*pi)
  pdf=exp(-.5*qf)/con
  pdf=pdf/(dnorm(x1)*dnorm(x2))
  pdf
}

#' @title Student copula density
#' @description Density at (u,v)
#' @param u vector of values in (0,1)
#' @param v vector of values in (0,1)
#' @param cpar copula parameter, -1< cpar<1
#' @param dfC degrees of freedom
#'
#' @author Pavel Krupskii
#' @return \item{out}{Vector of densities}
#' @examples
#' out = dbvtcop(0.3,0.5,-0.7,25)
#' @export
#'
dbvtcop=function(u,v,cpar,dfC)
{
  rh=cpar; nu=dfC
  con=exp(lgamma((nu+2.)/2.)-lgamma(nu/2.))/(pi*nu);
  con=con/sqrt(1.-rh*rh);
  ex=-(nu+2.)/2.;
  xt1=qt(u,nu); den1=dt(xt1,nu);
  xt2=qt(v,nu); den2=dt(xt2,nu);
  r11=1./(1-rh*rh); r12=-rh*r11;
  qf=xt1*xt1*r11+xt2*xt2*r11+2.*xt1*xt2*r12;
  tem=1.+qf/nu;
  pdf=con*(tem^ex)/(den1*den2);
  #lpdf=log(con)+ex*log(tem)-log(den1)-log(den2);
  pdf
}

#' @title Farlie-Gumbel-Morgenstern copula density, -1<= cpar<=
#' @description Density at (u,v)
#' @param u vector of values in (0,1)
#' @param v vector of values in (0,1)
#' @param cpar copula parameter > 0
#'
#' @author Pavel Krupskii
#' @return \item{out}{Vector of densities}
#' @examples
#' out = dfgm(0.3,0.5,0.2)
#' @export
#'
dfgm=function(u,v,cpar)
{
  pdf = 1+cpar*(1-2*u)*(1-2*v)
}

#' @title Normal density
#' @description Density at (x1,x2)
#' @param x1 vector of values
#' @param x2 vector of values
#' @param rh correlation parameter, -1< rh <1
#'
#' @author Pavel Krupskii
#' @return \item{out}{Vector of densities}
#' @examples
#' out = dbvn(0.3,0.5,-0.6)
#'
#'
dbvn=function(x1,x2,rh)
{ #qf=x1^2+x2^2-2*rh*x1*x2
  qf=x1^2+x2^2-2*rh*x1*x2
  qf=qf/(1-rh^2)
  con=sqrt(1-rh^2)*(2*pi)
  pdf=exp(-.5*qf)/con
  pdf
}

#' @title Normal density (version 2)
#' @description Density at (x1,x2)
#' @param x1 vector of values
#' @param x2 vector of values
#' @param rh correlation parameter, -1< rh <1
#'
#' @author Pavel Krupskii
#' @return \item{out}{Vector of densities}
#' @examples
#' out = dbvn2(0.3,0.5,-0.4)
#'
#'
dbvn2=function(x1,x2,rh)
{ qf=x1^2+x2^2-2*rh*x1*x2
qf=qf/(1-rh^2)
con=sqrt(1-rh^2)*(2*pi)
pdf=exp(-.5*qf)/con
pdf
}
