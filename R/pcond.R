#' @title Conditional cdf
#' @description This function computes the conditional cdf C(U|V) for a copula C
#' @param U values at which the cdf is evaluated
#' @param V value of the conditioning variable in (0,1)
#' @param family "gaussian" , "t" , "clayton" ,  "joe", "frank" , "fgm", gumbel", "plackett", "galambos", "huesler-reiss"
#' @param rot    rotation: 0 (default), 90, 180 (survival), or 270
#' @param cpar   copula parameter (vector)
#' @param dfC    degrees of freedom of the Student copula (default is NULL)
#'
#' @author Pavel Krupskii and Bruno N. Remillard, January 20, 2022
#' @return \item{p}{Conditional cdf}
#' @examples
#' p = pcond(0.1,0.2,"clayton",rot=270,cpar=0.87)
#' @export

pcond = function(U,V,family,rot=0,cpar,dfC=NULL){
  switch(family,
         "gaussian" = {

           p= pcondnor(U,V,cpar)
         },

         "t" = {
           p = pcondt(U,V,cpar,dfC)
         },

         "joe" = {
           if(rot==0) { p= pcondjoe(U,V,cpar) }
           else if (rot==90)
           {p = 1-pcondjoe(1-U,V,cpar)}
           else if (rot==180)
           { p= 1- pcondjoe(1-U,1-V,cpar)}
           else
           {p = pcondjoe(U,1-V,cpar)}
           },

         "frank" = {
            p= pcondfrk(U,V,cpar)
         },

         "clayton" = {
           if(rot==0) { p= pcondcla(U,V,cpar) }
           else if (rot==90)
           {p = 1-pcondcla(1-U,V,cpar)}
           else if (rot==180)
           { p= 1- pcondcla(1-U,1-V,cpar)}
           else
           {p = pcondcla(U,1-V,cpar)}

         },

         "frank" = {

           p= pcondfrk(U,V,cpar)
         },

         "gumbel" = {
           if(rot==0) { p=pcondgum(U,V,cpar) }
           else if(rot==90){p=1-pcondgum(1-U,V,cpar)}
           else if(rot==180){ p=1-pcondgum(1-U,1-V,cpar) }
           else {p = pcondgum(U,1-V,cpar)}
         },

         "plackett" = {
           p = pcondpla(U,V,cpar)
         },
         "galambos" = {
           if(rot==0) { p=pcondgal(U,V,cpar) }
           else if(rot==90){p=1-pcondgal(1-U,V,cpar)}
           else if(rot==180){ p=1-pcondgal(1-U,1-V,cpar) }
           else {p = pcondgal(U,1-V,cpar)}
         },
         "huesler-reiss" = {
           if(rot==0) { p=pcondhr(U,V,cpar) }
           else if(rot==90){p=1-pcondhr(1-U,V,cpar)}
           else if(rot==180){ p=1-pcondhr(1-U,1-V,cpar) }
           else {p = pcondhr(U,1-V,cpar)}
         }
  )
  p[U==0]=0
  p[U==1]=1
return(p)

}


#' @title Conditional Gaussian
#' @param u values at which the cdf is evaluated
#' @param v value of the conditioning variable in (0,1)
#' @param cpar  copula parameter
#'
#' @return \item{ccdf}{Conditional cdf}
#' @examples
#' pcondnor(0.5,0.6,0.6)
#' @export
#'
pcondnor<-function(u,v,cpar)
{
  rho = cpar
ccdf<-pnorm((qnorm(u)-rho*qnorm(v))/sqrt(1-rho^2))
ccdf[ v <= 0 | u <= 0 | v >= 1] <- 0
ccdf[ u == 1 ] <- 1
ccdf
}

#' @title Conditional Student
#' @description Conditional Student is Y2|Y1=y1 ~ t(nu+1,location=rho*y1, sigma(y1)), where here sigma^2 = (1-rho^2)(nu+y1^2)/(nu+1)
#' @param u values at which the cdf is evaluated
#' @param v value of the conditioning variable in (0,1)
#' @param cpar  copula parameter
#' @param dfC   degrees of freedom
#'
#' @return \item{ccdf}{Conditional cdf}
#' @examples
#' pcondt(0.5,0.6,0.6,15)
#' @export
#'
  pcondt=function(u,v,cpar,dfC)
  {
  rh=cpar; nu=dfC
  v[v==0]<-1.e-5  # temporary fix
  u[u==0]<-1.e-6
  v[v==1]<-1-1.e-5  # temporary fix
  u[u==1]<-1-1.e-6
  y1=qt(v,nu);
  y2=qt(u,nu);
  mv=rh*y1;
  s2=(1.-rh*rh)*(nu+y1*y1)/(nu+1.);
  ccdf=pt((y2-mv)/sqrt(s2),nu+1.);
  ccdf
  }

  #' @title Conditional Plackett (B2)
  #' @param u values at which the cdf is evaluated
  #' @param v value of the conditioning variable in (0,1)
  #' @param cpar  copula parameter >1
  #'
  #' @return \item{ccdf}{Conditional cdf}
  #' @examples
  #' pcondpla(0.5,0.6,2)
  #' @export
  #'
  pcondpla=function(u,v,cpar)
  { #if(cpar==1.) return(u)
    eta=cpar-1.;
    tem=1.+eta*(v+u); tem1=tem*tem-4.*cpar*eta*v*u;
    tem2=sqrt(tem1);
    ccdf=(eta*v+1.-(eta+2.)*u)/tem2;
    ccdf=.5*(1.-ccdf);
    #ifelse(cpar==1., u, ccdf)
    ccdf
  }

  #' @title Conditional Frank (B3)
  #' @param u values at which the cdf is evaluated
  #' @param v value of the conditioning variable in (0,1)
  #' @param cpar  copula parameter
  #'
  #' @return \item{ccdf}{Conditional cdf}
  #' @examples
  #' pcondfrk(0.5,0.6,2)
  #' @export
  #'
  pcondfrk=function(u,v,cpar)
  { #if(cpar==0.) return(u)
    cpar1=1.-exp(-cpar);
    tem=1.-exp(-cpar*v);
    ccdf=(1.-tem)/(cpar1/(1.-exp(-cpar*u))-tem);
    ccdf
  }

  #' @title Conditional Clayton
  #' @param u values at which the cdf is evaluated
  #' @param v value of the conditioning variable in (0,1)
  #' @param cpar  copula parameter
  #'
  #' @return \item{ccdf}{Conditional cdf}
  #' @examples
  #' pcondcla(0.5,0.6,2)
  #' @export
  #'
  pcondcla=function(u,v,cpar)
  {
    n = length(u)
    ccdf=rep(0,n)
    ind0=(cpar==0)
    ind1=(cpar>0)
    ccdf[ind0] = u[ind0]

    vv=v[ind1]
    a = vv^cpar
    uu=u[ind1]
    ccdf[ind1] = (1+a*(uu^(-cpar)-1))^(-1-1/cpar)

    ccdf[v==0]=1
    ccdf[u==0]=0
    ccdf

  }




  #' @title Conditional Joe (B5)
  #' @param u values at which the cdf is evaluated
  #' @param v value of the conditioning variable in (0,1)
  #' @param cpar  copula parameter
  #'
  #' @return \item{ccdf}{Conditional cdf}
  #' @examples
  #' pcondjoe(0.5,0.6,2)
  #' @export
  #'
   pcondjoe=function(u,v,cpar)
  {
  temu=(1.-u)^cpar;
  temv=(1.-v)^cpar;
  ccdf=1.+temu/temv-temu;
  ccdf=ccdf^(-1.+1./cpar);
  ccdf=ccdf*(1.-temu);
  ccdf
  }


  #' @title Conditional Gumbel (B6)
  #' @param u values at which the cdf is evaluated
  #' @param v value of the conditioning variable in (0,1)
  #' @param cpar  copula parameter >1
  #
  #' @return \item{ccdf}{Conditional cdf}
  #' @examples
  #' pcondgum(0.5,0.6,2)
  #' @export
  #'
  pcondgum=function(u,v,cpar)
  { v[v==0]<-1.e-7  # temporary fix (to be consistent with pgum
    v[v==1]<-1-1.e-7
  u[u==0]<-1.e-7
  #u[u==0]<-1.e-6
  x= -log(v); y= -log(u);
  tem1=x^cpar; tem2=y^cpar; svm=tem1+tem2; tem=svm^(1./cpar);
  ccdf=exp(-tem);
  ccdf=ccdf*(1+tem2/tem1)^(-1.+1./cpar);
  ccdf=ccdf/v;
  ccdf
  }


  #' @title Conditional Galambos (B7)
  #' @param u values at which the cdf is evaluated
  #' @param v value of the conditioning variable in (0,1)
  #' @param cpar  copula parameter
  #'
  #' @return \item{ccdf}{Conditional cdf}
  #' @examples
  #' pcondgal(0.5,0.6,2)
  #' @export
  #'
  pcondgal=function(u, v, cpar)
  { x= -log(v); y= -log(u);
  tem1=x^(-cpar); tem2=y^(-cpar); sm=tem1+tem2; tem=sm^(-1./cpar);
  ccdf=exp(-(x+y-tem));
  ccdf=ccdf*(1.-(1+tem2/tem1)^(-1.-1./cpar));
  ccdf=ccdf/v;
  ccdf
  }

  #' @title Conditional Huesler-Reiss (B8)
  #' @param u values at which the cdf is evaluated
  #' @param v value of the conditioning variable in (0,1)
  #' @param cpar  copula parameter >0
  #'
  #' @return \item{ccdf}{Conditional cdf}
  #' @examples
  #' pcondhr(0.5,0.6,2)
  #' @export
  #'
  pcondhr=function(u, v, cpar)
  { x= -log(v); y= -log(u);
  z=x/y; cpar1=1./cpar; lz=log(z);
  tem1=cpar1+.5*cpar*lz; tem2=cpar1-.5*cpar*lz;
  p1=pnorm(tem1); p2=pnorm(tem2);
  lcdf=-x*p1-y*p2;  cdf=exp(lcdf);
  ccdf=cdf*p1/v
  ccdf
  }

  #' @title Conditional FGM (B10)
  #' @param u probability
  #' @param v value of the conditioning variable in (0,1)
  #' @param cpar  copula parameter -1<=cpar<=1
  #' @return \item{ccdf}{Conditional cdf}
  #' @examples
  #' pcondfgm(0.5,0.6,0.9)
  #' @export
  #'

  pcondfgm=function(u, v, cpar)
  { tem=1+cpar*(1-2*v)*(1-u)
  u*tem
  }
