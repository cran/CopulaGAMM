#' @title Inverse conditional cdf
#'
#' @description This function computes the quantile of conditional cdf  C(U|v) for a copula C
#'
#' @param w probability
#' @param v value of the conditioning variable in (0,1)
#' @param family "gaussian" , "t" , "clayton" ,  "fgm",  "frank" , "gumbel", "plackett", "galambos", "huesler-reiss"
#' @param cpar   copula parameter (vector)
#' @param rot    rotation: 0 (default), 90, 180 (survival), or 270
#' @return \item{U}{Conditional quantile}
#' @author Pavel Krupskii and Bruno N. Remillard, January 20, 2022
#' @return \item{U}{Conditional quantile}
#'
#'
#' @examples
#' U = qcond(0.1,0.2,"gaussian",0.87)
#' @export

qcond = function(w,v,family,cpar,rot=0)
{
  switch(family,
         "gaussian" = {

           U= qcondnor(w,v,cpar)
         },

         "t" = {
           U = qcondt(w,v,cpar)
         },

         "joe" = {

           U= qcondjoe(w,v,cpar)
         },

         "clayton" = {

             if(rot==0)
               { U= qcondcla(w,v,cpar) }
             else if (rot==90)
             {U = 1- qcondcla(1-w,v,cpar)}
             else if (rot==180)
             { U= 1- qcondcla(1-w,1-v,cpar)}
             else
             {U = qcondcla(w,1-v,cpar)}
         },


         "frank" = {


           U= qcondfra(w,v,cpar)
         },

         "gumbel" = {
           if(rot==0)
             { U= qcondgum(w,v,cpar) }
           else if (rot==90)
           {U = 1- qcondgum(1-w,v,cpar)}
           else if (rot==180)
           { U= 1- qcondgum(1-w,1-v,cpar)}
           else
           {U = qcondgum(w,1-v,cpar)}
         },


         "plackett" = {
           U = qcondpla(w,v,cpar)
         },

         "huesler-reiss" ={

           if(rot==0)
           { U= qcondhr(w,v,cpar) }
           else if (rot==90)
           {U = 1- qcondhr(1-w,v,cpar)}
           else if (rot==180)
           { U= 1- qcondhr(1-w,1-v,cpar)}
           else
           {U = qcondhr(w,1-v,cpar)}
         },

         "fgm" = {

           if(rot==0)
           { U= qcondfgm(w,v,cpar) }
           else if (rot==90)
           {U = 1- qcondfgm(1-w,v,cpar)}
           else if (rot==180)
           { U= 1- qcondfgm(1-w,1-v,cpar)}
           else
           {U = qcondfgm(w,1-v,cpar)}
         },

         "galambos" = {

           if(rot==0)
           { U= qcondgal(w,v,cpar) }
           else if (rot==90)
           {U = 1- qcondgal(1-w,v,cpar)}
           else if (rot==180)
           { U= 1- qcondgal(1-w,1-v,cpar)}
           else
           {U = qcondgal(w,1-v,cpar)}
         }
  )
  return(U)
}


#' @title Inverse Gumbel
#' @param w probability
#' @param v value of the conditioning variable in (0,1)
#' @param th  copula parameter
#' @return \item{out}{Conditional quantile}
#'
qcondgum = function(w,v,th){
#func1=function(q,th) return(pcond(q,v,"gumbel",cpar=th)) 
func1=function(q,iconv) return(pcond(q,v[!iconv],"gumbel",cpar=th[!iconv])) 
out=invfunc(w,func1,th,tol=1e-8)
out
}


#' @title Inverse Student
#' @param w probability
#' @param v value of the conditioning variable in (0,1)
#' @param th  copula parameter
#' @return \item{out}{Conditional quantile}
#'
#for generating data from t copula
#th[1] - correlation; th[2] - df
qcondt = function(w,v,th){
if(is.vector(th)) th = matrix(th,nrow=1)
rho = th[,1]; nu = th[,2];
t2 = qt(v,nu);
tq = qt(w,nu+1);
snu = sqrt(nu+1);
const = tq*sqrt(1-rho^2)/snu;
out = rho*t2+const*sqrt(nu+t2^2);
out = pt(out,nu);
out
}

#' @title Inverse Gaussian
#' @param w probability
#' @param v value of the conditioning variable in (0,1)
#' @param th  copula parameter (correlation)
#' @return \item{out}{Conditional quantile}
#'
qcondnor = function(w,v,th){
n2 = qnorm(v);
nq = qnorm(w);
out = th*n2+nq*sqrt(1-th^2);
out = pnorm(out);
out
}

#' @title Inverse Joe
#' @param w probability
#' @param v value of the conditioning variable in (0,1)
#' @param th  copula parameter >-1
#' @return \item{out}{Conditional quantile}
#'

qcondjoe = function(w,v,th){
  #func1=function(q,v){ return(pcondjoe(q,v,th)) }
  func1=function(q,iconv){ return(pcondjoe(q,v[!iconv],th[!iconv])) }
  out=invfunc(w,func1,v,tol=1e-8)
  out
}

#' @title Inverse clayton
#' @param w probability
#' @param v value of the conditioning variable in (0,1)
#' @param th  copula parameter
#' @return \item{out}{Conditional quantile}
#'
qcondcla = function(w,v,th){
  a = -th/(1+th)
  #if(th==0){out = w}
  #else{
  out = (1+ v^(-th) *(w^a-1))^(-1/th)
  th0 = (th==0)
  th1 = (th>11)
  out[th0] = w[th0]
  out[th1] = v[th1]
  #out[th>11]=v[th>11]
  #}
  out
}


#' @title Inverse Plackett
#' @param w probability
#' @param v value of the conditioning variable in (0,1)
#' @param th  copula parameter
#' @return \item{out}{Conditional quantile}

#for generating data from Plackett copula
#th - param >-1
qcondpla = function(w,v,th){


  v1<- w*(1-w)
  A <- th+(th-1)^2*v1


  C <- v1* (1+(th-1)*v)^2
  B <- 2*(th-1)*v1* (1-(th+1)*v)-th
  out <- 0.5*(-B + sign(w-0.5)*sqrt(B^2-4*A*C))/A
  out
}


#' @title Inverse Frank
#' @param w probability
#' @param v value of the conditioning variable in (0,1)
#' @param th  copula parameter
#' @return \item{out}{Conditional quantile}
#'
qcondfra = function(w,v,th){
  x2 = exp(-th*v)

  a = exp(-th)
  #if(th==0){out = w}
  #else{
  out = -log(  (x2+w*(a-x2))/(x2+w*(1-x2)) )/th
  th0 = (th==0)
  out[th0] = w[th0]#}
  out
}




#' @title Inverse Galambos
#' @param w probability
#' @param v value of the conditioning variable in (0,1)
#' @param th  copula parameter >0
#' @return \item{out}{Conditional quantile}
#'
qcondgal = function(w,v,th){
  #func1=function(q,v){ return(pcondgal(q,v,th)) }
  func1=function(q,iconv){ return(pcondgal(q,v[!iconv],th[!iconv])) }
  out=invfunc(w,func1,v,tol=1e-8)
  out
}



#' @title Inverse  FGM (B10)
#' @param w probability
#' @param v value of the conditioning variable in (0,1)
#' @param th  copula parameter -1<=th<=1
#' @return \item{out}{Conditional quantile}
#'
qcondfgm=function(w,v,th)
{
  a=th*(1-2*v)
  b = -1-a
  discr=b*b-4*w*a
  2*w/(-b+sqrt(discr)  )
}


#' @title Inverse Huesler-Reiss
#' @param w probability
#' @param v value of the conditioning variable in (0,1)
#' @param th  copula parameter >0
#' @return \item{out}{Conditional quantile}
#'
qcondhr = function(w,v,th){
  #func1=function(q,v){ return(pcondhr(q,v,th)) }
  func1=function(q,iconv){ return(pcondhr(q,v[!iconv],th[!iconv])) }
  out=invfunc(w,func1,v,tol=1e-8)
  out
}

