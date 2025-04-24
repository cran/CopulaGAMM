#' @title Inverse function
#'
#' @description This function is used to get the inverse of a monotonic function on (0,1), depending on parameters, and using the bisection method
#'
#' @param q   Function value (can be a vector if func() supports a vector argument)
#' @param func Function of one argument to be inverted
#' @param th   Function parameters
#' @param lb  Lower bound for the possible values
#' @param ub  Upper bound for the possible values
#' @param tol Tolerance for the inversion
#' @param nbreak Maximum number of iterations (default is 40)
#' @author Pavel Krupskii
#' @return \item{out}{Inverse values}
#'
#'
#' @export




invfunc = function(q,func,th,lb=1e-12,ub=1-1e-12,tol=1e-8,nbreak=40){
ibreak=0
lq = length(q);
if(length(th)==1) th = rep(th[1],lq)
#lth= length(th);
v0 = rep(lb,lq);
v1 = rep(ub,lq);
out = rep(0, lq);
iconv = rep(FALSE, lq);
func0 = function(v,q,th) return(func(v,th)-q)
#f0 = func0(v0,q,th); f1 = func0(v1,q,th);
f0 = func0(v0,q,iconv); f1 = func0(v1,q,iconv);
l0 = (f0 > 0 & f1 > 0);
l1 = (f0 < 0 & f1 < 0);
iconv[l0 | l1] = TRUE;
tol0 = rep(1,lq); tol0[iconv] = 0;
while(prod(iconv) < 0.5 & ibreak < nbreak) {
  ibreak=ibreak+1
  v0a = v0[!iconv]; v1a = v1[!iconv];
  v2a = (v0a+v1a)/2; qa = q[!iconv];
  f0a = f0[!iconv]; f1a = f1[!iconv];
  tha=th[!iconv];
  #f2a = func0(v2a,qa,tha);
  f2a = func0(v2a,qa,iconv);
  i1 = (f0a*f2a < 0); i2 = (f1a*f2a < 0);
  v1a[i1] = v2a[i1]; v1[!iconv] = v1a;
  v0a[i2] = v2a[i2]; v0[!iconv] = v0a;

  f0[!iconv][i1] = f0a[i1];
  f0[!iconv][i2] = f2a[i2];
  f1[!iconv][i1] = f2a[i1];
  f1[!iconv][i2] = f1a[i2];

  tol0[!iconv] = abs(f2a);
  iconv = (abs(tol0) < tol)

}

out=(v0+v1)/2
out[l0] = lb; out[l1] = ub
out
}



