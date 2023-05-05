#' @title Link to copula parameter
#'
#' @description Computes the copula parameters given a linear combination of covariates.
#'
#' @param th      vector of linear combinations
#' @param family  copula family: "gaussian" , "t" , "clayton" ,   "claytonR" , "frank" , "gumbel", "gumbelR".
#'
#' @return \item{cpar}{Associated copula parameters}
#' @return \item{hder}{Derivative of link function}
#'
#' @references Krupskii, Nasri & Remillard (2023). On factor copula-based mixed regression models
#' @author Pavel Krupskii and Bruno N. Remillard, January 20, 2023
#' @examples
#' out = linkCop(-1,"gaussian")
#' @export

linkCop = function(th,family="clayton")
{
  minf = function(x){
    max(-0.999,x)
  }
  maxf = function(x){
    min(0.999,x)
  }

expf  = exp(-th);
expf2 = expf*expf
rho   = (1-expf2)/(1+expf2)
rho=sapply(rho,minf)
rho=sapply(rho,maxf)
tau = 1/(1+expf)
tau=sapply(tau,maxf)
switch(family,
       "t" = {
              cpar=rho;

              hder = 1-rho^2
             },
       "gaussian" = {
                     cpar=rho;
                     tau = 2*asin(rho)/pi
                     hder = 1-rho^2
                    },
       "fgm" = {
         cpar=rho;

         hder = 1-rho^2
       },
       "frank" = {
                    cpar = th;
                    hder = 1
                 },
       "clayton" = {
                     cpar=2*tau/(1-tau);
                     hder = 2/expf
                   },

       "gumbel" = { cpar= 1/(1-tau)
                    hder = 1/expf
                  }

          )

out = list(cpar = cpar, hder=hder)
return(out)
}
