#' @title Copula-based estimation of mixed regression models for continuous or discrete response
#'
#' @description This function computes the estimation of a copula-based  2-level hierarchical model.
#'
#' @param y       n x 1 vector of response variable (assumed continuous).
#' @param model   margins: "binomial" or "bernoulli","poisson", "nbinom" (Negative Binomial), "geometric", "multinomial", "gaussian" or "normal", "t", "laplace" , "exponential", "weibull".
#' @param family  copula family: "gaussian" (normal), "t" , "clayton" ,  "frank" , "fgm", gumbel".
#' @param rot     rotation: 0 (default), 90, 180 (survival), or 270
#' @param clu     variable of size n defining the clusters; can be a factor
#' @param xc      covariates of size n for the estimation of the copula, in addition to the constant; default is NULL.
#' @param xm      covariates of size n for the estimation of the mean of the margin, in addition to the constant; default is NULL.
#' @param start   starting point for the estimation; could be the ones associated with a Gaussian-copula model defined by lmer.
#' @param LB      lower bound for the parameters.
#' @param UB      upper bound for the parameters.
#' @param nq      number of nodes and weighted for Gaussian quadrature of the product of conditional copulas; default is 25.
#' @param dfC     degrees of freedom for a Student margin; default is NULL.
#' @param dfM     degrees of freedom for a Student margin; default is NULL for non-t distribution,
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
#' out.poisson=EstCopulaGAMM(y,model,family,rot=0,clu,xc,xm,start,LB,UB,nq=31,prediction=TRUE)
#' @export

EstCopulaGAMM=  function(y,model,family="clayton", rot = 0, clu,
                       xc=NULL,xm=NULL,start, LB, UB, nq=25,
                       dfC=NULL,dfM=NULL,offset=NULL, prediction=TRUE)
{

dismod=c("binomial", "bernoulli","poisson", "nbinom", "geometric", "multinomial")
contmod = c("gaussian","normal", "t", "laplace" , "exponential", "weibull")

if(family=="normal"){family="gaussian"}
out=NULL
if(model %in% dismod){
  out=EstDiscrete(y,model,family, rot, clu, xc,xm,start, LB, UB, nq,dfC,offset, prediction)
}else
{
  if(model %in% contmod){
  out = EstContinuous(y,model,family, rot, clu, xc,xm,start, LB, UB, nq,dfC,dfM, prediction)
  }
}
class(out) <- "CopulaGAMM"
return(out)
}
#function(
