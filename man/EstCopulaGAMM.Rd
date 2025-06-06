% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EstCopulaGAMM.R
\name{EstCopulaGAMM}
\alias{EstCopulaGAMM}
\title{Copula-based estimation of mixed regression models for continuous or discrete response}
\usage{
EstCopulaGAMM(
  y,
  model,
  family = "clayton",
  rot = 0,
  clu,
  xc = NULL,
  xm = NULL,
  start,
  LB,
  UB,
  nq = 25,
  dfC = NULL,
  dfM = NULL,
  offset = NULL,
  prediction = TRUE
)
}
\arguments{
\item{y}{n x 1 vector of response variable (assumed continuous).}

\item{model}{margins: "binomial" or "bernoulli","poisson", "nbinom" (Negative Binomial), "geometric", "multinomial", "gaussian" or "normal", "t", "laplace" , "exponential", "weibull".}

\item{family}{copula family: "gaussian" (normal), "t" , "clayton" ,  "frank" , "fgm", gumbel".}

\item{rot}{rotation: 0 (default), 90, 180 (survival), or 270}

\item{clu}{variable of size n defining the clusters; can be a factor}

\item{xc}{covariates of size n for the estimation of the copula, in addition to the constant; default is NULL.}

\item{xm}{covariates of size n for the estimation of the mean of the margin, in addition to the constant; default is NULL.}

\item{start}{starting point for the estimation; could be the ones associated with a Gaussian-copula model defined by lmer.}

\item{LB}{lower bound for the parameters.}

\item{UB}{upper bound for the parameters.}

\item{nq}{number of nodes and weighted for Gaussian quadrature of the product of conditional copulas; default is 25.}

\item{dfC}{degrees of freedom for a Student margin; default is NULL.}

\item{dfM}{degrees of freedom for a Student margin; default is NULL for non-t distribution,}

\item{offset}{offset (default is NULL)}

\item{prediction}{logical variable for prediction of latent variables V (default is TRUE).}
}
\value{
\item{coefficients}{Estimated parameters}

\item{sd}{Standard deviations of the estimated parameters}

\item{tstat}{T statistics for the estimated parameters}

\item{pval}{P-values of the t statistics for the estimated parameters}

\item{gradient}{Gradient of the log-likelihood}

\item{loglik}{Log-likelihood}

\item{aic}{AIC coefficient}

\item{bic}{BIC coefficient}

\item{cov}{Covariance matrix of the estimations}

\item{grd}{Gradients by clusters}

\item{clu}{Cluster values}

\item{Matxc}{Matrix of covariates defining the copula parameters, including a constant}

\item{Matxm}{Matrix of covariates defining the margin parameters, including a constant}

\item{V}{Estimated value of the latent variable by clusters (if prediction=TRUE)}

\item{cluster}{Unique clusters}

\item{family}{Copula family}

\item{thC0}{Estimated parameters of the copula by observation}

\item{thF}{Estimated parameters of the margins by observation}

\item{rot}{rotation}

\item{dfC}{Degrees of freedom for the Student copula}

\item{model}{Name of the margins}

\item{disc}{Discrete margin number}
}
\description{
This function computes the estimation of a copula-based  2-level hierarchical model.
}
\examples{
data(sim.poisson) #simulated data with Poisson margins
start=c(2,8,3,-1); LB =    c(-3,  3,  -7,  -6);UB=c( 7, 13,   13,   4)
y=sim.poisson$y; clu=sim.poisson$clu;
xc=sim.poisson$xc; xm=sim.poisson$xm
model = "poisson"; family="frank"
out.poisson=EstCopulaGAMM(y,model,family,rot=0,clu,xc,xm,start,LB,UB,nq=31,prediction=TRUE)
}
\references{
Krupskii, Nasri & Remillard (2023). On factor copula-based mixed regression models
}
\author{
Pavel Krupskii, Bouchra R. Nasri and Bruno N. Remillard
}
