% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/margins.R
\name{multinomcpdf}
\alias{multinomcpdf}
\title{Multinomial with p = 1/(1+exp(-th)) cdf/pdf and ders}
\usage{
multinomcpdf(z, th, x)
}
\arguments{
\item{z}{vector of responses taking values in 1,...,nL: as.number(z) if z is a factor!}

\item{th}{th is a n x (L-1) matrix of parameters, i.e., mpar = a=[a_1,1,...a_1,k2,a_2,1,...a_2,k2,... a_L-1,1... a_L-1,k2], and first level is the baseline.}

\item{x}{matrix of covariates (including the constant)}
}
\value{
\item{out}{Matrix of conditional cdf, derivative with respect to parameter, pdf, }
}
\description{
This function computes the cdf, pdf, and associated derivatives
}
\examples{
x=matrix(c(1,1,-1,-1,0,2),nrow=2)
z = c(1,3)
th = matrix(c(1,2,3,4,5,6),nrow=2)
out = multinomcpdf(z,th,x = x)
}
\author{
Pavel Krupskii and Bruno N. Remillard, January 20, 2022
}
