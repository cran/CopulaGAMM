% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coplik.R
\name{ffgmders}
\alias{ffgmders}
\title{Farlie-Gumbel-Morgenstern copula cdf/pdf and ders}
\usage{
ffgmders(u, v, cpar, du = FALSE)
}
\arguments{
\item{u}{vector of values in (0,1)}

\item{v}{conditioning variable in (0,1)}

\item{cpar}{copula parameter in [-1,1]}

\item{du}{logical value (default = FALSE) for the derivative of the copula density with respect to u}
}
\value{
\item{out}{Matrix of conditional cdf, pdf, and derivatives with respect to parameter}
}
\description{
Derivatives  C(u|v), C'_dl(u|v), c(u,v), c'_dl(u,v), c'_u(u,v) for the linking copula
}
\examples{
out = ffgmders(0.3,0.5,2,TRUE)
}
\author{
Pavel Krupskii and Bruno N. Remillard, January 20, 2022
}
