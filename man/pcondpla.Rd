% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pcond.R
\name{pcondpla}
\alias{pcondpla}
\title{Conditional Plackett (B2)}
\usage{
pcondpla(u, v, cpar)
}
\arguments{
\item{u}{values at which the cdf is evaluated}

\item{v}{value of the conditioning variable in (0,1)}

\item{cpar}{copula parameter >1}
}
\value{
\item{ccdf}{Conditional cdf}
}
\description{
Conditional Plackett (B2)
}
\examples{
pcondpla(0.5,0.6,2)
}
