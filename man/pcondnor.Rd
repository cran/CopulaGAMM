% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pcond.R
\name{pcondnor}
\alias{pcondnor}
\title{Conditional Gaussian}
\usage{
pcondnor(u, v, cpar)
}
\arguments{
\item{u}{values at which the cdf is evaluated}

\item{v}{value of the conditioning variable in (0,1)}

\item{cpar}{copula parameter}
}
\value{
\item{ccdf}{Conditional cdf}
}
\description{
Conditional Gaussian
}
\examples{
pcondnor(0.5,0.6,0.6)
}
