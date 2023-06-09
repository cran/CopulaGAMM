% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/margins.R
\name{nbinom1cpdf}
\alias{nbinom1cpdf}
\title{Negative binomial cdf/pdf and ders}
\usage{
nbinom1cpdf(z, th)
}
\arguments{
\item{z}{vector of responses}

\item{th}{th[,1] is size > 0 and th[,2] is mean > 0; size  does not have to be integer}
}
\value{
\item{out}{Matrix of conditional cdf, derivative with respect to parameter, pdf, }
}
\description{
This function computes the cdf, pdf, and associated derivatives
}
\examples{
out = nbinom1cpdf(0,c(1,0.5))
}
\author{
Pavel Krupskii and Bruno N. Remillard, January 20, 2022
}
