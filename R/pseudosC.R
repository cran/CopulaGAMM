#' @title Estimation cdf, left-continuous cdf, and pseudo-observations
#'
#' @description This function estimates the empirical cdf, its left limit, and pseudo-observations for a univatiate vector y
#'
#' @param y univariate data
#'
#' @author Bruno N. Remillard, January 20, 2022
#' @return \item{Fn}{Emprirical cdf}
#' @return \item{Fm}{Left-contniuous cdf}
#' @return \item{U}{Pseudo-obsevations}
#'
#'
#' @examples
#' y = rpois(100,2)
#' out=pseudosC(y)
#'
#' @export
#'



pseudosC <- function(y)
{
  n = length(y)

  out0 = .C("stats",
            as.double(y),
            as.integer(n),
            Fn = double(n),
            Fm = double(n),
            PACKAGE = "CopulaGAMM"
  )

  out0$U = n*out0$Fn/(n+1)
  return(out0)


}
