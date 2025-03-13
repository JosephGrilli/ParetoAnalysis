#' Analytical solution to pareto type 1 parameters
#'
#' @description Returns parameters following the type 1 pareto form, using analytical solution from Aban, Meerschaert & Panorska (2006).
#' @param x Data values used in parameter estimation.
#' @param w Data weights used in estimation. Default = 1.
#' @param r The breakpoint for the Hill estimator. Must be between 1 and n-1.
#' @return Returns parameter estimates for alpha and sigma (and the minimum and maximum for x) for the truncated generalised pareto.
#' @examples
#' Pfunction(Value,w=Weights,r=length(Value)-1)

Pfunction <- function(x,w=1,r=NULL){
  if(is.null(r)){r=length(x)-1}
  if(r<1|r>(length(x)-1)){return(print("r is not specified in range between 1 and length(x)-1"))}
  n <- length(x)
  N <- sum(w)
  gamma.pt1.hat1 <- min(x)
  alpha.pt1.hat1 <- N/(sum(w*log(x))-N*log(gamma.pt1.hat1))

  alpha.pt1.hat2 <- sum(w[1:r])/(sum(w[1:r]*(log(x[1:r]) - log(x[r+1]))))
  gamma.pt1.hat2 <- min(x)
  return(c("alpha"=alpha.pt1.hat2,"gamma"=gamma.pt1.hat2))
}
