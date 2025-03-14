#' Estimates truncated generalised pareto parameters
#'
#' @description Returns estimated parameters following the truncated generalised pareto form.
#' @param x Data values used in parameter estimation.
#' @param w Data weights used in estimation. Default = 1.
#' @param gamma Sets the lower threshold parameter for the generalised pareto. Default is NULL, setting the threshold to the minimum of x.
#' @param nu Sets the upper threshold parameter for the truncated generalised pareto. Default is NULL, setting the threshold to the minimum of x.
#' @param OptimSearch The algorithm used in the optim function. If "Global", uses SANN in optim function. If "Local", uses L-BFGS-B in optim function. If "Default", uses the form from Zwijnenburg, Grilli & Engelbrecht (2022).
#' @return Returns parameter estimates for alpha and sigma (and the minimum and maximum for x) for the truncated generalised pareto.
#' @export
#' @examples
#' \dontrun{
#' TGPfunction(Value,Weights,min(Value),max(Value),OptimSearch="Default")
#' }

TGPfunction <- function(x,w=1,gamma=NULL,nu=NULL,OptimSearch="Default"){
  if (is.null(gamma)){gamma <- min(x)}
  if (is.null(nu)){nu <- max(x)}
  LL <- function(param, alpha, sigma){
    alpha <- param[1]
    sigma <- param[2]
    logl <-  sum(w*(log(1/sigma) - ((1+alpha)/alpha)*log((sigma +alpha*(x-gamma))/sigma) - log(1-(((sigma + alpha*(nu-gamma))/sigma)^(-1/alpha)))))
    return(-logl)
  }

  if (OptimSearch=="Default"){
    Optimresult <- LikelihoodOptimisation(LL, x=x, w=w, noParams=2, method="Global")
  } else {Optimresult <- LikelihoodOptimisation(LL, x=x, w=w, noParams=2, method=OptimSearch)}
  return(c("alpha"=Optimresult$alpha,"sigma"=Optimresult$sigma,"gamma"=Optimresult$gamma,"nu"=Optimresult$nu))
}
