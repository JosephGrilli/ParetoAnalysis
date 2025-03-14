#' Estimates generalised pareto parameters
#'
#' @description Returns estimated parameters following the generalised pareto form.
#' @param x Data values used in parameter estimation.
#' @param w Data weights used in estimation. Default = 1.
#' @param gamma Sets the lower threshold parameter for the generalised pareto. Default is NULL, setting the threshold to the minimum of x.
#' @param OptimSearch The algorithm used in the optim function. If "Global", uses SANN in optim function. If "Local", uses L-BFGS-B in optim function. If "Default", uses the form from Zwijnenburg, Grilli & Engelbrecht (2022).
#' @return Returns parameter estimates for alpha and sigma (and the minimum and maximum for x) for the generalised pareto.
#' @export
#' @examples
#' \dontrun{
#' GPfunction(x=Value,w=Weights,gamma=min(Value),OptimSearch="Default")
#' }

GPfunction <- function(x,w=1,gamma=NULL,OptimSearch="Default"){
  if (is.null(gamma)){gamma <- min (x)}
  LL <- function(param, alpha, sigma){
    alpha <- param[1]
    sigma <- param[2]
    logl <-  sum(w*(log(1/sigma) - ((1+alpha)/alpha)*log((sigma +alpha*(x-gamma))/sigma)))
    # Testing bias rank correction (Vermeulen)
    # In Vermeulen, the rank is corrected by -1/N, which means that the top ranked position becomes 0.
    # This is needed to make the last observation the top of the distribution (i.e. CDF=1)
    # This has not been included in the maximum likelihood function, but maybe should be.
    # This would be done in the pdf part of the likeliood, so (x*(f(x)))^w becomes (x*(f(x)-1/N))^w = (x*f(x)-x/N))^w.
    # This means adding -(x*log(n)) to the log likelihood.
    # logl <-  sum(w*(log(1/sigma) - ((1+alpha)/alpha)*log((sigma +alpha*(x-gamma))/sigma))-(x*log(n)))

    # Need to test this by using the observations with the new parameter estimates and calculating the CDF. This should bring the CDF value
    # of the greatest value closer to 1 (using the rankEquation form). It does for my first test, but it is a very marginal improvment.
    # However, this is in simulated data so in real data it might be more important
    return(-logl)
  }
  if (OptimSearch=="Default"){
    Optimresult <- LikelihoodOptimisation(LL, x=x, w=w, noParams=2, method="Local")
  } else {Optimresult <- LikelihoodOptimisation(LL, x=x, w=w, noParams=2, method=OptimSearch)}
  return(c("alpha"=Optimresult$alpha,"sigma"=Optimresult$sigma,"gamma"=Optimresult$gamma,"nu"=Optimresult$nu))
}
