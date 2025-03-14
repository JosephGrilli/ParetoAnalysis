#' Estimates truncated type 1 pareto parameters
#'
#' @description Returns estimated parameters following the truncated type 1 pareto form.
#' @param x Data values used in parameter estimation.
#' @param w Data weights used in estimation. Default = 1.
#' @param gamma Sets the lower threshold parameter for the type 1 pareto. Default is NULL, setting the threshold to the minimum of x.
#' @param nu Sets the upper threshold parameter for the truncated type 1 pareto. Default is NULL, setting the threshold to the minimum of x.
#' @param OptimSearch The algorithm used in the optim function. If "Global", uses SANN in optim function. If "Local", uses L-BFGS-B in optim function. If "Default", uses the form from Zwijnenburg, Grilli & Engelbrecht (2022).
#' @return Returns parameter estimates for alpha and sigma (and the minimum and maximum for x) for the truncated type 1 pareto.
#' @export
#' @examples
#' \dontrun{
#' TPfunction(Value,Weights,min(Value),max(Value),OptimSearch="Default")
#' }

TPfunction <- function(x,w=1,gamma=NULL,nu=NULL,OptimSearch="Default"){
  if (is.null(gamma)){gamma <- min(x)}
  if (is.null(nu)){nu <- max(x)}
  N <- sum(w)
  # Truncated Pareto Type 1 from Aban, Meerschaert & Panorska (2006), known thresholds
  LL <- function(param, alpha){
    alpha <- param[1]
    logl <-  N*log(alpha) + N*alpha*log(gamma) - N*log(1-(gamma/nu)^alpha) - (alpha)*sum(w*log(x))
    return(-logl)
  }

  if (OptimSearch=="Default"){
    Optimresult <- LikelihoodOptimisation(LL, x=x, w=w,1, method="Local")
  } else {Optimresult <- LikelihoodOptimisation(LL, x=x, w=w, noParams=2, method=OptimSearch)}

  alpha.tp1.hat1 <- Optimresult$alpha
  gamma.tp1.hat1 <- Optimresult$gamma
  nu <- Optimresult$nu


  # Truncated Pareto Type 1 from Aban, Meerschaert & Panorska (2006), unknown thresholds
  LL <- function(param, alpha){
    alpha <- param[1]
    logl <-  sum(w[1:r])*log(alpha) + sum(w[1:r])*alpha*log(x[r+1]) - sum(w[1:r])*log(1-(x[r+1]/nu)^alpha) - (alpha)*sum(w[1:r]*log(x[1:r])) # Truncated
    return(-logl)
  }

  if (OptimSearch=="Default"){
    Optimresult <- LikelihoodOptimisation(LL, x=x, w=w,1, method="Local")
  } else {Optimresult <- LikelihoodOptimisation(LL, x=x, w=w, noParams=2, method=OptimSearch)}

  alpha.tp1.hat2 <- Optimresult$alpha

  # Calculate true threshold from Aban, Meerschaert & Panorska (2006) equation 6
  gamma.tp1.hat2 <- (sum(w[1:r])^(1/alpha.tp1.hat2))*(x[r+1])*((sum(w)-(sum(w)-sum(w[1:r]))*((x[r+1]/nu)^alpha.tp1.hat2))^(-1/alpha.tp1.hat2))

  # Keep new set of parameters if they are finite.
  if (!is.finite(gamma.tp1.hat2)){
    gamma.tp1.hat2 <- gamma.tp1.hat1
    alpha.tp1.hat2 <- alpha.tp1.hat1
  }
  return(c("alpha"=alpha.tp1.hat2,"gamma"=gamma.tp1.hat2,"nu"=nu))
}
