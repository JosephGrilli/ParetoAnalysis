#' (Log) Likelihood optimisation function for Pareto
#'
#' @description This function finds parameters that maximise the likelihood function using pareto pdfs and survey data.
#' @param LL A function that specifies the (negative) log likelihood derived from the pdf. Parameters are named to make problem easier. Values are named x and weights are named w. Example for Generalised Pareto: LL <- function(param, alpha, sigma){
#' alpha <- param[1]
#' sigma <- param[2]
#' logl <-  sum(w(log(1/sigma) - ((1+alpha)/alpha)log((sigma +alpha(x-gamma))/sigma)))
#' return(-logl)
#' }
#' @param x Data values used in parameter estimation.
#' @param w Data weights used in estimation. Default = 1.
#' @param noParams The number of parameters being optimised over.
#' @param method The algorithm used in the search function. If "Global", uses SANN in optim function. If "Local", uses L-BFGS-B in optim function.
#' @return Returns parameter estimates for alpha and sigma (and the minimum and maximum for x)
#' @examples
#' LikelihoodOptimisation(LL, x=x, w=w,1, method="Local")

# Call likelihood functions (bounded gradient local maximum, and global simulated annealing)
LikelihoodOptimisation <- function(LL, x, w=1, noParams, method = c("Global","Local")){

  # mu[[alpha, gamma]]
  parameters <- c(rep(3, noParams))

  if (method=="Local"){
    result <- optim(parameters, LL, method = "L-BFGS-B", lower = 0.0001)[["par"]]
  }

  # Global maximum solution  for smaller samples to get better results - for larger samples this takes too long
  if (method=="Global"){
    if (length(x)<=40000) {
      print("Sample size less than 40,000. Applying SANN Global Optimisation.")
    } else {print("Sample size more than 40,000. SANN Global Optimisation algorithm may perform slowly at this size - consider using OptimSearch=\"Local\"")}
    result <- optim(parameters, LL, method = "SANN")[["par"]]
  }

  output <- list("alpha"=result[1], "sigma"=result[2], "gamma"=min(x), "nu"=max(x))
  return(output)
}
