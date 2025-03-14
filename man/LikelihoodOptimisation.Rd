% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LikelihoodOptimisation.R
\name{LikelihoodOptimisation}
\alias{LikelihoodOptimisation}
\title{(Log) Likelihood optimisation function for Pareto}
\usage{
LikelihoodOptimisation(LL, x, w = 1, noParams, method = c("Global", "Local"))
}
\arguments{
\item{LL}{A function that specifies the (negative) log likelihood derived from the pdf. Parameters are named to make problem easier. Values are named x and weights are named w.}

\item{x}{Data values used in parameter estimation.}

\item{w}{Data weights used in estimation. Default = 1.}

\item{noParams}{The number of parameters being optimised over.}

\item{method}{The algorithm used in the search function. If "Global", uses SANN in optim function. If "Local", uses L-BFGS-B in optim function.}
}
\value{
Returns parameter estimates for alpha and sigma (and the minimum and maximum for x)
}
\description{
This function finds parameters that maximise the likelihood function using pareto pdfs and survey data.
}
\examples{
\dontrun{
LL <- function(param, alpha, sigma){
alpha <- param[1]
sigma <- param[2]
logl <-  sum(w(log(1/sigma) - ((1+alpha)/alpha)log((sigma +alpha(x-gamma))/sigma)))
return(-logl)
}
LikelihoodOptimisation(LL, x=x, w=w,1, method="Local")
}
}
