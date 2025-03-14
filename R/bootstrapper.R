#' Bootstrapped variance estimates
#'
#' @description Calculates confidence bounds using repeated estimates of likelihood function with bootstrapped draws from dataset.
#' @param LL  A function that specifies the (negative) log likelihood derived from the pdf. Parameters are named to make problem easier. Values are named x and weights are named w.
#' @param x Data values used in parameter estimation.
#' @param w Data weights used in estimation. Default = 1.
#' @param noParams The number of parameters being optimised over in LL (either alpha or alpha and sigma.
#' @param method The algorithm used in the search function. If "Global", uses SANN in optim function. If "Local", uses L-BFGS-B in optim function. Default is local.
#' @param loops Number of loops used to construct the distribution the variances are calculated from. Default is 600.
#' @param alpha Value of alpha parameter in pareto distribution.
#' @param sigma Value of sigma parameter in pareto distribution.
#' @param gamma Value of gamma parameter. If NULL, set to lower threshold or minimum of x.
#' @param nu Value of gamma parameter. If NULL, set to upper threshold or maximum of x.
#' @return Returns matrix of mean, variance, bias, and bias significance for alpha, sigma, gamma, and nu.
#'
#' @import stats
#' @export
#'
#' @examples
#' \dontrun{
#' LL <- function(param, alpha, sigma){
#' alpha <- param[1]
#' sigma <- param[2]
#' logl <-  sum(w(log(1/sigma) - ((1+alpha)/alpha)log((sigma +alpha(x-gamma))/sigma)))
#' return(-logl)
#' }
#' bootstrapper(LL, x=Value, w=Weights, noParams=2, method="Local", loops=600, alpha=1.8, sigma=1.2, gamma=min(Value), nu=max(Value))
#' }

bootstrapper <- function(LL, x, w=1, noParams, method="Local", loops=600, alpha, sigma, gamma=NULL, nu=NULL){
  n <- length(x)
  if (is.null(gamma)){gamma <- min(x)}
  if (is.null(nu)){nu <- max(x)}

  observations <- seq(1,n)
  totalDraw <- sample(observations,n*loops,replace=TRUE)
  sampleDraw <- matrix(totalDraw, n, loops)
  resu <- matrix(NA,loops,4)
  for (i in 1:loops){
    SampleX <- x[sampleDraw[,i]]
    SampleW <- w[sampleDraw[,i]]

    SampleX <- SampleX[order(SampleX, decreasing = TRUE)]
    SampleW <- SampleW[order(SampleX, decreasing = TRUE)]

    BootLace <- LikelihoodOptimisation(LL, x=SampleX, w=SampleW, noParams=noParams, method=method)
    resu[i,] <- unlist(BootLace)
  }
  colnames(resu) <- c("alpha","sigma","gamma","nu")

  colnames(resu) <- c("alpha","sigma","gamma","nu")

  # Prevent infinite values from affecting draws
  resu <- resu[!is.infinite(rowSums(resu, na.rm=TRUE)),]

  resu2 <- colMeans(resu)
  resu2 <- rbind(resu2,diag(var(resu)))
  resu2 <- rbind(resu2, c(mean(resu[,1]-alpha),
                          mean(resu[,2]-sigma),
                          mean(resu[,3]-gamma),
                          mean(resu[,4]-nu)))
  resu2 <- rbind(resu2, c(t.test(resu[,1]-alpha)["p.value"],
                          t.test(resu[,2]-sigma)["p.value"],
                          t.test(resu[,3]-gamma)["p.value"],
                          t.test(resu[,4]-nu)["p.value"]))
  rownames(resu2) <- c("mean","variance","bias","bias significance")
  return(resu2)
}
