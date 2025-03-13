#' Hill parameter estimator
#'
#' @description Loops over dataset to find Hill parameter r that maximises the fit to the Pareto Type 1 distribution.
#' @param x Data values used in parameter estimation.
#' @param w Data weights used in estimation. Default = 1.
#' @return Length of observation of x that should be used to combat bias in pareto Type 1 estimate.
#' @examples
#' HillLoopfunction(Value,Weights)

HillLoopfunction <- function(x,w){
  simu <- data.frame("Value"=x,"Weights"=w)
  df.1 <- c(2,0)
  for (r in round(n*0.05, digits=0):round(n*0.95, digits=0)) {
    alpha.pt1.hat2 <- sum(simu$Weights[1:r])/(sum(simu$Weights[1:r]*(log(simu$Value[1:r]) - log(simu$Value[r+1]))))

    CDF <- 1 - (simu$Value[r+1]/simu$Value[1:r])^alpha.pt1.hat2

    dist <- aggregate(Weights ~ Value, simu[1:r,], sum)
    dist <- dist[order(dist$Value,decreasing=FALSE),]
    dist$ECDF <- (cumsum(dist$Weights)/sum(dist$Weights))

    ECDF <-dist$ECDF[match(simu$Value[1:r],dist$Value)]

    df <- c(-(cor(c(CDF), c(ECDF))^2), r)
    if(!is.na(df[1])){
      if (df[1]<=df.1[1]) {df.1 <-df}
    }
  }
  r <- df.1[2]
  return(r)
}
