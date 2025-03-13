#' Create samples from estimated top tail using calibrate approach.
#'
#' @description Generate new observations from a specified distribution using the randraw sampling approach in Cantarella, Neri & Ranalli, 2021.
#' @param x Data values used in parameter estimation.
#' @param w Data weights used in estimation. Default = 1.
#' @param id ID values for observations from database to help with linking. Default is NA.
#' @param inputDemoName Names for variables in inputDemo. If NULL (default), set names to column names of inputDemo.
#' @param inputDemo Matrix (nxm) of demographic variables of households used to calibrate weights. Must be same length as x and w.
#' @param alpha Value of alpha parameter in pareto distribution.
#' @param gamma Value of gamma parameter. If NULL, set to lower threshold or minimum of x.
#' @param sigma Value of sigma parameter in generalised pareto distribution.
#' @param specification Defines distribution CDF used for estimating population. Can be "Type 1" (default)  or "Generalized".
#' @param loopmax Maximum number of iterations. Default is 2000.
#' @return Returns data frame of x, w, and id adjusted to match pareto and calibration criteria.
#' @examples
#' calibrateSample(x=Value,w=Weight,id=NA,inputDemoName=c("Gender","Age","Employment"),inputDemo=data[,c("Gender","Age","Employment")],alpha=1.8,gamma=NULL,sigma=NA,specification="Type 1",loopmax=2000)

calibrateSample <- function(x,w=1,id=NA,inputDemoName=NULL,inputDemo=inputDemo,alpha=1.8,gamma=NULL,sigma=NA,specification="Type 1",loopmax=2000){
  if (is.null(gamma)){gamma <- min(x)}
  if (is.null(inputDemoName)){inputDemoName <- colnames(inputDemo)}
  simu <- data.frame(x,w,id,inputDemo)
  colnames(simu) <- c("Value", "Weights", "id",inputDemoName)
  converge <- FALSE

  loop <- 1
  new.obs <- simu[simu$Value>=gamma,c("Value", "Weights", "id")]

  # Set demographic totals
  inputDemo<- new.obs[,inputDemoName]
  aux <- calibVars(inputDemo)
  aux[is.na(aux)] <- 0
  aux.total <- colSums(aux*new.obs$Weights, na.rm = TRUE)
  new.obs$Weights.new <- new.obs$Weights
  new.obs$ECDF <- -1
  if (specification=="Type 1"){
    new.obs$CDF  <- (1 - (gamma/new.obs$Value)^alpha) # Should I scale so CDF and ECDF are same range?/(1 - (gamma/max(new.obs$Value))^alpha)
    new.obs$PDF <- (alpha/new.obs$Value)*((gamma/new.obs$Value)^alpha)
  }
  if (specification=="Generalized") {
    new.obs$CDF  <- 1-((sigma + alpha*(new.obs$Value) - alpha*gamma)/sigma)^(-1/alpha) #
    new.obs$PDF <- (1/sigma)*(((sigma +(alpha*new.obs$Value) -(alpha*gamma))/sigma)^((-1-alpha)/alpha))
  }

  while (converge=="FALSE") {
    # Scale Weights by CDF densities so matches proposed distribution
    # CDF/ECDF * weight: Problem that at x_0, CDF(x_0)=0, so would remove all
    # these weights. However, could do for CDF(x(i+1)), and for x(max),
    # this would be CDF(1) and so not adjusted. However, then also have the
    # problem for multiple entries of same value in complex sample.
    # They would all have the same density, but for CDF(x)=0.1, if we had
    # two observations, both would get scaled to this level.
    # Therefore, while the KS-test may be ok, here, would need to scale the increases.
    # Therefore, need to scale the weight by the CDF(x(i+1)), and by its proportion of
    # the total sum of weights at that level

    # Method to calculate ECDF without for loop
    dist <- aggregate(cbind(Weights.new,CDF) ~ Value, new.obs, sum)
    dist <- dist[order(dist$Value, decreasing=FALSE),]
    dist$ECDF <- (cumsum(dist$Weights.new)/sum(dist$Weights.new))
    new.obs <- new.obs[order(new.obs$Value, decreasing = FALSE),]
    new.obs$ECDF <- merge(dist, new.obs, by.x = "Value", by.y = "Value")[,"ECDF.x"]
    new.obs$scale <- new.obs$CDF/new.obs$ECDF
    # new.obs$scale <- (new.obs$PDF/sum(new.obs$PDF))/(new.obs$Weights.new/sum(new.obs$Weights.new)) # Switched from CDF, else cumulative changes are
    # factored into last observation. However, causes problems that need to aggregate weight, since sum pdf is not 1.

    new.obs$Weights.new <- new.obs$scale * new.obs$Weights.new
    # Retain original demographic totals
    new.obs$g <- calibWeights(aux, new.obs$Weights.new, aux.total)
    new.obs$Weights.new <- new.obs$g * new.obs$Weights.new
    print(paste0(loop,"... gap: ", weighted.mean(x=new.obs$g, w=new.obs$Weights.new)-1))
    converge <- abs((weighted.mean(x=new.obs$g, w=new.obs$Weights.new)-1)) <= 0.0001
    loop <- loop + 1
    if (loop>=loopmax) {
      converge <- TRUE
      print("Loop Max iteration reached.")
    }
  }
  new.obs <- new.obs[order(new.obs$Value, decreasing = TRUE),]
  if (specification=="Type 1"){new.obs$zipf <- (gamma/new.obs$Value)^alpha}
  if (specification=="Generalized") {new.obs$zipf <- ((sigma + alpha*new.obs$Value - alpha*gamma)/sigma)^(-1/alpha) }
  return(new.obs)
}
