#' Create samples from estimated top tail using randraw sampling.
#'
#' @description Generate new observations from a specified distribution using the randraw sampling approach in Zwijnenburg, Grilli & Engelbrecht, 2022.
#' @param x Data values used in parameter estimation.
#' @param w Data weights used in estimation. Default = 1.
#' @param id ID values for observations from database to help with linking. Default is NA.
#' @param m2 Estimated number of draws required to complete pareto distribution. Output from estimatedParetoPopFunction function. If NULL (default), is calculated using parameter values, specified form and maximum(x), following approach in estimatedParetoPopFunction.
#' @param popu Estimated total number of observations in pareto distribution. Output from estimatedParetoPopFunction function. If NULL (default), is calculated using parameter values, specified form and maximum(x), following approach in estimatedParetoPopFunction.
#' @param alpha Value of alpha parameter in pareto distribution.
#' @param gamma Value of gamma parameter. If NULL, set to lower threshold or minimum of x.
#' @param sigma Value of sigma parameter in generalised pareto distribution.
#' @param Sampling Method of drawing values from the uncovered region. Can be "Deterministic" (method in Zwijnenburg, Grilli & Engelbrecht, 2022) or "Inverse" for inverse sampling.
#' @param specification Defines distribution CDF used for estimating population. Can be "Type 1" (default)  or "Generalized".
#' @return Returns data frame with original values, weights, and the new values calculated for each observation.
#' @export
#' @examples
#' \dontrun{
#' specification(x=Value,w=Weights,id=hid,m2=m2,popu=popu,alpha=1.8,gamma=100,sigma=NA,Sampling="Inverse",specification="Type 1")
#' }

randrawSample <- function(x,w=1,id=NA,m2=NULL,popu=NULL,alpha=1.8,gamma=100,sigma=NA,Sampling="Inverse",specification="Type 1"){
  if (is.null(gamma)){gamma <- min(x)}
  simu <- data.frame(x,w,id)
  colnames(simu) <- c("Value", "Weights", "id")


  if (is.null(popu)){
    if (specification=="Type 1"){popu <- sum(simu$Weights)/(1-(gamma/max(simu$Value))^alpha)}
    if (specification=="Generalized"){popu <- sum(simu$Weights)/(1-((sigma + alpha*max(simu$Value) - alpha*gamma)/sigma)^(-1/alpha))}
  }
  if (is.null(m2)){m2 <- floor(popu - sum(simu$Weights))}
  # Draw random observations from uniform distribution, and randomise across draws
  new.obs <- simu[,c("Value", "Weights", "id")]
  new.obs.2 <- aggregate(Weights ~ Value, new.obs, sum)
  new.obs.2$type <- "Original"
  newsmat <- new.obs.2

  if (m2>=1) {
    synths <- rep(0, m2)

    # Draw observations from CDF (deterministic)
    if (Sampling=="Deterministic"){
      if (specification=="Type 1") {synths[1] <- ((max(simu$Value)^(-alpha))-((gamma^-alpha)*(popu^-1)))^(-1/alpha)}
      if (specification=="Generalized") {synths[1] <- gamma + (sigma/alpha)*((((((sigma + alpha*max(simu$Value) - alpha*gamma)/sigma)^(-1/alpha))  - (1/popu))^(-alpha))-1)}
      if (m2>=2) {
        for (i in 2:m2) {
          if (specification=="Type 1") {synths[i] <- ((synths[(i-1)]^(-alpha))-((gamma^-alpha)*(popu^-1)))^(-1/alpha)}
          if (specification=="Generalized") {synths[i] <- gamma + (sigma/alpha)*((((((sigma + alpha*synths[(i-1)] - alpha*gamma)/sigma)^(-1/alpha))  - (1/popu))^(-alpha))-1)}
        }
      }
    }

    # Draw observations via inverse sampling
    if (Sampling=="Inverse"){
      #synths[,1] <-runif(m2, min = max(robust.population$CDF[robust.population$Value==max(robust.population$Value)]), max=1)
      if (specification=="Type 1") {m1m <- 1-(gamma/(simu$Value))^alpha}
      if (specification=="Generalized") {m1m <- 1-((sigma + alpha*(simu$Value) - alpha*gamma)/sigma)^(-1/alpha)}
      synths[,1] <-runif(m2, min=m1m, max=1)
      synths <- synths[order(synths)]
      if (specification=="Type 1") {synths <- gamma*(1-synths)^(-1/alpha)}
      if (specification=="Generalized") {synths <- gamma + (sigma/alpha)*( (1-synths[,1])^(-alpha) -1)  }
    }

    # Generate values for Zipf plot using new sample
    synths <- data.frame("Value"=synths, "Weights"=1, "type"="Synths")
    newsmat <- rbind(newsmat, synths)
  }

  # Calculate the proportion each observation is of total weight and allocate average of that share to them from new values
  newsmat <- newsmat[order(newsmat$Value, decreasing = TRUE),]
  newsmat$Prop <- (cumsum(newsmat$Weights)/sum(newsmat$Weights))
  newsmat$sumw <- cumsum(newsmat$Weights)

  new.obs.2 <- new.obs.2[order(new.obs.2$Value, decreasing = TRUE),]
  new.obs.2$Prop <- (cumsum(new.obs.2$Weights)/sum(new.obs.2$Weights))
  colnames(newsmat) <- c("value","weight","type","Prop","sumw")

  new.obs.2$nValue <- vapply(new.obs.2$Prop, top_percent_svy, numeric(1), data = newsmat[,c("weight", "value","sumw")])
  new.obs.2$dValue <- c(new.obs.2$nValue[1], diff(new.obs.2$nValue))
  new.obs.2$sValue <- new.obs.2$dValue/new.obs.2$Weights # since is based on share of total weights, need to do here first

  # Match new values to old ones and calculate the per weight value
  # If have multiple observations of same value, need to split the total value between the two, then work out the value per weight

  new.obs$sValue <- new.obs.2$sValue[match(new.obs$Value,new.obs.2$Value)]
  new.obs$sValue[is.na(new.obs$sValue)] <- new.obs$Value[is.na(new.obs$sValue)] # LIS is sometimes missing the matches - this is a contingency measure to stop NAs from reducing values. However, never get this error on my end. There should be a successful one-to-many match, but LIS doesn"t?
  colnames(new.obs) <- c("Value.o","Weights.new","id","Value")

  new.obs <- new.obs[order(new.obs$Value, decreasing = TRUE),]
  return(new.obs)
}
