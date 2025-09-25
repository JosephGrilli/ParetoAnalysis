#' Create samples from estimated top tail using synthetic sampling.
#'
#' @description Generate new observations from a specified distribution using the synthetic sampling approach in Engel et al. (2022).
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
#' @return Returns data frame with values, weights, ids, status as from original data or synthetic, and zipft plot point.
#' @examples
#' synthSample(x=Value,w=Weights,id=hid,m2=m2,popu=popu,alpha=1.8,gamma=100,sigma=NA,Sampling="Inverse",specification="Type 1")

synthSample <- function(x,w=1,id=NA,m2=NULL,popu=NULL,alpha=1.8,gamma=NULL,sigma=NA,Sampling="Inverse",specification="Type 1"){
  if (is.null(gamma)){gamma <- min(x)}
  simu <- data.frame(x,w,id)
  colnames(simu) <- c("Value", "Weights", "id")


  if (is.null(popu)){
    if (specification=="Type 1"){popu <- sum(simu$Weights)/(1-(gamma/max(simu$Value))^alpha)}
    if (specification=="Generalized"){popu <- sum(simu$Weights)/(1-((sigma + alpha*max(simu$Value) - alpha*gamma)/sigma)^(-1/alpha))}
  }
  if (is.null(m2)){m2 <- floor(popu - sum(simu$Weights))}

  # Find values to be attached to these missing observations so that they fill out the distribution
  new.obs <- simu[,c("Value", "Weights", "id")]
  new.obs$type <- "Original"
  if (m2>=1) {
    # Draw observations from CDF (deterministic)
    if (Sampling=="Deterministic"){
      synths <- matrix(data=0, nrow = m2, ncol = 2)

      if (specification=="Type 1") {synths[1,1] <- ((max(simu$Value)^(-alpha))-((gamma^-alpha)*(popu^-1)))^(-1/alpha)}
      if (specification=="Generalized") {synths[1,1] <- gamma + (sigma/alpha)*((((((sigma + alpha*max(simu$Value) - alpha*gamma)/sigma)^(-1/alpha))  - (1/popu))^(-alpha))-1)  }

      if (m2>=2) {
        for (i in 2:m2) {
          if (specification=="Type 1") {synths[i,1] <- ((synths[(i-1),1]^(-alpha))-((gamma^-alpha)*(popu^-1)))^(-1/alpha)}
          if (specification=="Generalized") {synths[i,1] <- gamma + (sigma/alpha)*((((((sigma + alpha*synths[(i-1),1] - alpha*gamma)/sigma)^(-1/alpha))  - (1/popu))^(-alpha))-1)}
        }
      }
    }

    # Draw observations via inverse sampling
    if (Sampling=="Inverse"){
      synths <- matrix(data=0, nrow = m2, ncol = 2)
      #synths[,1] <-runif(m2, min = max(robust.population$CDF[robust.population$Value==max(robust.population$Value)]), max=1)
      if (specification=="Type 1") {m1m <- 1-(gamma/(simu$Value))^alpha}
      if (specification=="Generalized") {m1m <- 1-((sigma + alpha*(simu$Value) - alpha*gamma)/sigma)^(-1/alpha)}
      synths[,1] <-runif(m2, min=m1m, max=1)
      synths <- synths[order(synths[,1]),]
      if (specification=="Type 1") {synths[,1] <- gamma*(1-synths[,1])^(-1/alpha)}
      if (specification=="Generalized") {synths[,1] <- gamma + (sigma/alpha)*( (1-synths[,1])^(-alpha) -1)}
    }

    # Generate values for Zipf plot using new sample
    synths[,2] <- 1
    synths <- as.data.frame(synths)
    colnames(synths) <- c("Value", "Weights")
    synths$id <- NA #New part for ids
    synths$type <- "Synths"
    new.obs <- rbind(new.obs, synths)
  }
  new.obs <- new.obs[order(new.obs$Value, decreasing = TRUE),]
  if (specification=="Type 1") {new.obs$zipf <- (gamma/new.obs$Value)^alpha }
  if (specification=="Generalized") {new.obs$zipf <- ((sigma + alpha*new.obs$Value - alpha*gamma)/sigma)^(-1/alpha)}

  if (!is.null(id)) {
    colnames(new.obs) <- c("Value", "Weights.new", "id", "type", "zipf")
  }
  if (is.null(id)) {
    colnames(new.obs) <- c("Value", "Weights.new", "type", "zipf")
  }
  return(new.obs)
}
