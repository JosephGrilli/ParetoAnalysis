#' Estimate total pareto population
#'
#' @description Function extrapolates the weights in proportion of CDF covered by data (m1= sum(weights)) to estimate the total population (m=m1/F(xmax)) and population of uncovered area (m2).
#' @param x Data values used in parameter estimation.
#' @param w Data weights used in estimation. Default = 1.
#' @param alpha Value of alpha parameter in pareto distribution.
#' @param gamma Value of gamma parameter. If NULL, set to lower threshold or minimum of x.
#' @param sigma Value of sigma parameter in generalised pareto distribution.
#' @param specification Defines distribution CDF used for estimating population. Can be "Type 1" (default)  or "Generalized".
#' @return Returns a table of total population in pareto tail (including origianl values and the newly estimated m2 population), the total population of the tail (popu), and the number of added points (Lm2).
#' @export
#' @examples
#' \dontrun{
#' estimatedParetoPopFunction(x=Value,w=Weights,alpha=1.8,gamma=min(Value),sigma=1.2,specification="Generalized")
#' }

estimatedParetoPopFunction <- function(x,w,alpha,gamma=NULL,sigma=NA,specification="Type 1"){
  if (is.null(gamma)){gamma <- min(x)}
  # Calculate population from tail not contained in survey (m2)
  # m = m(1) + m(2)
  # m(1) = sum(weights)
  # m(1) = m[F(xmax)]
  # m(2) = m - m(1)

  # Calculate population from tail not contained in survey (m2)
  if (specification=="Type 1"){popu <- sum(w)/(1-(gamma/max(x))^alpha)}
  if (specification=="Generalized"){popu <- sum(w)/(1-((sigma + alpha*max(x) - alpha*gamma)/sigma)^(-1/alpha))}

  Lm2 <- floor(popu - sum(w))

  robust.population <- data.frame("Value"=x,"Weights"=w)
  robust.population <- aggregate(Weights ~ Value, robust.population, sum)
  rownames(robust.population) <- NULL
  robust.population <- na.omit(robust.population)
  robust.population <- robust.population[order(robust.population$Value, decreasing = FALSE),]
  robust.population$m1 <- cumsum(robust.population$Weights)

  if (specification=="Type 1"){robust.population$CDF <- 1-(gamma/(robust.population$Value))^alpha}
  if (specification=="Generalized"){robust.population$CDF <- 1-((sigma + alpha*(robust.population$Value) - alpha*gamma)/sigma)^(-1/alpha)}

  robust.population$popu <- robust.population$m1/robust.population$CDF
  robust.population$m2 <- robust.population$popu - robust.population$m1
  n <- nrow(robust.population)
  robust.population <- robust.population[round(n*0.01, digits=0):round(n*1, digits=0),] # Stop extreme values
  robust.population <- na.omit(robust.population)
  df_m1 <- data.frame("Value"=robust.population$Value,"Weights"=robust.population$Weights,
                      "population"=robust.population$m1,"CDF"=robust.population$CDF,"type"="m1", row.names = NULL)
  df_m2 <- data.frame("Value"=robust.population$Value,"Weights"=robust.population$Weights,
                      "population"=robust.population$m2,"CDF"=robust.population$CDF,"type"="m2", row.names = NULL)
  robust.population <- rbind(df_m1,df_m2)
  rownames(robust.population) <- NULL

  result <- structure(
    list(
      table = robust.population,
      population = popu,
      newObs = Lm2
    ),
    row.names = NULL
  )

  return(result)
}
