#' Calculates confidence bounds for pseudo-maximum likelihood
#'
#' @description Calculates confidence bounds for pseudo-maximum likelihood using sandwich variance-covariance matrix (Geyer (2013): https://www.stat.umn.edu/geyer/5601/notes/sand.pdf).
#' @param x Data values used in parameter estimation.
#' @param w Data weights used in estimation. Default = 1.
#' @param EXPRESS expression function giving the log likelihood function. This is differentiated with respect to ParamNames. Example for generalised pareto: expression(w(log(1/sigma) - ((1+alpha)/alpha)log((sigma +alpha(x-100.0129))/sigma)))
#' @param ParamNames The names of the parameters used in the expression. Example for generalised pareto: c("alpha", "sigma")
#' @param ParamEstimates Estimates obtained for the parameters. Following ordering of ParamNames
#' @param gamma The lower threshold for the pareto function. Entered into the function so expression can be evaluated.
#' @param Nu The upper threshold for the pareto function. Entered into the function so expression can be evaluated.
#' @return Returns confidence bands for each of the variables, and the variance-covariance matrix
#' @examples
#' SandwichMaker(x=Value, w=Weights, EXPRESS=expression(w(log(1/sigma)-((1+alpha)/alpha)log((sigma+alpha*(x-gamma.gep.hat1))/sigma))), ParamNames=c("alpha", "sigma"),ParamEstimates=c(1.8,1.2),gamma=100,Nu=1000000)

SandwichMaker <- function(x, w=1, EXPRESS, ParamNames, ParamEstimates,gamma=100,Nu=NA){
  # x: Value inputs
  # w: Weight inputs
  # EXPRESS: Maximum Likelihood Equation, the first derivative of the likelihood function, which is then differentiated with respect to the parameters to construct score and hessian matrices. Enter as expression() input.
  #     EXPRESS example: EXPRESS <- expression(w*(log(1/sigma) - ((1+alpha)/alpha)*log((sigma +alpha*(x-100.0129))/sigma))) # Generalised Pareto
  # ParamNames: Vector of strings of names of parameters in EXPRESS, that function is differentiated with respect to.
  #     ParamNames example: ParamNames <- c("alpha", "sigma")
if (is.na(Nu)){nu=max(x)}
  # ParamEstimates: Estimates for values of parameters which are constructing bands on
  if (exists("ScoreNames")){rm(ScoreNames)}
  if (exists("HessianNames")){rm(HessianNames)}

  # Calculate the score matrix (first derivatives that we are optimizing)
  for (i in ParamNames) {
    df <- D(EXPRESS, i)
    assign(paste0("Score.",i),df)
    if (!exists("ScoreNames")){ScoreNames <- c()}
    ScoreNames <- c(ScoreNames, paste0("Score.",i))
  }

  # Calculate the Hessian matrix (second derivatives that measure the rate of change away from optimal)
  for (i in ParamNames) {
    for (j in ScoreNames){
      df <- D(get(j), i)
      assign(paste0("Hessian.",gsub("Score.","",j),".",i),df)
      if (!exists("HessianNames")){HessianNames <- c()}
      HessianNames <- c(HessianNames,paste0("Hessian.",gsub("Score.","",j),".",i))
    }
  }

  # Variance estimator
  vfun <- function(ParamEstimates, ParamNames, x, w) {
    if (length(ParamEstimates)==length(ParamNames)){
      # Evaluate Score Functions
      for (i in 1:length(ParamNames)){
        val <- ParamEstimates[i]
        assign(ParamNames[i],val)
      }
      # Estimate Variance-Covariance Matrix Components
      VCM <- expand.grid(rep(list(ScoreNames), 2))
      VCM[,3] <- NA
      for (i in 1:nrow(VCM)){
        VCM[i,3] <- sum(eval(get(paste(VCM[i,1]))) * eval(get(paste(VCM[i,2]))))
      }

      # Construct Variance-Covariance Matrix
      matrix(c(VCM[,3]), length(ParamEstimates), length(ParamEstimates))
    }
  }

  # This is penalising for misspecification (in MLE, is just I)
  jfun <- function(ParamEstimates, ParamNames, x, w) {
    if (length(ParamEstimates)==length(ParamNames)){

      for (i in 1:length(ParamNames)){
        val <- ParamEstimates[i]
        assign(ParamNames[i],val)
      }

      # Estimate Penalty Matrix Components
      VCM <- expand.grid(list(HessianNames))
      VCM[,2] <- NA
      for (i in 1:nrow(VCM)){
        VCM[i,2] <- sum(-eval(get(paste(VCM[i,1]))) )
      }

      # Construct Penalty Matrix Matrix
      matrix(c(VCM[,2]), length(ParamEstimates), length(ParamEstimates))
    }
  }

  Vhat <- vfun(ParamEstimates, ParamNames, x, w)
  Jhat <- jfun(ParamEstimates, ParamNames, x, w)

  # This is the sandwich
  Mhat <- solve(Jhat) %*% Vhat %*% solve(Jhat)
  crit <- qnorm((1 + 0.95)/2)

  results <- list()

  for (i in 1:length(ParamNames)) {
    results[[ParamNames[i]]] <- ParamEstimates[i] +  crit * c(-1, 1) * sqrt(Mhat[i, i])
  }
  results[["VCM"]] <- Mhat

  return(results)
}
