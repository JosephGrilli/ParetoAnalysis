#' Hardcoded confidence bounds for pseudo-maximum likelihood
#'
#' @description Calculates confidence bounds for pseudo-maximum likelihood using sandwich variance-covariance matrix (Geyer (2013): https://www.stat.umn.edu/geyer/5601/notes/sand.pdf).
#' @param x Data values used in parameter estimation.
#' @param w Data weights used in estimation. Default = 1.
#' @param gamma Lower threshold parameter for pareto. Default is minimum of x.
#' @param nu Upper threshold parameter for pareto. Default is maximum of x.
#' @param specification Evaluate sandwich estimate for "Type 1" or "Generalized" (hard coded because of calling gamma).
#' @param truncated Evaluate sandwich estimate for truncated (TRUE) or untruncated (FALSE) (hard coded because of calling gamma).
#' @param ParamNames The names of the parameters used in the expression. Example for generalised pareto: c("alpha", "sigma")
#' @param ParamEstimates Estimates obtained for the parameters. Following ordering of ParamNames
#' @return Returns confidence bands for each of the variables, and the variance-covariance matrix
#' @export
#' @examples
#' \dontrun{
#' HSandwichMaker(x=Value, w=Weights,gamma=100,nu=max(Value) specification="Generalized", truncated=FALSE, ParamNames=c("alpha", "sigma"),ParamEstimates=c(1.8,1.2))
#' }

HSandwichMaker <- function(x, w=1,gamma=NULL,nu=NULL, specification, truncated, ParamNames, ParamEstimates){
  # x: Value inputs
  # w: Weight inputs
  # EXPRESS: Maximum Likelihood Equation, the first derivative of the likelihood function, which is then differentiated with respect to the parameters to construct score and hessian matrices. Enter as expression() input.
  #     EXPRESS example: EXPRESS <- expression(w*(log(1/sigma) - ((1+alpha)/alpha)*log((sigma +alpha*(x-100.0129))/sigma))) # Generalised Pareto
  # ParamNames: Vector of strings of names of parameters in EXPRESS, that function is differentiated with respect to.
  #     ParamNames example: ParamNames <- c("alpha", "sigma")
  # ParamEstimates: Estimates for values of parameters which are constructing bands on
  if (exists("ScoreNames")){rm(ScoreNames)}
  if (exists("HessianNames")){rm(HessianNames)}
  if (is.null(gamma)){gamma <- min(x)}
  if (is.null(nu)){nu <- max(x)}

  for (i in 1:length(ParamNames)){
    df <- ParamEstimates[i]
    assign(ParamNames[i],df)
    rm(df)
  }

  # Calculate the score matrix (first derivatives that we are optimizing)
  # for (i in ParamNames) {
  #   df <- D(EXPRESS, i)
  #   assign(paste0("Score.",i),df)
  #   if (!exists("ScoreNames")){ScoreNames <- c()}
  #   ScoreNames <- c(ScoreNames, paste0("Score.",i))
  # }

  # Do this for each combination of specifications and truncation (take the expression and names from the main code, run it through here line by line, then store the output)
  if(specification=="Type 1" & truncated==FALSE){
    Score.alpha <- expression(w * (1/alpha + (log(gamma) - log(x))))
  }
  if(specification=="Type 1" & truncated==TRUE){
    Score.alpha <- expression(w * (1/alpha + (log(gamma) - log(x)) + (gamma/nu)^alpha *log((gamma/nu))/(1 - (gamma/nu)^alpha)))
    }
  if(specification=="Generalized" & truncated==FALSE){
    Score.alpha <- expression(-(w * ((1/alpha - (1 + alpha)/alpha^2) * log((sigma + alpha * (x - gamma))/sigma) + ((1 + alpha)/alpha) * ((x - gamma)/sigma/((sigma + alpha * (x - gamma))/sigma)))))
    Score.sigma <- expression(-(w * (1/sigma^2/(1/sigma) + ((1 + alpha)/alpha) * ((1/sigma - (sigma + alpha * (x - gamma))/sigma^2)/((sigma + alpha * (x - gamma))/sigma)))))
  }
  if(specification=="Generalized" & truncated==TRUE){
    Score.alpha <- expression(-(w * ((1/alpha - (1 + alpha)/alpha^2) * log((sigma + alpha * (x - gamma))/sigma) + ((1 + alpha)/alpha) * ((x - gamma)/sigma/((sigma + alpha * (x - gamma))/sigma)) -
                                       (((nu - gamma)/sigma)^(-1/alpha) + alpha * (((nu - gamma)/sigma)^(-1/alpha) * (log(((nu - gamma)/sigma)) * (1/alpha^2))))/(1 - (1 + alpha * ((nu - gamma)/sigma)^(-1/alpha))))))
    Score.sigma <- expression(-(w * (1/sigma^2/(1/sigma) + ((1 + alpha)/alpha) * ((1/sigma - (sigma + alpha * (x - gamma))/sigma^2)/((sigma + alpha * (x - gamma))/sigma)) + alpha * (((nu - gamma)/sigma)^((-1/alpha) - 1) *
                                     ((-1/alpha) * ((nu - gamma)/sigma^2)))/(1 - (1 + alpha * ((nu - gamma)/sigma)^(-1/alpha))))))
  }

  # Calculate the Hessian matrix (second derivatives that measure the rate of change away from optimal)
  # for (i in ParamNames) {
  #   for (j in ScoreNames){
  #     df <- D(get(j), i)
  #     assign(paste0("Hessian.",gsub("Score.","",j),".",i),df)
  #     if (!exists("HessianNames")){HessianNames <- c()}
  #     HessianNames <- c(HessianNames,paste0("Hessian.",gsub("Score.","",j),".",i))
  #   }
  # }

  # Do this for each combination of specifications and truncation (take the expression and names from the main code, run it through here line by line, then store the output)
  if(specification=="Type 1" & truncated==FALSE){Hessian.alpha.alpha <- expression(-(w * (1/alpha^2)))}
  if(specification=="Type 1" & truncated==TRUE){
    Hessian.alpha.alpha <- expression(w * ((gamma/nu)^alpha * log((gamma/nu)) * log((gamma/nu))/(1 - (gamma/nu)^alpha) + (gamma/nu)^alpha *
                                                                                         log((gamma/nu)) * ((gamma/nu)^alpha * log((gamma/nu)))/(1 - (gamma/nu)^alpha)^2 - 1/alpha^2))
  }
  if(specification=="Generalized" & truncated==FALSE){
    Hessian.alpha.alpha <- expression(-(w * ((1/alpha - (1 + alpha)/alpha^2) * ((x - gamma)/sigma/((sigma +
                                                                                                      alpha * (x - gamma))/sigma)) - (1/alpha^2 + (1/alpha^2 -
                                                                                                                                                     (1 + alpha) * (2 * alpha)/(alpha^2)^2)) * log((sigma + alpha *
                                                                                                                                                                                                      (x - gamma))/sigma) + ((1/alpha - (1 + alpha)/alpha^2) *
                                                                                                                                                                                                                               ((x - gamma)/sigma/((sigma + alpha * (x - gamma))/sigma)) -
                                                                                                                                                                                                                               ((1 + alpha)/alpha) * ((x - gamma)/sigma * ((x -
                                                                                                                                                                                                                                                                              gamma)/sigma)/((sigma + alpha * (x - gamma))/sigma)^2)))))
    Hessian.sigma.alpha <- expression(-(w * ((1/alpha - (1 + alpha)/alpha^2) * ((1/sigma - (sigma +
                                                                                              alpha * (x - gamma))/sigma^2)/((sigma + alpha *
                                                                                                                                (x - gamma))/sigma)) - ((1 + alpha)/alpha) * ((x -
                                                                                                                                                                                 gamma)/sigma^2/((sigma + alpha * (x - gamma))/sigma) +
                                                                                                                                                                                (1/sigma - (sigma + alpha * (x - gamma))/sigma^2) *
                                                                                                                                                                                ((x - gamma)/sigma)/((sigma + alpha * (x - gamma))/sigma)^2))))
    Hessian.alpha.sigma <- expression(-(w * ((1/alpha - (1 + alpha)/alpha^2) * ((1/sigma - (sigma +
                                                                                              alpha * (x - gamma))/sigma^2)/((sigma + alpha *
                                                                                                                                (x - gamma))/sigma)) - ((1 + alpha)/alpha) * ((x -
                                                                                                                                                                                 gamma)/sigma^2/((sigma + alpha * (x - gamma))/sigma) +
                                                                                                                                                                                (x - gamma)/sigma * (1/sigma - (sigma + alpha *
                                                                                                                                                                                                                  (x - gamma))/sigma^2)/((sigma + alpha * (x -
                                                                                                                                                                                                                                                             gamma))/sigma)^2))))
    Hessian.sigma.sigma <- expression(w * (((1 + alpha)/alpha) * ((1/sigma^2 + (1/sigma^2 - (sigma +
                                                                                               alpha * (x - gamma)) * (2 * sigma)/(sigma^2)^2))/((sigma +
                                                                                                                                                    alpha * (x - gamma))/sigma) + (1/sigma - (sigma +
                                                                                                                                                                                                alpha * (x - gamma))/sigma^2) * (1/sigma - (sigma +
                                                                                                                                                                                                                                              alpha * (x - gamma))/sigma^2)/((sigma + alpha *
                                                                                                                                                                                                                                                                                (x - gamma))/sigma)^2) + (2 * sigma/(sigma^2)^2/(1/sigma) -
                                                                                                                                                                                                                                                                                                            1/sigma^2 * (1/sigma^2)/(1/sigma)^2)))
  }
  if(specification=="Generalized" & truncated==TRUE){
    Hessian.alpha.alpha <- expression(-(w * ((1/alpha - (1 + alpha)/alpha^2) * ((x - gamma)/sigma/((sigma + alpha * (x - gamma))/sigma)) - (1/alpha^2 + (1/alpha^2 -
                                     (1 + alpha) * (2 * alpha)/(alpha^2)^2)) * log((sigma + alpha * (x - gamma))/sigma) + ((1/alpha - (1 + alpha)/alpha^2) *
                                     ((x - gamma)/sigma/((sigma + alpha * (x - gamma))/sigma)) - ((1 + alpha)/alpha) * ((x - gamma)/sigma * ((x - gamma)/sigma)/((sigma + alpha * (x - gamma))/sigma)^2)) -
                                     ((((nu - gamma)/sigma)^(-1/alpha) * (log(((nu - gamma)/sigma)) * (1/alpha^2)) + ((((nu - gamma)/sigma)^(-1/alpha) *
                                     (log(((nu - gamma)/sigma)) * (1/alpha^2))) + alpha * (((nu - gamma)/sigma)^(-1/alpha) * (log(((nu -
                                     gamma)/sigma)) * (1/alpha^2)) * (log(((nu - gamma)/sigma)) * (1/alpha^2)) - ((nu - gamma)/sigma)^(-1/alpha) *
                                     (log(((nu - gamma)/sigma)) * (2 * alpha/(alpha^2)^2)))))/(1 - (1 + alpha * ((nu - gamma)/sigma)^(-1/alpha))) +
                                     (((nu - gamma)/sigma)^(-1/alpha) + alpha * (((nu - gamma)/sigma)^(-1/alpha) * (log(((nu - gamma)/sigma)) *
                                     (1/alpha^2)))) * (((nu - gamma)/sigma)^(-1/alpha) + alpha * (((nu - gamma)/sigma)^(-1/alpha) *
                                     (log(((nu - gamma)/sigma)) * (1/alpha^2))))/(1 - (1 + alpha * ((nu - gamma)/sigma)^(-1/alpha)))^2))))
    Hessian.sigma.alpha <- expression(-(w * ((1/alpha - (1 + alpha)/alpha^2) * ((1/sigma - (sigma + alpha * (x - gamma))/sigma^2)/((sigma + alpha *
                                      (x - gamma))/sigma)) - ((1 + alpha)/alpha) * ((x - gamma)/sigma^2/((sigma + alpha * (x - gamma))/sigma) +
                                      (1/sigma - (sigma + alpha * (x - gamma))/sigma^2) * ((x - gamma)/sigma)/((sigma + alpha * (x - gamma))/sigma)^2) +
                                      (((((nu - gamma)/sigma)^((-1/alpha) - 1) * ((-1/alpha) * ((nu - gamma)/sigma^2))) + alpha * (((nu - gamma)/sigma)^((-1/alpha) -
                                      1) * (log(((nu - gamma)/sigma)) * (1/alpha^2)) * ((-1/alpha) * ((nu - gamma)/sigma^2)) + ((nu -
                                      gamma)/sigma)^((-1/alpha) - 1) * (1/alpha^2 * ((nu - gamma)/sigma^2))))/(1 - (1 + alpha *
                                      ((nu - gamma)/sigma)^(-1/alpha))) + alpha * (((nu - gamma)/sigma)^((-1/alpha) - 1) * ((-1/alpha) *
                                      ((nu - gamma)/sigma^2))) * (((nu - gamma)/sigma)^(-1/alpha) + alpha * (((nu - gamma)/sigma)^(-1/alpha) * (log(((nu -
                                      gamma)/sigma)) * (1/alpha^2))))/(1 - (1 + alpha * ((nu - gamma)/sigma)^(-1/alpha)))^2))))
    Hessian.alpha.sigma <- expression(-(w * ((1/alpha - (1 + alpha)/alpha^2) * ((1/sigma - (sigma + alpha * (x - gamma))/sigma^2)/((sigma + alpha *
                                      (x - gamma))/sigma)) - ((1 + alpha)/alpha) * ((x - gamma)/sigma^2/((sigma + alpha * (x - gamma))/sigma) +
                                      (x - gamma)/sigma * (1/sigma - (sigma + alpha * (x - gamma))/sigma^2)/((sigma + alpha * (x -
                                      gamma))/sigma)^2) + ((alpha * (((nu - gamma)/sigma)^(-1/alpha) * ((nu - gamma)/sigma^2/((nu - gamma)/sigma) *
                                      (1/alpha^2)) + ((nu - gamma)/sigma)^((-1/alpha) - 1) * ((-1/alpha) * ((nu - gamma)/sigma^2)) * (log(((nu -
                                      gamma)/sigma)) * (1/alpha^2))) + ((nu - gamma)/sigma)^((-1/alpha) - 1) * ((-1/alpha) * ((nu - gamma)/sigma^2)))/(1 -
                                      (1 + alpha * ((nu - gamma)/sigma)^(-1/alpha))) + (((nu - gamma)/sigma)^(-1/alpha) + alpha * (((nu -
                                      gamma)/sigma)^(-1/alpha) * (log(((nu - gamma)/sigma)) * (1/alpha^2)))) * (alpha * (((nu - gamma)/sigma)^((-1/alpha) -
                                      1) * ((-1/alpha) * ((nu - gamma)/sigma^2))))/(1 - (1 + alpha * ((nu - gamma)/sigma)^(-1/alpha)))^2))))
    Hessian.sigma.sigma <- expression(w * (alpha * (((nu - gamma)/sigma)^((-1/alpha) - 1) * ((-1/alpha) * ((nu - gamma) * (2 * sigma)/(sigma^2)^2)) +
                                      ((nu - gamma)/sigma)^(((-1/alpha) - 1) - 1) * (((-1/alpha) - 1) * ((nu - gamma)/sigma^2)) * ((-1/alpha) *
                                      ((nu - gamma)/sigma^2)))/(1 - (1 + alpha * ((nu - gamma)/sigma)^(-1/alpha))) + alpha * (((nu - gamma)/sigma)^((-1/alpha) -
                                      1) * ((-1/alpha) * ((nu - gamma)/sigma^2))) * (alpha * (((nu - gamma)/sigma)^((-1/alpha) - 1) * ((-1/alpha) *
                                      ((nu - gamma)/sigma^2))))/(1 - (1 + alpha * ((nu - gamma)/sigma)^(-1/alpha)))^2 + (((1 + alpha)/alpha) *
                                      ((1/sigma^2 + (1/sigma^2 - (sigma + alpha * (x - gamma)) * (2 * sigma)/(sigma^2)^2))/((sigma + alpha * (x - gamma))/sigma) +
                                      (1/sigma - (sigma + alpha * (x - gamma))/sigma^2) * (1/sigma - (sigma + alpha * (x - gamma))/sigma^2)/((sigma +
                                      alpha * (x - gamma))/sigma)^2) + (2 * sigma/(sigma^2)^2/(1/sigma) - 1/sigma^2 * (1/sigma^2)/(1/sigma)^2))))
  }



  # Variance estimator
  # vfun <- function(ParamEstimates, ParamNames, x, w) {
  #   if (length(ParamEstimates)==length(ParamNames)){
  #     # Evaluate Score Functions
  #     for (i in 1:length(ParamNames)){
  #       val <- ParamEstimates[i]
  #       assign(ParamNames[i],val)
  #     }
  #     # Estimate Variance-Covariance Matrix Components
  #     VCM <- expand.grid(rep(list(ScoreNames), 2))
  #     VCM[,3] <- NA
  #     for (i in 1:nrow(VCM)){
  #       VCM[i,3] <- sum(eval(get(paste(VCM[i,1]))) * eval(get(paste(VCM[i,2]))))
  #     }
  #
  #     # Construct Variance-Covariance Matrix
  #     matrix(c(VCM[,3]), length(ParamEstimates), length(ParamEstimates))
  #   }
  # }
  #
  # # This is penalising for misspecification (in MLE, is just I)
  # jfun <- function(ParamEstimates, ParamNames, x, w) {
  #   if (length(ParamEstimates)==length(ParamNames)){
  #
  #     for (i in 1:length(ParamNames)){
  #       val <- ParamEstimates[i]
  #       assign(ParamNames[i],val)
  #     }
  #
  #     # Estimate Penalty Matrix Components
  #     VCM <- expand.grid(list(HessianNames))
  #     VCM[,2] <- NA
  #     for (i in 1:nrow(VCM)){
  #       VCM[i,2] <- sum(-eval(get(paste(VCM[i,1]))) )
  #     }
  #
  #     # Construct Penalty Matrix Matrix
  #     matrix(c(VCM[,2]), length(ParamEstimates), length(ParamEstimates))
  #   }
  # }
  #
  # Vhat <- vfun(ParamEstimates, ParamNames, x, w)
  # Jhat <- jfun(ParamEstimates, ParamNames, x, w)

  # Do this for each combination of specifications and truncation (take the expression and names from the main code, run it through here line by line, then store the output)
  if(specification=="Type 1"){
    Vhat <- matrix(data=(sum(eval(Score.alpha)*eval(Score.alpha))),1,1)
    Jhat <- matrix(data=sum(-eval(Hessian.alpha.alpha)),1,1)
  }
  if(specification=="Generalized"){
    Vhat <- matrix(data=c((sum(eval(Score.alpha)*eval(Score.alpha))),
                          (sum(eval(Score.sigma)*eval(Score.alpha))),
                          (sum(eval(Score.alpha)*eval(Score.sigma))),
                          (sum(eval(Score.sigma)*eval(Score.sigma)))),2,2)
    Jhat <- matrix(data=(c(sum(-eval(Hessian.alpha.alpha)),
                           sum(-eval(Hessian.sigma.alpha)),
                           sum(-eval(Hessian.alpha.sigma)),
                           sum(-eval(Hessian.sigma.sigma)))),2,2)
  }

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
