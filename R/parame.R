#' Model parameter value code
#'
#' @description Part of the centralised approach code, returns parameter values based on estimates and selected model form.
#' @param output.estimates The table output from the pareto code providing parameter value estimates, test values, and model-selection criterion
#' @param specification The model specification parameters from output.estimates will be taken from ("Type 1" or "Generalized").
#' @param top.trunc The inclusion of the upper truncation parameter from output.estimates (TRUE or FALSE).
#' @param calibrate.ALPHA Used to set alpha parameter if not using results. Default is NULL.
#' @param calibrate.SIGMA Used to set sigma parameter if not using results. Default is NULL.
#' @return Returns the preferred specification parameter estimates for alpha, gamma, sigma, and nu.
#' @examples
#' parame(output.estimates=output.estimates,specification="Type 1",top.trunc=TRUE,calibrate.ALPHA=NULL,calibrate.SIGMA=NULL)

parame <- function(output.estimates,specification="Type 1",top.trunc=TRUE,calibrate.ALPHA=NULL,calibrate.SIGMA=NULL){

  if (specification=="Type 1") {
    if (top.trunc==TRUE) {
      alpha <- output.estimates["Truncated Pareto Type 1","alpha"]
      gamma <- output.estimates["Truncated Pareto Type 1","gamma"]
      sigma <- output.estimates["Truncated Pareto Type 1","sigma"]
      nu    <- output.estimates["Truncated Pareto Type 1","nu"]
    }
    if (top.trunc==FALSE) {
      alpha <- output.estimates["Pareto Type 1","alpha"]
      gamma <- output.estimates["Pareto Type 1","gamma"]
      sigma <- output.estimates["Pareto Type 1","sigma"]
      nu    <- output.estimates["Pareto Type 1","nu"]
    }
    if (!is.null(calibrate.SIGMA)){
      sigma <- calibrate.SIGMA
    }
    if (!is.null(calibrate.ALPHA)){
      alpha <- calibrate.ALPHA
    }
  }

  if (specification=="Generalized") {
    if (top.trunc==TRUE) {
      alpha <- output.estimates["Truncated Generalized Pareto","alpha"]
      gamma <- output.estimates["Truncated Generalized Pareto","gamma"]
      sigma <- output.estimates["Truncated Generalized Pareto","sigma"]
      nu    <- output.estimates["Truncated Generalized Pareto","nu"]
    }
    if (top.trunc==FALSE) {
      alpha <- output.estimates["Generalized Pareto","alpha"]
      gamma <- output.estimates["Generalized Pareto","gamma"]
      sigma <- output.estimates["Generalized Pareto","sigma"]
      nu    <- output.estimates["Generalized Pareto","nu"]
    }
    if (!is.null(calibrate.SIGMA)){
      sigma <- calibrate.SIGMA
    }
    if (!is.null(calibrate.ALPHA)){
      alpha <- calibrate.ALPHA
    }

  }
  return(c("alpha"=alpha,"sigma"=sigma,"gamma"=gamma,"nu"=nu))
}
