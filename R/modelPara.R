#' Model selection code
#'
#' @description Part of the centralised approach code, taking the output table to identify the model specification and whether there is upper truncation.
#' @param output.estimates The table output from the pareto code providing parameter value estimates, test values, and model-selection criterion
#' @param model.selection The method used to select the model from the arms approach (see Zwijnenburg, Grilli & Engelbrecht, 2022), akaike information criterion, bayesian information criterion, or log likelihood. Default = "BIC".
#' @param ksp The level of significant p-value used to test fit of data to proposed distribution. Default is 0.01.
#' @param ttp The level of significant p-value used to test for evidence of finite upper truncation point. Default is 0.1.
#' @return Returns the selected specification (from "Type 1" or "Generalized") and if the distribution has an upper truncation value (TRUE or FALSE)
#' @examples
#' modelPara(output.estimates=output.estimates,model.selection="Arms",ksp=0.01,ttp=0.1)

modelPara <- function(output.estimates,model.selection="BIC",ksp = 0.01,ttp = 0.1){

  specification <- " "
  top.trunc <- FALSE

  if (model.selection=="Arms"){
    # Identify the model specification and parameter values based on test results
    Arm1 <- output.estimates$`KS Test p-value`[1] > ksp
    Arm2 <- output.estimates$`KS Test p-value`[2] > ksp
    Arm3 <- output.estimates$`KS Test p-value`[3] > ksp
    Arm4 <- output.estimates$`KS Test p-value`[4] > ksp
    Arm5 <- output.estimates$`Truncation Test p-value`[3] > ttp
    Arm6 <- output.estimates$`Truncation Test p-value`[4] > ttp

    if (Arm1 && Arm3 && Arm5) {
      specification <- "Type 1"
      top.trunc <- FALSE
    } else if (Arm1 && !Arm3) {
      specification <- "Type 1"
      top.trunc <- FALSE
    } else if (Arm1 && Arm3 && !Arm5) {
      specification <- "Type 1"
      top.trunc <- TRUE
    } else if (!Arm1 && Arm2 && !Arm4 && Arm6) {
      specification <- "Generalized"
      top.trunc <- FALSE
    } else if (!Arm1 && Arm2 && !Arm4 && !Arm6) {
      specification <- "Generalized"
      top.trunc <- FALSE
    } else if (!Arm1 && Arm2 && Arm4 && !Arm6) {
      specification <- "Generalized"
      top.trunc <- TRUE
    } else if (!Arm1 && Arm2 && Arm4 && Arm6) {
      specification <- "Generalized"
      top.trunc <- FALSE
    } else if (!Arm1 && !Arm2 && !Arm4 && !Arm3) {
      print("Error: No Pareto retained by KS Test. Potential misspecification. Suggest changing x_0.")
      specification==" "
      top.trunc <- NA
    } else if (!Arm1 && !Arm2 && !Arm4 && Arm6) {
      specification <- "Generalized"
      top.trunc <- TRUE
      print("Is Truncated Generalized. Logical Inconsistency that cannot reject null untruncated Generalized against the alternative of truncated Generalized. Potential that Truncation is very low or parameter estimate for untruncated model is very poorly defined.")
    } else if (!Arm1 && !Arm2 && Arm3 && !Arm4 && !Arm5) {
      specification <- "Type 1"
      top.trunc <- TRUE
    } else if (!Arm1 && !Arm2 && Arm3 && !Arm4 && Arm5) {
      specification <- "Type 1"
      top.trunc <- TRUE
      print("Is Truncated Type 1. Logical Inconsistency that cannot reject null untruncated Type 1 against the alternative of truncated Type 1. Potential that Truncation is very low or parameter estimate for untruncated model is very poorly defined.")
    } else if (!Arm1 && !Arm2 && !Arm3 && Arm4 && !Arm5) {
      specification <- "Generalized"
      top.trunc <- TRUE
    } else if (!Arm1 && !Arm2 && Arm3 && Arm4 && !Arm5) {
      specification <- "Type 1"
      top.trunc <- TRUE
    }

  }
  if (model.selection=="Loglikelihood"){
    winner <- which(output.estimates[,13]==max(output.estimates[,13]))
    if (winner==1|winner==3){specification <- "Type 1"}
    if (winner==2|winner==4){specification <- "Generalized"}
    if (winner==1|winner==2){top.trunc <- FALSE}
    if (winner==3|winner==4){top.trunc <- TRUE}
  }
  if (model.selection=="AIC"){
    winner <- which(output.estimates[,14]==min(output.estimates[,14]))
    if (winner==1|winner==3){specification <- "Type 1"}
    if (winner==2|winner==4){specification <- "Generalized"}
    if (winner==1|winner==2){top.trunc <- FALSE}
    if (winner==3|winner==4){top.trunc <- TRUE}
  }
  if (model.selection=="BIC"){
    winner <- which(output.estimates[,15]==min(output.estimates[,15]))
    if (winner==1|winner==3){specification <- "Type 1"}
    if (winner==2|winner==4){specification <- "Generalized"}
    if (winner==1|winner==2){top.trunc <- FALSE}
    if (winner==3|winner==4){top.trunc <- TRUE}
  }
  return(c("specification"=specification,"top.trunc"=top.trunc))
}
