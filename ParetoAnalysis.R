#' Runs centralised apprach (Zwijnenburg, Grilli & Engelbrecht, 2022) Pareto code
#'
#' @description Estimates pareto values from input data (or simulated data) with option to simulate data from estimated pareto using multiple approaches.
#' @param inputValues Vector (nx1) of values held by household
#' @param inputWeights: Vector (nx1) of weights of households
#' @param inputid: Vector (nx1) of ids of households
#' @param inputDemo: Matrix (nxm) of demographic variables of households
#' @param inputDemoName: Vector(1xm) of names of demographic variables
#' @param x_0: Threshold value (numeric) or strings "top10", "top5", "top1" to select top tail of input data. Any string input takes form "top##", with the number ## being extracted.
#' @param type: Simulation Option for Pareto Type 1 ("Type 1"), Generalized Pareto ("Generalized"), or EU-SILC ("eusilc") data
#' @param trun: Simulation Option for Untruncated ("Untruncated"), Truncated ("Truncated"), or randomly truncated ("randTruncated") data.
#   Untruncated: Leaves simulated data as is
#   Truncated: Truncates data between 15% and 85% values
#   randTruncated: Truncates between a lower (0,0.5) and upper (lower, 1) proportions
#' @param method: Draws from estimated Pareto distribution using added rich households ("Synths"), adjusted weights ("Calibrate"), or adjusted values ("ranDraw")
#   Synths: Draws households from space above maximum survey value, calculated estimated population of missing area.
#   Calibrate: Iteratively draw weights that match the empirical CDF to the theoretical CDF while also retaining demographic totals
#   RanDraw: Draw sum of weights in total values from estimated distribution, order, and allocate to existing households
#   None: Skips the correction stage and only runs as an estimation
#' @param simnum: Simulation Option for number of observations to be drawn from Pareto Type 1 or Generalized Pareto distributions
#' @param pre.specification: Adjustment Option to overwrite goodness-of-fit results and select distribution type ("Type 1" or "Generalized")
#' @param pre.top.trunc: Adjustment Option to overwrite truncation test and select if truncated (TRUE or FALSE)
#' @param graphs: TRUE/FALSE option to produce graphs
#' @param ksp: Kolomogorov-Smirnov test statistic p-value level of significance to reject the null hypothesis that the data is drawn from the distribution
#' @param ttp: Truncation test statistics p-value level of significance to reject null hypothesis that the distribution is not truncated
#' @param upperT: Option for simulated Truncated ("Truncated") data to vary the percentage truncated at top of the distribution
#' @param lowerT: Option for simulated Truncated ("Truncated") data to vary the percentage truncated at bottom of the distribution
#' @param sim.ALPHA: Option for simulated  data to vary the Pareto Shape Parameter
#' @param sim.SIGMA: Option for simulated  data to vary the Pareto Scale Parameter
#' @param loopmax: Maximum number of iterations in the Calibrate iterative step
#' @param calibrate.ALPHA: Option for setting Pareto Shape Parameter when using pre.specification
#' @param calibrate.SIGMA: Option for setting Pareto Scale Parameter when using pre.specification
#' @param HillEstimator: Option for setting Hill Estimator share [0,1] of data used in estimating shape parameters. If NULL, it is calculated as either 0.5, or optimised if sample is small enough. Default 1 (i.e. no Hill Estimator adjustment)
#' @param OptimSearch: Sets the search style used for likelihood optimisation. "Default" uses existing options, or can set all to "Global", "Local", or "Both" (runs global then local).
#' @param seeded: Sets the seed to have fixed random draws (default: 20200916). Can therefore also use to alter seed. If set to NULL, then seed is not set (so can be randomised outside of the code, or doesnt overwrite existing seed setting.)
#' @param Sampling: Method of drawing new observations in the tail. Default option is "Deterministic" (calculate value from CDF where next observation would occur). Alternative is "Inverse" (make m2 random draws in the uncovered region of the tail)
#' @param model.selection: Method for choosing the winning model. Default option is "Arms", which is a logical decision tree to choose the smallest model, favouring trunction, from those retained in KS test.
#                  Alternative methods: "Loglikelihood", "AIC", "BIC"; will choose based on the favoured model by these metrics conditional on a specified model being chosen (i.e. it lacks the capabilities
#                  to select no model). Generally, BIC is the preferred of these as it heavily penalises redundant parameters and so will better delineate whether a truncation parameter is needed (the
#                  penalty increases from 2 to log(N), which is larger for N>7)
#' @param simulatePopulation: If TRUE (default), households are simulated in the estimated distribution. If FALSE, code ends with Pareto estimates output.

#' @return Returns estimated values, variances, goodness of fit and information criterion for pareto model. If simulatePopulation=TRUE, also returns estimated data sampled from pareto using method. If graphs=TRUE, also provides test graphical output.
#' 
#' @import lmtest
#' @import ggplot2
#' @import simPop
#' @export
#' 
#' @examples
#' \dontrun{
#' ParetoAnalysis(inputValues=Value,inputWeights=Weight,id=hid)
#' ParetoAnalysis(type="Type 1",trun = "Truncated")
#' ParetoAnalysis(x_0 = 100, type="Type 1",trun = "Truncated", simnum=10000, upperT = 0.85, lowerT=0.15,sim.ALPHA=1.8,sim.SIGMA=1.2)
#' ParetoAnalysis(inputValues=Value,inputWeights=Weight,id=hid, pre.specification="Generalized", pre.top.trunc="TRUE", calibrate.ALPHA=1.8,calibrate.SIGMA=1.2)
#' ParetoAnalysis(inputValues=Value,inputWeights=Weight,id=hid,x_0="top5",method="Synths", graphs=TRUE, ksp=0.05, ttp = 0.05, HillEstimator=1, OptimSearch="Local", Sampling="Inverse", model.selection = "BIC", simulatePopulation = FALSE)
#' }

ParetoAnalysis <- function(inputValues = NULL, inputWeights = NULL, inputid = NULL, inputDemo = NULL, inputDemoName = NULL, x_0 = NULL, type, trun = "Untruncated", method = "RanDraw", simnum=1000, pre.specification = NULL, pre.top.trunc = NULL, graphs = FALSE, ksp = 0.01, ttp = 0.1, upperT = 0.85, lowerT=0.15, sim.ALPHA = 1.8, sim.SIGMA = 1.2, loopmax = 2000, calibrate.ALPHA=NULL, calibrate.SIGMA=NULL, HillEstimator=NULL, OptimSearch="Default", seeded=NULL, Sampling="Deterministic", model.selection = c("Arms", "Loglikelihood", "AIC", "BIC"), simulatePopulation = TRUE) {

  ################################################################################
  ### Option block
  ################################################################################
  # Fixing options
  if (all(length(model.selection)==1 & model.selection %in% c("Arms", "Loglikelihood", "AIC", "BIC"))) {
    model.selection <- match.arg(model.selection, c("Arms", "Loglikelihood", "AIC", "BIC"), several.ok = FALSE)
  } else {
    print("model.selection invalid. Default selected.")
    model.selection <- "Arms"
  }

  # Survey data
  if (!is.null(inputValues)) {
    simu <- as.data.frame(inputValues)
    colnames(simu) <- c("Value") # Renaming the first column to "Value"
    simu$Weights <- if(!is.null(inputWeights)){
      inputWeights} else {1} # Setting Weights to inputWeights if available, else 1
    simu$id <- if(!is.null(inputid)){
      inputid} else {NA} # Setting id to inputid if available, else NA

    if (!is.null(inputDemo)) {
      simu <- cbind(simu, inputDemo)
      if (!is.null(inputDemoName)) {
        colnames(simu) <- c("Value", "Weights", inputDemoName)
      }
    }
    simu <- simu[order(simu$Value, decreasing = TRUE),] # Sorting the data frame by "Value" in descending order

    if (is.null(x_0)) {
      x_0 <- "top10" # Setting x_0 to "top10" if it is NULL
    }
    if (is.numeric(x_0)) {
      simu <- simu[simu$Value >= x_0,] # Filtering rows where "Value" is greater than or equal to x_0
    }
    if (is.character(x_0)) {
      if (!grepl("top",x_0)){return(print("Threshold is not numeric or specifying 'top' then a number. Please correct."))}
      simu$cumsum <- cumsum(simu$Weights) # Creating a new column "cumsum" with cumulative sum of "Weights"
      threshold <- as.numeric(gsub("top","",x_0))/100 # Remove top specification, convert to share number
      simu <- simu[simu$cumsum < threshold * sum(simu$Weights),] # Applying threshold based on cumulative sum
      x_0 <- min(simu$Value) # Setting x_0 to the minimum value in the remaining data
      simu <- simu[, colnames(simu) != c("cumsum")] # Removing the "cumsum" column
    }

    # Removing nonpositive thresholds (function is undefined)
    if (x_0 <= 0) {
      x_0 <- min(simu$Value[simu$Value > 0]) # Setting x_0 to the minimum positive value
      simu <- simu[simu$Value >= x_0,] # Filtering rows where "Value" is greater than or equal to x_0
    }
  }

  # Threshold value
  if (is.null(x_0)) {
    x_0 <- 100 # Setting x_0 to 100 if it is NULL
  }

   ################################################################################
  ### Data Section
  ################################################################################
  
if(!is.null(seeded)){
    set.seed(seeded)
  }

  weird <- scales::trans_new("signed_log", transform=function(x) sign(x)*log(abs(x)), inverse=function(x) sign(x)*exp(abs(x)))

  if (is.null(inputValues)) {
    y <- runif(simnum)
    simu <- data.frame()
    if (type == "Generalized") {
      # Simulate Generalized Pareto distributed data
      alpha <- sim.ALPHA
      sigma <- sim.SIGMA
      simu <- as.data.frame( x_0 + (sigma/alpha)*((1-y)^-alpha) - (sigma/alpha) )
      colnames(simu) <- "Value"
      simu$Weights <- 1
    }

    if (type == "Type 1") {
      # Simulate Pareto Type 1 distributed data
      alpha <- sim.ALPHA
      simu <- as.data.frame(x_0*(1-y)^(-1/alpha))
      colnames(simu) <- "Value"
      simu$Weights <- 1
    }

    if (type == "eusilc"){
      # EU-SILC sample data from laeken
      required_packages <- c("laeken")
      new_packages <- required_packages[!(required_packages %in% installed.packages())]
      if (length(new_packages)) { install.packages(new_packages, dependencies = TRUE, quiet = TRUE)}
      sapply(required_packages, suppressWarnings(suppressMessages(require)), quietly = TRUE, character.only = TRUE)

      data(eusilc)
      simu <- data.frame(eusilc$eqIncome, eusilc$db090, eusilc$rb090, eusilc$hsize, eusilc$age, eusilc$pl030)
      colnames(simu) <- c("Value", "Weights", "gender", "hsize", "age.group", "economic.status")
    }

    # Order data
    simu <- data.frame(simu[order(simu$Value, decreasing = TRUE),])

    if (type == "eusilc"){
      # Select top 10% of EU-SILC income distribution
      simu$cumsum <- cumsum(simu$Weights)
      simu <- simu[simu$cumsum<0.1*sum(simu$Weights),]
      x_0 <- min(simu$Value)
      simu <- simu[,colnames(simu)!=c("cumsum")]
      simu$age.group <- floor(simu$age.group/5) - 2
      simu$age.group[simu$age.group<1] <- -1
      simu$age.group[simu$age.group>17] <- 17
      inputDemoName <- c("gender", "hsize", "age.group", "economic.status")
      inputDemo <- simu[,c("gender", "hsize", "age.group", "economic.status")]
    }

    if (trun == "randTruncated") {
      # Truncated simulated distribution at random points
      lower <- runif(1,0,0.5)
      upper <- runif(1,lower,1)
      simu <- simu[(round(lower*nrow(simu),0)):(round(upper*nrow(simu),0)),]
    }

    if (trun == "Truncated") {
      # Truncated simulated distribution at random points
      # Reversed order (upper to lower) because data is descending. Therefore, to trim the top %, need to retain the rest
      lower <- (1-lowerT)*nrow(simu)
      upper <- (1-upperT)*nrow(simu)
      if (upper<1) {upper==1}
      if (lower>nrow(simu)) {lower==nrow(simu)}
      #lower <- 0.15*nrow(simu)
      #upper <- 0.85*nrow(simu)
      simu <- simu[upper:lower,]
    }
    simu$id <- NA
  }
  simu <- simu[simu$Value>=x_0,] # The Pareto distribution support is [x_0, infty), with x_0>0

  simu$Rank <- cumsum(simu$Weights)/sum(simu$Weights)
  if (graphs==TRUE) {
    plot0 <- ggplot(simu, aes(x=log(Value), y=(Rank), size = Weights)) + geom_jitter(alpha = 0.5) +
      ggtitle("Zipf Plot: Data Points") + scale_y_log10() + ylab("1-F(x)")
  }

  ################################################################################
  ### Pareto Type 1
  ################################################################################

  # Weighted MLE (from Aban, Meerschaert & Panorska, 2006) analytical Solution
  # gamma.pt1.hat1 <- min(simu$Value)
  # alpha.pt1.hat1 <- N/(sum(simu$Weights*log(simu$Value))-N*log(gamma.pt1.hat1))

  N <- sum(simu$Weights)
  n <- nrow(simu)

  # Hill estimator performs better than MLE, but is still not bias/robust to selection of k
  # Weighted Hill Estimator (from Aban, Meerschaert & Panorska, 2006)
  # Selecting r to minimise KS Test p-value (i.e. maximise the goodness of fit) conditional on the specification.
  # Could also try mean excess function Clauset, Shalizi & Newman (2009) equation (3.11) propose KS test method,
  # but here would need a measure of the best from a likelihood ratio of selected contenders

  if (is.null(HillEstimator)){
    # Adding loops for smaller samples to get better results - for larger samples this takes too long
    if (nrow(simu)<=40000) {
      print("Sample size less than 40,000. Applying Hill Estimator Loops.")
      r <- HillLoopfunction(simu$Value,simu$Weights)
    } else {r <- n-1}
  } else {r <- n*HillEstimator}
  r <- min(r, n-1) # Condition that Hill Estimator must remove the smallest observation so that gamma = value[r+1] is defined. Alternatively, could skip the Hill Estimator or add an extra "gamma" observation at value[r+1] to replicate non-adjusted form.

  res <- Pfunction(x=simu$Value,w=simu$Weights,r=r)

  alpha.pt1.hat2 <- res["alpha"]
  gamma.pt1.hat2 <- res["gamma"]

  logL1 <- sum(simu$Weights)*log(alpha.pt1.hat2) + sum(simu$Weights)*alpha.pt1.hat2*log(gamma.pt1.hat2) - (alpha.pt1.hat2+1)*sum(simu$Weights*log(simu$Value))
  vermeulen <- lm(log(cumsum(simu$Weights) - 0.5) ~ log(simu$Value))[1]  # Vermeulen (2014) alpha estimate


  Pt1Sandwich <- SandwichMaker(simu$Value, simu$Weights, gamma=gamma.pt1.hat2, nu=NULL, expression(w*(log(alpha)-log(x)+alpha*(log(gamma.pt1.hat2)-log(x)))), c("alpha"),c(alpha.pt1.hat2))


  # # Specific bootstrap for Type 1 Analytical solution
  #   n <- length(simu$Value)
  #   observations <- seq(1,n)
  #   totalDraw <- sample(observations,n*loops,replace=TRUE)
  #   sampleDraw <- matrix(totalDraw, n, loops)
  #   resu <- matrix(NA,loops,4)
  #   for (i in 1:loops){
  #     SampleX <- simu$Value[sampleDraw[,i]]
  #     SampleW <- simu$Weights[sampleDraw[,i]]
  #
  #     SampleX <- SampleX[order(SampleX, decreasing = TRUE)]
  #     SampleW <- SampleW[order(SampleX, decreasing = TRUE)]
  #
  #     BootLace <- c(sum(SampleW[1:r])/(sum(SampleW[1:r]*(log(SampleX[1:r]) - log(SampleX[r+1])))), NA, min(SampleX), max(SampleX))
  #     resu[i,] <- unlist(BootLace)
  #   }
  #   colnames(resu) <- c("alpha","sigma","gamma","nu")
  #
  #   # Prevent infinite values from affecting draws
  #   resu <- resu[!is.infinite(rowSums(resu, na.rm=TRUE)),]
  #
  #   resu2 <- colMeans(resu)
  #   resu2 <- rbind(resu2,diag(var(resu)))
  #   resu2 <- rbind(resu2, c(mean(resu[,1]-alpha.pt1.hat2),
  #                         NA,
  #                         mean(resu[,3]-gamma.pt1.hat2),
  #                         mean(resu[,4]-max(simu$Value))))
  #   resu2 <- rbind(resu2, c(t.test(resu[,1]-alpha.pt1.hat2)["p.value"],
  #                           NA,
  #                           t.test(resu[,3]-gamma.pt1.hat2)["p.value"],
  #                           t.test(resu[,4]-max(simu$Value))["p.value"]))
  #   rownames(resu2) <- c("mean","variance","bias","bias significance")

  ################################################################################
  ### Pareto Type 1 Plots
  ################################################################################

  UPareto1Plots <- PPlots(function(x){(gamma.pt1.hat2*(1-x)^(-1/alpha.pt1.hat2))},
                          function(x){(gamma.pt1.hat2/x)^alpha.pt1.hat2},
                          "Pareto Type 1",
                          simu$Value, simu$Weights,
                          ksp,TRUE)

  pvalue.pt1 <- UPareto1Plots$pValue

  if (graphs==TRUE){
    plot1 <- UPareto1Plots$Zipf
    plot01 <- UPareto1Plots$KSplot
    plot001 <- UPareto1Plots$Residplot
  }


  ################################################################################
  ### Generalized Pareto Function
  ################################################################################

  res <- GPfunction(simu$Value,simu$Weights,min(simu$Value),OptimSearch=OptimSearch)

  alpha.gep.hat1 <- res["alpha"]
  sigma.gep.hat1 <- res["sigma"]
  gamma.gep.hat1 <- res["gamma"]
  nu <- res["nu"]

  logL2 <- sum(simu$Weights*(log(1/sigma.gep.hat1) - ((1+alpha.gep.hat1)/alpha.gep.hat1)*log((sigma.gep.hat1 +alpha.gep.hat1*(simu$Value-gamma.gep.hat1))/sigma.gep.hat1)))

  GepSandwich <- SandwichMaker(simu$Value, simu$Weights, expression(w*(log(1/sigma)-((1+alpha)/alpha)*log((sigma+alpha*(x-gamma.gep.hat1))/sigma))), c("alpha", "sigma"),c(alpha.gep.hat1,sigma.gep.hat1))

  ################################################################################
  ### Generalized Pareto Plots
  ################################################################################

  UGeneralPlots <- PPlots(function(x){(gamma.gep.hat1 + (sigma.gep.hat1/alpha.gep.hat1)*((1-x)^-alpha.gep.hat1) - (sigma.gep.hat1/alpha.gep.hat1) )},
                          function(x){(1+(alpha.gep.hat1*(x - gamma.gep.hat1)/sigma.gep.hat1))^(-1/alpha.gep.hat1)},
                          "Generalised Pareto",
                          simu$Value, simu$Weights,
                          ksp,TRUE)

  pvalue.gep <- UGeneralPlots$pValue

  if (graphs==TRUE){
    plot2 <- UGeneralPlots$Zipf
    plot02 <- UGeneralPlots$KSplot
    plot002 <- UGeneralPlots$Residplot
  }



  ################################################################################
  ### Truncated Type 1 Pareto
  ################################################################################

  res <- TPfunction(simu$Value,simu$Weights,min(simu$Value),max(simu$Value),OptimSearch=OptimSearch)

  alpha.tp1.hat2 <- res["alpha"]
  gamma.tp1.hat2 <- res["gamma"]
  nu <- res["nu"]

  C <- (r/n)*((simu$Value[r+1])^alpha.pt1.hat2)
  p1 <- exp(-n*C*(nu^(-alpha.pt1.hat2)))
  if (p1<= ttp) {
    print("Null hypothesis rejected, Truncated Pareto type 1 in preferred to Untruncated Pareto Type 1.")
  }
  if (p1> ttp) {
    print("Null hypothesis that Pareto type 1 is untruncated distribution cannot be rejected.")
  }

  logL3 <- sum(simu$Weights)*log(alpha.tp1.hat2) + sum(simu$Weights)*alpha.tp1.hat2*log(gamma.tp1.hat2) - sum(simu$Weights)*log(1-(gamma.tp1.hat2/nu)^alpha.tp1.hat2) - (alpha.tp1.hat2+1)*sum(simu$Weights*log(simu$Value))

  Tp1Sandwich <- SandwichMaker(simu$Value, simu$Weights, expression(w*(log(alpha)-log(x)+alpha*(log(gamma.pt1.hat2)-log(x))-log(1-(gamma.tp1.hat2/nu)^alpha))), c("alpha"),c(alpha.tp1.hat2))

  ################################################################################
  ### Truncated Type 1 Pareto Plots
  ################################################################################

  TPareto1Plots <- PPlots(function(x){gamma.tp1.hat2*(1-x*(1-(gamma.tp1.hat2/nu)^alpha.tp1.hat2))^(-1/alpha.tp1.hat2)},
                          function(x){(((gamma.tp1.hat2/x)^alpha.tp1.hat2)-((gamma.tp1.hat2/nu)^alpha.tp1.hat2))/(1-(gamma.tp1.hat2/nu)^alpha.tp1.hat2)},
                          "Truncated Pareto Type 1",
                          simu$Value, simu$Weights,
                          ksp,TRUE)

  pvalue.tp1 <- TPareto1Plots$pValue

  if (graphs==TRUE){
    plot3 <- TPareto1Plots$Zipf
    plot03 <- TPareto1Plots$KSplot
    plot003 <- TPareto1Plots$Residplot
  }



  ################################################################################
  ### Truncated Generalized Pareto
  ################################################################################

  res <- TGPfunction(simu$Value,simu$Weights,min(simu$Value),max(simu$Value),OptimSearch=OptimSearch)

  alpha.gtp.hat1 <- res["alpha"]
  sigma.gtp.hat1 <- res["sigma"]
  gamma.gtp.hat1 <- res["gamma"]
  nu <- res["nu"]

  # Derived from survival function and q-test in Aban, Meerschaert & Panorska (2006)
  p2 <- exp(-n*(((sigma.gtp.hat1+alpha.gep.hat1*nu-alpha.gep.hat1*gamma.gtp.hat1)/sigma.gtp.hat1)^(-1/alpha.gep.hat1)))
  if (p2<= ttp) {
    print("Null hypothesis rejected, Truncated Generalized Pareto in preferred to Untruncated Generalized Pareto.")
  }
  if (p2> ttp) {
    print("Null hypothesis that Generalized Pareto is untruncated distribution cannot be rejected.")
  }

  logL4 <- sum(simu$Weights*(log(1/sigma.gtp.hat1) - ((1+alpha.gtp.hat1)/alpha.gtp.hat1)*log((sigma.gtp.hat1 +alpha.gtp.hat1*(simu$Value-gamma.gtp.hat1))/sigma.gtp.hat1) - log(1-(((sigma.gtp.hat1 + alpha.gtp.hat1*(nu-gamma.gtp.hat1))/sigma.gtp.hat1)^(-1/alpha.gtp.hat1)))))

  GtpSandwich <- SandwichMaker(simu$Value, simu$Weights, expression(w*(log(1/sigma)-((1+alpha)/alpha)*log((sigma+alpha*(x-gamma.gtp.hat1))/sigma) - log(1-(1+alpha*((nu-gamma.gtp.hat1)/sigma)^(-1/alpha))))), c("alpha", "sigma"),c(alpha.gtp.hat1,sigma.gtp.hat1))

  ################################################################################
  ### Truncated Generalized Pareto Plots
  ################################################################################

  TGeneralPlots <- PPlots(function(x){gamma.gtp.hat1 + (sigma.gtp.hat1/alpha.gtp.hat1)*((1-x*(1-(((sigma.gtp.hat1+(alpha.gtp.hat1*nu)-(alpha.gtp.hat1*gamma.gtp.hat1))/sigma.gtp.hat1)^(-1/alpha.gtp.hat1))))^(-alpha.gtp.hat1)-1)},
                          function(x){((((sigma.gtp.hat1+alpha.gtp.hat1*x -alpha.gtp.hat1*gamma.gtp.hat1)/sigma.gtp.hat1)^(-1/alpha.gtp.hat1)) - (((sigma.gtp.hat1+alpha.gtp.hat1*nu -alpha.gtp.hat1*gamma.gtp.hat1)/sigma.gtp.hat1)^(-1/alpha.gtp.hat1)))/(1-((sigma.gtp.hat1+alpha.gtp.hat1*nu -alpha.gtp.hat1*gamma.gtp.hat1)/sigma.gtp.hat1)^(-1/alpha.gtp.hat1))},
                          "Truncated Generalised Pareto",
                          simu$Value, simu$Weights,
                          ksp,TRUE)

  pvalue.gtp <- TGeneralPlots$pValue

  if (graphs==TRUE){
    plot4 <- TGeneralPlots$Zipf
    plot04 <- TGeneralPlots$KSplot
    plot004 <- TGeneralPlots$Residplot
  }



  ################################################################################
  ### Output matrix to return
  ################################################################################

  output.estimates <- as.data.frame(matrix(data = NA, nrow=4, ncol=15))
  rownames(output.estimates) <- c("Pareto Type 1", "Generalized Pareto", "Truncated Pareto Type 1", "Truncated Generalized Pareto")
  colnames(output.estimates) <- c("alpha", "Variance (alpha)", "Upper CI (alpha)", "Lower CI (alpha)" , "gamma", "sigma", "Variance (sigma)", "Upper CI (sigma)", "Lower CI (sigma)", "nu", "KS Test p-value","Truncation Test p-value", "Log Likelihoods", "AIC", "BIC")
  output.estimates[1,] <- c(alpha.pt1.hat2, Pt1Sandwich$VCM[1,1], Pt1Sandwich$alpha[2], Pt1Sandwich$alpha[1] , gamma.pt1.hat2,             NA,                   NA,                   NA,                   NA, NA, pvalue.pt1, NA, logL1,(2*2)-(2*logL1),(2*log(N))-(2*logL1))
  output.estimates[2,] <- c(alpha.gep.hat1, GepSandwich$VCM[1,1], GepSandwich$alpha[2], GepSandwich$alpha[1] , gamma.gep.hat1, sigma.gep.hat1, GepSandwich$VCM[2,2], GepSandwich$sigma[2], GepSandwich$sigma[1], NA, pvalue.gep, NA, logL2,(2*3)-(2*logL2),(3*log(N))-(2*logL2))
  output.estimates[3,] <- c(alpha.tp1.hat2, Tp1Sandwich$VCM[1,1], Tp1Sandwich$alpha[2], Tp1Sandwich$alpha[1] , gamma.tp1.hat2,             NA,                   NA,                   NA,                   NA, nu, pvalue.tp1, p1, logL3,(2*3)-(2*logL3),(3*log(N))-(2*logL3))
  output.estimates[4,] <- c(alpha.gtp.hat1, GtpSandwich$VCM[1,1], GtpSandwich$alpha[2], GtpSandwich$alpha[1] , gamma.gtp.hat1, sigma.gtp.hat1, GtpSandwich$VCM[2,2], GtpSandwich$sigma[2], GtpSandwich$sigma[1], nu, pvalue.gtp, p2, logL4,(2*4)-(2*logL4),(4*log(N))-(2*logL4))

  # In calculating the BIC, the STATA user guide BIC note (pages 2 and 3) state that either the number of observations or the sum of weights can be used:
  # That is a deep question. If the observations really are independent, then you should use N = M.
  # If the observations within group are not just correlated but are duplicates of one another, and they
  # had to be so, then you should use M = G. Between those two extremes, you should probably
  # use a number between N and G, but determining what that number should be from measured
  # correlations is difficult. Using N = M is conservative in that, if anything, it overweights complexity.
  # Conservativeness, however, is subjective, too: using N = G could be considered more conservative
  # in that fewer constraints are being placed on the data.
  # When the estimated correlation is high, our reaction would be that using N = G is probably more
  # reasonable. Our first reaction, however, would be that using BIC to compare models is probably a
  # misuse of the measure.
  # Stata uses N = M. An informal survey of web-based literature suggests that N = M is the
  # popular choice.

  # In accordance with this, N=M is used. This is because of the use of importance weights, meaning the household represents others, rather
  # than the exact same household being observed on multiple occassions.

  ################################################################################
  ### Selecting Model
  ################################################################################
  print(output.estimates)

  # End point if not simulating households from the population
  if (simulatePopulation==FALSE){
    print(output.estimates)
    if (graphs == TRUE) { output <- list(output.estimates,
                                         plot0, plot1, plot2, plot3, plot4,
                                         plot01, plot02, plot03, plot04,
                                         plot001, plot002, plot003, plot004)
    names (output) <- c("Estimates",
                        "Data Points","Plot Check 1","Plot Check 2","Plot Check 3","Plot Check 4",
                        "KS Test 1","KS Test 2","KS Test 3","KS Test 4",
                        "Residual Plot 1","Residual Plot 2","Residual Plot 3","Residual Plot 4")
    }
    if (graphs != TRUE) { output <- list(output.estimates)
    names (output) <- c("Estimates")
    }
    return(output)
  }

  res <- modelPara(output.estimates=output.estimates,model.selection=model.selection,ksp=ksp,ttp=ttp)
  specification <- res["specification"]
  top.trunc <- res["top.trunc"]

  # If ARMS rejects all specifications, exit output
  if (model.selection=="Arms"&specification==" "){
    print("Error: No Pareto retained by KS Test. Either use model selection LogLikelihood, AIC, or BIC to select best performing
          model by information criterion or likelihood race. Alternatively, test for other thresholds that support pareto distributions,
          or specify a model form using pre.specification and pre.top.trunc and parameters using calibrate.ALPHA and calibrate.SIGMA.")
    if (graphs == TRUE & is.null(pre.specification) & is.null(pre.top.trunc)) {
      output <- list(plot0, plot1, plot2, plot3, plot4,
                     plot01, plot02, plot03, plot04)
      names(output) <- c("Data Points", "Plot Check 1", "Plot Check 2", "Plot Check 3", "Plot Check 4",
                         "KS Test 1", "KS Test 2", "KS Test 3", "KS Test 4")
    } else {output <- NULL}
    return(output)
  }

  print(paste("Specification:", specification))
  print(paste("Top Truncation:", top.trunc))

  # Options to overwrite identification method and specify distribution you would like
  if (!is.null(pre.specification)) {
    if (pre.specification=="Type 1"|pre.specification=="Generalized") {
      specification <- pre.specification
      print("Overwriting specification test results with prespecified distribution. Proceeding to top tail adjustment...")
    }
  }
  if (!is.null(pre.top.trunc)) {
    if (pre.top.trunc==TRUE|pre.top.trunc==FALSE) {
      top.trunc <- pre.top.trunc
      print("Overwriting truncation test results with prespecified truncation. Proceeding to top tail adjustment...")
    }
  }

  # Name distribution that is being used
  if (top.trunc==TRUE) { print(paste0("Fitting to Truncated ", specification, " Pareto Distribution."))}
  if (top.trunc==FALSE) { print(paste0("Fitting to Untruncated ", specification, " Pareto Distribution."))}

  ################################################################################
  ### Adjust Tail
  ################################################################################

  output.totals <- c()
  new.obs <- c()

  res <- parame(output.estimates=output.estimates,specification=specification,top.trunc=top.trunc,
                calibrate.ALPHA=calibrate.ALPHA,calibrate.SIGMA=calibrate.SIGMA)
  alpha <- res["alpha"]
  sigma <- res["sigma"]
  gamma <- res["gamma"]
  nu <- res["nu"]

  print(paste0("Tail Totals before adjustment: Total number of households: ", sum(simu$Weights[simu$Value>=gamma], na.rm=TRUE), "; Total Value: ", sum(simu$Weights[simu$Value>=gamma]*simu$Value[simu$Value>=gamma], na.rm=TRUE)))
  print(paste0("Tail Totals before adjustment: Total number of observations: ", length(simu$Weights[simu$Value>=gamma])))


  res <- estimatedParetoPopFunction(x=simu$Value,w=simu$Weights,alpha=alpha,gamma=gamma,sigma=sigma,specification=specification)
  robust.population <- res[[1]]
  popu <- res[[2]]
  m2 <- res[[3]]

  if (graphs==TRUE) {
    plot6 <- ggplot(robust.population, aes(x=log(Value), y=population, group = type, colour = type, fill = type)) +
      geom_area() +
      ggtitle("Population Estimation (1% to 100% of observations)") +
      ylab("Total Population/Population Composition")
  }

  if (method=="Synths") {
    new.obs <- synthSample(x=simu$Value,w=simu$Weights,id=simu$id,m2=m2,popu=popu,alpha=alpha,gamma=gamma,sigma=NA,Sampling=Sampling,specification=specification)
    if (graphs==TRUE){
      plot5 <- ggplot(new.obs, aes(x=log(Value), y=zipf, size=Weights.new, group=type, color=type, shape=type)) +
        geom_jitter(alpha=0.5) +
        scale_color_manual(values=c("black", "blue")) +
        scale_y_log10() + ggtitle("Zipf plot with synthetic households") +
        ylab("1-F(x)")
    }
  }
  if (method=="Calibrate"){
    new.obs <- calibrateSample(x=simu$Value,w=simu$Weights,id=simu$id,inputDemoName=inputDemoName,inputDemo=inputDemo,alpha=alpha,gamma=gamma,sigma=sigma,specification=specification,loopmax=loopmax)
    if (graphs==TRUE) {
      plot5 <- ggplot(new.obs, aes(x=log(Value), y=(Weights.new/Weights), size=Weights.new)) +
        geom_jitter(alpha=0.5) +
        geom_abline(intercept = 1, slope = 0, color = "red") +
        scale_color_manual(values=c("black", "blue")) +
        ggtitle("Size of weight adjustments to fit Pareto Tail") +
        ylab("Ratio of New Weights to Original Weights")
    }
  }
  if (method=="RanDraw") {
    new.obs <- randrawSample(x=simu$Value,w=simu$Weights,id=simu$id,m2=m2,popu=popu,alpha=alpha,gamma=gamma,sigma=sigma,Sampling=Sampling,specification=specification)
    if (graphs==TRUE) {
      plot5 <- ggplot(new.obs, aes(x=log(Value), y=(Value/Value.o), size=Weights.new)) +
        geom_jitter(alpha=0.5) +
        geom_abline(intercept = 1, slope = 0, color = "red") +
        scale_color_manual(values=c("black", "blue")) +
        ggtitle("Size of value adjustments to fit Pareto Tail") +
        ylab("Ratio of New Values to Original Values")
    }
  }



  ################################################################################
  ### Plotting CDF vs ECDF
  ################################################################################


  if (graphs==TRUE) {
    bins = 100
    if (method=="Calibrate"|method=="RanDraw") { distplot1 <- data.frame(c(seq(gamma, 2*max(new.obs$Value), length.out=bins)))}
    if (method=="Synths") { distplot1 <- data.frame(c(seq(gamma, max(new.obs$Value), length.out=bins)))}
    colnames(distplot1) <- "Value.high"
    distplot1$Value.low <- c(gamma, distplot1$Value.high[1:(length(distplot1$Value.high)-1)])
    distplot1 <- distplot1[2:nrow(distplot1),]
    if (specification=="Generalized") {distplot1$Weights.new <- (1-((sigma + alpha*(distplot1$Value.high) - alpha*gamma)/sigma)^(-1/alpha)) -
      (1-((sigma + alpha*(distplot1$Value.low) - alpha*gamma)/sigma)^(-1/alpha))}
    if (specification=="Type 1")      {distplot1$Weights.new <- (1-(gamma/(distplot1$Value.high))^alpha) -
      (1-(gamma/(distplot1$Value.low))^alpha)}
    distplot1$Weights.new <- distplot1$Weights.new *popu
    distplot1$type <- "Theoretical"

    if (method=="Calibrate"|method=="RanDraw") { distplot2 <- data.frame(c(seq(gamma, 2*max(new.obs$Value), length.out=bins)))}
    if (method=="Synths") { distplot2 <- data.frame(c(seq(gamma, max(new.obs$Value), length.out=bins)))}
    colnames(distplot2) <- "Value.high"
    distplot2$Value.low <- c(gamma, distplot2$Value.high[1:(length(distplot2$Value.high)-1)])
    distplot2 <- distplot2[2:nrow(distplot2),]
    for (i in 1:length(distplot2$Value.high)) {
      distplot2$Weights.new[i] <- sum(new.obs$Weights.new[new.obs$Value<=distplot2$Value.high[i] & new.obs$Value>distplot2$Value.low[i]])
    }
    distplot2$type <- "Post-Adjustment"

    if (method=="Calibrate"|method=="RanDraw") { distplot4 <- data.frame(c(seq(gamma, 2*max(new.obs$Value), length.out=bins)))}
    if (method=="Synths") { distplot4 <- data.frame(c(seq(gamma, max(new.obs$Value), length.out=bins)))}
    colnames(distplot4) <- "Value.high"
    distplot4$Value.low <- c(gamma, distplot4$Value.high[1:(length(distplot4$Value.high)-1)])
    distplot4 <- distplot4[2:nrow(distplot4),]
    for (i in 1:length(distplot4$Value.high)) {
      distplot4$Weights.new[i] <- sum(simu$Weights[simu$Value<=distplot4$Value.high[i] & simu$Value>distplot4$Value.low[i]])
    }
    distplot4$type <- "Pre-Adjustment"



    distplot3 <- distplot2
    distplot3$type <- "Difference"
    distplot3$Weights.new <- abs(distplot1$Weights.new - distplot2$Weights.new)
    distplot3$mid <- distplot3$Value.low # (distplot3$Value.high + distplot3$Value.low)/2
    distplot3$width <- (distplot3$Value.high - distplot3$Value.low)

    densplot <- rbind(distplot1, distplot2, distplot4)
    densplot$mid <- (densplot$Value.high + densplot$Value.low)/2
    densplot$width <- (densplot$Value.high - densplot$Value.low)


    plot7 <- ggplot(densplot, aes(x=mid, y=(Weights.new+1),
                                  colour=factor(type, levels= c("Pre-Adjustment","Post-Adjustment","Theoretical", "Difference")),
                                  fill=factor(type, levels= c("Pre-Adjustment","Post-Adjustment","Theoretical", "Difference")))) +
      geom_bar(stat="identity", width = densplot$width, position="dodge") +
      scale_color_manual(values=c("red", "dodgerblue2", "dodgerblue4", "gray38")) +
      scale_fill_manual(values=c("red", "dodgerblue2", "dodgerblue4", "gray38")) +
      geom_step(stat="identity",data=distplot3) +
      geom_vline(xintercept=max(simu$Value), linetype="dashed", color = "Red") +
      scale_y_log10() +
      ggtitle(paste0("Histogram Binwidth: ",round(((max(densplot$Value.high)-gamma)/bins), digits=0),
                     "\nTheoretical Population: ",round(sum(densplot$Weights.new[densplot$type=="Theoretical"]), digits=0),
                     "\nPre-Adjustment Population: ",round(sum(densplot$Weights.new[densplot$type=="Pre-Adjustment"]), digits=0),
                     "\nPost-Adjustment Population: ",round(sum(densplot$Weights.new[densplot$type=="Post-Adjustment"]), digits=0))) +
      ylab("Population: Weights + 1") + theme(legend.title=element_blank())
  }

  if (method=="Synths") {
    output.totals <- as.data.frame(matrix(data=NA, nrow = 4, ncol = 2))
    rownames(output.totals) <- c("Original","Synths","Total","Proportion")
    colnames(output.totals) <- c("Total Value","Weights")
    output.totals["Original","Total Value"] <- sum(new.obs$Value[new.obs$type=="Original"]*new.obs$Weights[new.obs$type=="Original"])
    output.totals["Original","Weights"] <- sum(new.obs$Weights[new.obs$type=="Original"])
    output.totals["Synths","Total Value"] <- sum(new.obs$Value[new.obs$type=="Synths"]*new.obs$Weights[new.obs$type=="Synths"])
    output.totals["Synths","Weights"] <- sum(new.obs$Weights[new.obs$type=="Synths"])
    output.totals["Total","Total Value"] <- sum(new.obs$Value*new.obs$Weights)
    output.totals["Total","Weights"] <- sum(new.obs$Weights)
    output.totals["Proportion","Total Value"] <- paste0(round(100 * output.totals["Synths","Total Value"]/output.totals["Total","Total Value"],digits=3),"%")
    output.totals["Proportion","Weights"] <- paste0(round(100 * output.totals["Synths","Weights"]/output.totals["Total","Weights"],digits=3),"%")
  }

  if (method=="Calibrate") {
    output.totals <- as.data.frame(matrix(data=NA, nrow = 3, ncol = 2))
    rownames(output.totals) <- c("Original","Adjusted","Growth")
    colnames(output.totals) <- c("Total Value","Weights")
    output.totals["Original","Total Value"] <- sum(new.obs$Value*new.obs$Weights)
    output.totals["Original","Weights"] <- sum(new.obs$Weights)
    output.totals["Adjusted","Total Value"] <- sum(new.obs$Value*new.obs$Weights.new)
    output.totals["Adjusted","Weights"] <- sum(new.obs$Weights.new)
    output.totals["Growth","Total Value"] <- paste0(round(100 * (output.totals["Adjusted","Total Value"]/output.totals["Original","Total Value"]-1),digits=3),"%")
    output.totals["Growth","Weights"] <- paste0(round(100 * output.totals["Adjusted","Weights"]/output.totals["Original","Weights"],digits=3),"%")
  }


  if (method=="RanDraw") {
    output.totals <- as.data.frame(matrix(data=NA, nrow = 3, ncol = 2))
    rownames(output.totals) <- c("Original","Adjusted","Growth")
    colnames(output.totals) <- c("Total Value","Weights")
    output.totals["Original","Total Value"] <- sum(new.obs$Value.o*new.obs$Weights)
    output.totals["Original","Weights"] <- sum(new.obs$Weights)
    output.totals["Adjusted","Total Value"] <- sum(new.obs$Value*new.obs$Weights)
    output.totals["Adjusted","Weights"] <- sum(new.obs$Weights)
    output.totals["Growth","Total Value"] <- paste0(round(100 * output.totals["Adjusted","Total Value"]/output.totals["Original","Total Value"],digits=3),"%")
    output.totals["Growth","Weights"] <- paste0(round(100 * output.totals["Adjusted","Weights"]/output.totals["Original","Weights"],digits=3),"%")
  }

  # Now need to sample from correct distribution.
  # Number missing households: (1-output$`KS Test 1`$data[output$`KS Test 1`$data$Type=="CDF",][1,2])*sum(output$`Data Points`$data$Weights)
  if (method!="None") {
    Output.Obs <- data.frame(new.obs$Value)
    colnames(Output.Obs) <- "Value"
    if(method=="Calibrate") {
      Output.Obs$Weights <- new.obs$Weights.new
      Output.Obs <- cbind(Output.Obs, inputDemo)
      Output.Obs$Weights.o <- new.obs$Weights}
    if(method=="Synths"|method=="RanDraw") {Output.Obs$Weights <- new.obs$Weights}

    print(paste0("Tail Totals after adjustment: Total number of households: ", sum(new.obs$Weights.new, na.rm=TRUE), "; Total Value: ", sum(new.obs$Weights.new*new.obs$Value, na.rm=TRUE)))
    print(paste0("Tail Totals after adjustment: Total number of observations: ", length(new.obs$Weights.new)))
  }

  ################################################################################
  ### Function Output, dependent on options
  ################################################################################

  print(output.estimates)
  print(output.totals)

  if (graphs == TRUE) { output <- list(output.estimates, output.totals, Output.Obs, new.obs,
                                       plot0, plot1, plot2, plot3, plot4,
                                       plot5, plot6, plot7,
                                       plot01, plot02, plot03, plot04,
                                       plot001, plot002, plot003, plot004)
  names (output) <- c("Estimates","Tail Sums","Observations", "Detailed Observations",
                      "Data Points","Plot Check 1","Plot Check 2","Plot Check 3","Plot Check 4",
                      "Adjustment Plot", "Population Estimates", "Density plots",
                      "KS Test 1","KS Test 2","KS Test 3","KS Test 4",
                      "Residual Plot 1","Residual Plot 2","Residual Plot 3","Residual Plot 4")
  }
  if (graphs != TRUE) { output <- list(output.estimates, output.totals, new.obs)
  names (output) <- c("Estimates","Tail Sums", "Detailed Observations")
  }
  return(output)
}
