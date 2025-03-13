#' Zipf and Kolmogorov-Smirnov graphical tests for Pareto fit
#'
#' @description Function tests the fit of input data to a specified distribution.
#' @param ValueEquation Single variable function defining the inverse CDF function for the distribution of interest. Used to sample data. Example for Pareto Type 1: function(x){(gamma.pt1.hat2*(1-x)^(-1/alpha.pt1.hat2))}.
#' @param RankEquation Single variable function defining the complementary cumulative distribution function (CCDF). Used to rank x from 0 to 1.  Example for Pareto Type 1: function(x){(gamma.pt1.hat2/x)^alpha.pt1.hat2}.
#' @param TypeName Name of distribution.
#' @param inputData_Values Data values you wish check to match distribution used for RankEquation.
#' @param inputData_Weights Weights attached to data points. Default = 1.
#' @param ksp The level of significant p-value used to test fit of data to proposed distribution. Default is 0.01.
#' @param graphs If TRUE, graphs of tests are generated.
#' @return Returns the Kolmogorov-Smirnov test p-value. When graphs=TRUE, also returns Zipf, Kolmogorov-Smirnov, and residual plots.
#'
#' @import stats
#'
#' @examples
#' PPlots(function(x){(gamma.pt1.hat2*(1-x)^(-1/alpha.pt1.hat2))},function(x){(gamma.pt1.hat2/x)^alpha.pt1.hat2},"Pareto Type 1",simu$Value, simu$Weights,ksp,TRUE)

# Generate Zipf and KSplots
PPlots <- function(ValueEquation, RankEquation, TypeName=NA, inputData_Values, inputData_Weights=1, ksp= 0.01, graphs=TRUE){


  if (graphs==TRUE) {
    # Added min() in case top case have multiple obs
    yplot <- seq(0,1-min(inputData_Weights[inputData_Values==
                                             max(inputData_Values)]/sum(inputData_Weights)),length.out=10000)

    df.2 <- ValueEquation(x=yplot)
    df.2 <- as.data.frame(df.2[order(df.2, decreasing = TRUE)])
    colnames(df.2) <- "Value"
    df.2$Rank <- RankEquation(x=df.2$Value)
    df.2$type <- as.character(TypeName)
    df.2$Weights <- 1

    df.3 <- data.frame("Value"=inputData_Values, "Weights"=inputData_Weights)
    df.3 <- df.3[order(df.3$Value,decreasing=TRUE),]

    df.3$type <- "Data"
    df.3$Rank <- cumsum(df.3$Weights)/sum(df.3$Weights)

    df <- rbind(df.2, df.3)
    df$Value <- log(df$Value)

    resid <- df.3[,c("Value","Rank","Weights")]
    colnames(resid) <- c("Data","Rank","Weights")
    resid$Estimated <- ValueEquation(1-resid$Rank)
    resid$Error <- resid$Data-resid$Estimated

    model <- lm(Data~Estimated, data=resid)
    bpTest <- lmtest::bptest(model)
    bpTest <- data.frame("Label"=paste0("Null Hypothesis of Homoskedastic Variance: ", bpTest$p.value>0.05), "X"=-Inf, "Y"=-Inf)


    plot001 <- ggplot(resid, aes(x=Data, y=Error, size=Weights)) +
      geom_jitter(alpha=0.5) +
      geom_hline(yintercept = 0, color="black", alpha=0.6) +
      geom_hline(yintercept = stats::weighted.mean(x=resid$Error,w=resid$Weights), color="blue", linetype=2) +
      ggtitle(paste0("Residual Plot: ", TypeName)) +
      xlab("Data Values") + ylab("Data Value - Estimated Value")

    plot001 <- plot001 +
      annotate("label", x=bpTest$X, y=bpTest$Y,hjust=0,vjust=0, size=2.3, label=bpTest$Label)

    rm(df.2, df.3)

    plot1 <- ggplot(df, aes(x=(Value), y=(Rank), size=Weights, group=type, color=type)) +
      geom_jitter(data = df[df$type=="Data",], alpha=0.5) +
      geom_line(data = df[df$type==TypeName,], size = 1.4) +
      scale_color_manual(values=c("black", "red")) +
      ggtitle(paste0("Zipf Plot: ", TypeName)) +
      scale_y_log10() + ylab("1-F(x)")


  }


  # Test derived from Krieger & Pfeffermann (1997), equation 18
  KS <- data.frame("Weights"=inputData_Weights,"Value"=inputData_Values)
  KS$CDF <- 1 - RankEquation(x=KS$Value)
  KS <- KS[order(KS$Value, decreasing=TRUE),]

  dist <- aggregate(Weights ~ Value, KS, sum)
  dist <- dist[order(dist$Value, decreasing=FALSE),]
  dist$ECDF <- (cumsum(dist$Weights)/sum(dist$Weights))
  KS <- KS[order(KS$Value, decreasing = FALSE),]

  KS$ECDF <- dist$ECDF[match(KS$Value,dist$Value)]

  D <- max(abs(KS$ECDF-KS$CDF))
  D.test <- ks.test(KS$ECDF, KS$CDF)

  if (graphs==TRUE) {
    KSgraph <- rbind(data.frame("Value"=KS$Value,"Density"=KS$CDF,"Type"="CDF"),
                     data.frame("Value"=KS$Value,"Density"=KS$ECDF,"Type"="ECDF"))
    plot01 <- ggplot(KSgraph, aes(x=log(Value), y=Density, group = Type, color = Type)) +
      geom_line(size=1.5) + ggtitle(paste(TypeName)) + ylab("F(x)")
  }

  pvalue.pt1 <- D.test$p.value

  # If p.value < ksp, reject the null hypothesis that the data matches the proposed form.
  if (pvalue.pt1>ksp) {
    print(paste0(TypeName, " retained."))
  }

  if (pvalue.pt1<=ksp) {
    print(paste0(TypeName," rejected."))
  }
  if (graphs==TRUE){
    output <- list("pValue"=pvalue.pt1,"Zipf"=plot1,"KSplot"=plot01,"Residplot"=plot001)
  } else {output <- list("pValue"=pvalue.pt1)}
  return(output)
}
