#' Top Percentage of Survey
#'
#' @description Function calculates the percentage of the total value above proportion x, included for weighted data.
#' @param x Proportion of data to calculate the top percentage value from. Should be between 0 (minimum) and 1 (maximum).
#' @param value The series that you wish to calculate the top x share of.
#' @param weight The weight of observations in the data. Default = 1.
#' @return Value of the top x share of values.
#' @examples
#' top_percent_svy(0.5,value=runif(100,0,1000),weight=round(runif(100,1,100),0))

top_percent_svy <- function(x, value, weight=1) {
  if(x>1|x<0){
    return(print("Value for x falls outside of [0,1] range"))
  }
  data <- data.frame("value"=value,"weight"=weight)
  data[order(data$value,decreasing=TRUE),]
  data$sumw <- cumsum(data$weight)

  Wpercentile <- sum(data$weight) * x

  .id <- data$sumw <= Wpercentile
  result <- sum(data$weight[.id] * data$value[.id])

  if (sum(data$weight[.id]) < Wpercentile) {
    idx <- sum(.id)
    result <- result + (Wpercentile - data$sumw[idx]) * data$value[idx + 1]
  }
  result <- as.numeric(result)
  return(result)
}
