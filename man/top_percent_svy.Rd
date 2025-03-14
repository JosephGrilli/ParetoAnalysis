% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/top_percent_svy.R
\name{top_percent_svy}
\alias{top_percent_svy}
\title{Top Percentage of Survey}
\usage{
top_percent_svy(x, value, weight = 1)
}
\arguments{
\item{x}{Proportion of data to calculate the top percentage value from. Should be between 0 (minimum) and 1 (maximum).}

\item{value}{The series that you wish to calculate the top x share of.}

\item{weight}{The weight of observations in the data. Default = 1.}
}
\value{
Value of the top x share of values.
}
\description{
Function calculates the percentage of the total value above proportion x, included for weighted data.
}
\examples{
top_percent_svy(0.5,value=runif(100,0,1000),weight=round(runif(100,1,100),0))
}
