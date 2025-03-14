
# ParetoGrilli

<!-- badges: start -->
<!-- badges: end -->

The goal of ParetoGrilli is to provide a download source for the pareto analysis tool developed for the centralised approach.

## Installation

You can install the development version of ParetoGrilli like so:

``` r
devtools::install_github("JosephGrilli/ParetoAnalysis")
```

## Example

This is a basic example which shows you how to run the example cases for the Pareto analysis.If you have data for values, the code can be run in the most simple form (with all other options set to defaults):

``` r
library(ParetoGrilli)
ParetoAnalysis(inputValues=Value)
```

This can be expanded to be calculated for weighted observations:

``` r
library(ParetoGrilli)
ParetoAnalysis(inputValues=Value,inputWeights=Weight)
```

Furthermore, using x_0, pre.specifications, pre.top.trunc, calibrate.ALPHA, calibrate.SIGMA, the results can be overwritten to run the process using exogeneously determined parameters:

``` r
library(ParetoGrilli)
ParetoAnalysis(inputValues=Value,inputWeights=Weight,id=hid, x_0=1000, pre.specification="Generalized", pre.top.trunc="TRUE", calibrate.ALPHA=1.8,calibrate.SIGMA=1.2)
```


Alternatively, data can be simulated from a Pareto distribution. This was used to demonstrate how data drawn from a Pareto Type 1 DGP and then truncated can be estimated using the new methods to recover the untruncated parameter:

``` r
library(ParetoGrilli)
ParetoAnalysis(type="Type 1",trun = "Truncated")
```

The parameters used in the simulation of the data can also be altered:

``` r
library(ParetoGrilli)
ParetoAnalysis(x_0 = 100, type="Type 1",trun = "Truncated", simnum=10000, upperT = 0.85, lowerT=0.15,sim.ALPHA=1.8,sim.SIGMA=1.2)
```

There are many other possible parameter variations to produce different outcomes. The code has the ability to estimate parameters, draw sample populations from this estimated distribution, simulate data from a distribution, and estimate distributions using calibrated parameters. Other factors, such as the p-values for the truncation tests, Hill Estimator parameters, model selection crieria, sampling methods, and more can be altered by changing the default options:

``` r
library(ParetoGrilli)
ParetoAnalysis(inputValues=Value,inputWeights=Weight,id=hid,x_0="top5",method="Synths", graphs=TRUE, ksp=0.05, ttp = 0.05, HillEstimator=1, OptimSearch="Local", Sampling="Inverse", model.selection = "AIC", simulatePopulation = FALSE)
```
