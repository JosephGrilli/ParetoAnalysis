
# ParetoGrilli

<!-- badges: start -->
<!-- badges: end -->

The goal of ParetoGrilli is to provide a download source for the pareto analysis tool developed for the centralised approach. The methodology is described in [Zwijnenburg, Grilli & Engelbrecht (2022)](https://iariw.org/wp-content/uploads/2022/07/Jorret-Joseph-Pao-IARIW-2022.pdf) and [Zwijnenburg & Grilli (2024)](https://iariw.org/wp-content/uploads/2024/08/Draft_1_1_JZ3.pdf)

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

The full set of parameters, descriptions, and any relevent explanations is included in the code, but also printed below for completeness:
```r
inputValues # Vector (nx1) of values held by household.
inputWeights # Vector (nx1) of weights of households.
inputid # Vector (nx1) of ids of households.
inputDemo # Matrix (nxm) of demographic variables of households.
inputDemoName # Vector(1xm) of names of demographic variables.
x_0 # Threshold value (numeric) or strings "top10", "top5", "top1" to select top tail of input data. Any string input takes form "top##", with the number ## being extracted.
type # Simulation Option for Pareto Type 1 ("Type 1"), Generalized Pareto ("Generalized"), or EU-SILC ("eusilc") data.
trun # Simulation Option for Untruncated ("Untruncated"), Truncated ("Truncated"), or randomly truncated ("randTruncated") data.
     #   Untruncated: Leaves simulated data as is.
     #   Truncated: Truncates data between 15% and 85% values.
     #   randTruncated: Truncates between a lower (0,0.5) and upper (lower, 1) proportions.
method # Draws from estimated Pareto distribution using added rich households ("Synths"), adjusted weights ("Calibrate"), adjusted values ("ranDraw"), or skips the estimation ("None").
       #   Synths: Draws households from space above maximum survey value, calculated estimated population of missing area.
       #   Calibrate: Iteratively draw weights that match the empirical CDF to the theoretical CDF while also retaining demographic totals.
       #   RanDraw: Draw sum of weights in total values from estimated distribution, order, and allocate to existing households.
       #   None: Skips the correction stage and only runs as an estimation.
simnum # Simulation Option for number of observations to be drawn from Pareto Type 1 or Generalized Pareto distributions.
pre.specification # Adjustment Option to overwrite goodness-of-fit results and select distribution type ("Type 1" or "Generalized").
pre.top.trunc # Adjustment Option to overwrite truncation test and select if truncated (TRUE or FALSE).
graphs # TRUE/FALSE option to produce graphs.
ksp # Kolomogorov-Smirnov test statistic p-value level of significance to reject the null hypothesis that the data is drawn from the distribution.
ttp # Truncation test statistics p-value level of significance to reject null hypothesis that the distribution is not truncated.
upperT # Option for simulated Truncated ("Truncated") data to vary the percentage truncated at top of the distribution.
lowerT # Option for simulated Truncated ("Truncated") data to vary the percentage truncated at bottom of the distribution.
sim.ALPHA # Option for simulated data to vary the Pareto Shape Parameter. Will overwrite estimated parameter value if set (default is NULL).
sim.SIGMA # Option for simulated data to vary the Pareto Scale Parameter. Will overwrite estimated parameter value if set (default is NULL).
loopmax # Maximum number of iterations in the Calibrate iterative step.
calibrate.ALPHA # Option for setting Pareto Shape Parameter when using pre.specification.
calibrate.SIGMA # Option for setting Pareto Scale Parameter when using pre.specification.
HillEstimator # Option for setting Hill Estimator share [0,1] of data used in estimating shape parameters. If NULL, it is calculated as either 0.5, or optimised if sample is small enough. Default 1 (i.e. no Hill Estimator adjustment).
OptimSearch # Sets the search style used for likelihood optimisation. "Default" uses existing options, or can set all to "Global", "Local", or "Both" (runs global then local).
seeded # Sets the seed to have fixed random draws (default: 20200916). Can therefore also use to alter seed. If set to NULL, then seed is not set (so can be randomised outside of the code, or doesnt overwrite existing seed setting.).
Sampling # Method of drawing new observations in the tail. Default option is "Deterministic" (calculate value from CDF where next observation would occur). Alternative is "Inverse" (make m2 random draws in the uncovered region of the tail).
model.selection # Method for choosing the winning model. Default option is "Arms", which is a logical decision tree to choose the smallest model, favouring trunction, from those retained in KS test.
                # Alternative methods: "Loglikelihood", "AIC", "BIC"; will choose based on the favoured model by these metrics conditional on a specified model being chosen (i.e. it lacks the capabilities
                # to select no model). Generally, BIC is the preferred of these as it heavily penalises redundant parameters and so will better delineate whether a truncation parameter is needed (the
                # penalty increases from 2 to log(N), which is larger for N>7).
simulatePopulation # If TRUE (default), households are simulated in the estimated distribution. If FALSE, code ends with Pareto estimates output.
```
