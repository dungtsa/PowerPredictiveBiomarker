# test

The tool is for power calculation of predictive biomarker. Calculation requires a series of parameters to to determine subgroup proportion and subgroup censoring rate. Depending on study type (prospective or retrospective study), parameter setting is different (check the reference).

## Features

* The shiny applictaion for power calculation of predictive biomarker generates a statistical plan to justify the sample size

## Installation

Simply run the following from an R console:

```r
if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("rstudio/PowerPredictiveBiomarker")
```

## Getting Started

```r
require("PowerPredictiveBiomarker")
PowerPredictiveBiomarker.shiny()
```
![snapshot of shiny app](inst/img/worldcup.png)

