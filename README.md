# test

The tool is for power calcuation of predictive predictive biomarker. Calculaiton requires a series of parameters to to determine subgroup proportion and subgroup censoring rate. Depending on study type (prospective or retrospective study), parameter setting is different (check the reference).

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
![World cup timeline](inst/img/worldcup.png)

