## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, eval = FALSE)
options(rmarkdown.html_vignette.check_title = FALSE)

## ----install1, eval = FALSE---------------------------------------------------
#  install.packages("psbcSpeedUp")

## ----install2, eval = FALSE---------------------------------------------------
#  #install.packages("remotes")
#  remotes::install_github("ocbe-uio/psbcSpeedUp")

## ----results='hide', warning=FALSE--------------------------------------------
#  # Load the example dataset
#  data("exampleData", package = "psbcSpeedUp")
#  p <- exampleData$p
#  q <- exampleData$q
#  survObj <- exampleData[1:3]
#  
#  # Set hyperparameters (see help file for specifying more hyperparameters)
#  mypriorPara <- list('eta0'=0.02, 'kappa0'=1, 'c0'=2, 'r'=10/9, 'delta'=1e-05,
#  'lambdaSq'=1, 'sigmaSq'= runif(1, 0.1, 10), 'beta.prop.var'=1, 'beta.clin.var'=1)
#  
#  # run Bayesian Lasso Cox
#  library("psbcSpeedUp")
#  set.seed(123)
#  fitBayesCox <- psbcSpeedUp(survObj, p=p, q=q, hyperpar=mypriorPara,
#  nIter=1000, burnin=500, outFilePath="/tmp")

## ----fig.width=5, fig.height=8------------------------------------------------
#  plot(fitBayesCox)

## ----fig.width=6, fig.heigh=5-------------------------------------------------
#  plotBrier(fitBayesCox, times = 80)

## -----------------------------------------------------------------------------
#  predict(fitBayesCox, type = c("cumhazard", "survival"))

