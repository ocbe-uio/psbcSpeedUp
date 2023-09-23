# psbcSpeedUp

This is a C++ speed-up and extended version of the R-pakcage [psbcGroup](https://CRAN.R-project.org/package=psbcGroup).
It implements the Bayesian Lasso Cox model ([Lee et al., 2011](https://doi.org/10.2202/1557-4679.1301)) and the Bayesian Lasso Cox with mandatory variables ([Zucknick et al., 2015](https://doi.org/10.1002/bimj.201400160)).
Bayesian Lasso Cox models with other shrinkage and group priors ([Lee et al., 2015](https://doi.org/10.1002/sam.11266)) are to be implemented later on.

## Installation

Install the latest development version from GitHub

```r
#install.packages("remotes")
remotes::install_github("ocbe-uio/psbcSpeedUp")
```

## Examples

### Run a Bayesian Lasso Cox with mandatory variables

Data set `exampleData` consists of six components:
survival times `t`,
event status `di`,
covariates `x`,
number of genomics variables `p`,
number of clinical variables `q` and
true effects of covariates `beta_true`.
See `?exampleData` for more information of the data.
To run a Bayesian Lasso Cox model for variable selection of the first $p$ genomics variables, one can specify the hyperparameter `groupInd = 1:p`.

```r
# Load the example dataset
data("exampleData", package = "psbcSpeedUp")
p <- exampleData$p
q <- exampleData$q
survObj <- exampleData[1:3]

# Set hyperparameters (see help file for specifying more hyperparameters)
mypriorPara <- list('groupInd'=1:p, 'eta0'=0.02, 'kappa0'=1, 'c0'=2, 'r'=10/9, 'delta'=1e-05, 'lambdaSq'=1, 'sigmaSq'= runif(1, 0.1, 10),
'beta.prop.var'=1, 'beta.clin.var'=1)

# run Bayesian Lasso Cox
library("psbcSpeedUp")
set.seed(123)
fitBayesCox <- psbcSpeedUp(survObj, p=p, q=q, hyperpar=mypriorPara,
nIter=1000, burnin=500, outFilePath="/tmp")
```

```
Running MCMC iterations ...
[##################################################] 100%
DONE, exiting!
```

### Plot posterior estimates of regression cofficients

The function `psbcSpeedUp::plot()` can show the posterior mean and 95% credible intervals of regression coefficients.

```r
plot(fitBayesCox)
```

![](https://github.com/ocbe-uio/psbcSpeedUp/blob/main/README_plot_beta.png)<!-- -->

## References

> Kyu Ha Lee, Sounak Chakraborty, Jianguo Sun (2011).
> Bayesian variable selection in semiparametric proportional hazards model for high dimensional survival data.
> _The International Journal of Biostatistics_, 7:1. DOI: [10.2202/1557-4679.1301](https://doi.org/10.2202/1557-4679.1301).

> Kyu Ha Lee, Sounak Chakraborty, Jianguo Sun (2015).
> Survival prediction and variable selection with simultaneous shrinkage and grouping priors.
> _Statistical Analysis and Data Mining_, 8:114-127. DOI:[10.1002/sam.11266](https://doi.org/10.1002/sam.11266).

> Manuela Zucknick, Maral Saadati, Axel Benner (2015).
> Nonidentical twins: Comparison of frequentist and Bayesian lasso for Cox models.
> _Biometrical Journal_, 57:959-981. DOI:[10.1002/bimj.201400160](https://doi.org/10.1002/bimj.201400160).
