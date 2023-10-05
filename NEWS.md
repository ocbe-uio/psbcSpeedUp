# TODO in next versions

* use Bandicoot: Armadillo C++ library for GPU
* remove Rcpp::List object, instead writing all output objects into files
* allow multiple chains and omp for parallelisation
* add R functions for feature (stability) selection
* add R functions for survival predictions
* extend C++ source code for implementing other shrinkage and group priors

# psbcSpeedUp 2.0.4

* Improved C++ code to speedup 10-15%
* Added configuration files for checking omp, compilers etc.

# psbcSpeedUp 2.0.3

* First CRAN version
* Added adaptve jumping rule
* Fixed bug for q=0 to allowed variable selection for all covariates by default

# psbcSpeedUp 2.0.2

* Optimised the partition of the time axis for gamma process prior

# psbcSpeedUp 2.0.1

* Set armadillio random values and RNG

# psbcSpeedUp 2.0

* First GitHub commit version
