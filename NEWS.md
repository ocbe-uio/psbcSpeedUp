# psbcSpeedUp 2.0.7

* Set default `eta0` and `kappa0` via (Weibull) parametric survival model
* Fixed a bug in function `PSBC::updateTau_GL_cpp()`

# psbcSpeedUp 2.0.6

* Fixed a bug in function `PSBC::updateBH_cpp()`
* Other minor changes

# psbcSpeedUp 2.0.5

* Added R functions for survival predictions, e.g. time-dependent Brier scores, integrated Brier score

# psbcSpeedUp 2.0.4

* For Linux machines, installed GitHub v2.0.4 seems twice faster than CRAN v2.0.4 due to `omp_set_*` in drive.cpp to solve CRAN pre-test NOTE `CPU-elapse times ratio`
* Improved C++ code to speedup 10-20%
* Declared functions inside PSBC class
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
