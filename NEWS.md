<div style="text-align: left;">

# psbcSpeedUp 2.1.0 (2025-09-30)

* Replace `pugixml.cpp` by `Rcpp::List`
* Convert line endings in configure.ac to LF
* Require C++17
* Remove redundant C++ code

# psbcSpeedUp 2.0.7 (2024-07-01)

* Set default `eta0` and `kappa0` via (Weibull) parametric survival model
* Fixed a bug in function `PSBC::updateTau_GL_cpp()`

# psbcSpeedUp 2.0.6 (2024-03-21)

* Fixed a bug in function `PSBC::updateBH_cpp()`
* Other minor changes

# psbcSpeedUp 2.0.5 (2023-12-08)

* Added R functions for survival predictions, e.g. time-dependent Brier scores, integrated Brier score

# psbcSpeedUp 2.0.4 (2023-10-12)

* For Linux machines, installed GitHub v2.0.4 seems twice faster than CRAN v2.0.4 due to `omp_set_*` in drive.cpp to solve CRAN pre-test NOTE `CPU-elapse times ratio`
* Improved C++ code to speedup 10-20%
* Declared functions inside PSBC class
* Added configuration files for checking omp, compilers etc.

# psbcSpeedUp 2.0.3 (2023-10-02)

* First CRAN version
* Added adaptive jumping rule
* Fixed bug for q=0 to allowed variable selection for all covariates by default

# psbcSpeedUp 2.0.2 (2023-09-23) (GitHub Only)

* Optimised the partition of the time axis for gamma process prior

# psbcSpeedUp 2.0.1 (2023-09-19) (GitHub Only)

* Set armadillo random values and RNG

# psbcSpeedUp 2.0 (2023-09-16) (GitHub Only)

* First GitHub commit version

</div>