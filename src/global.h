#ifndef GLOBAL
#define GLOBAL

// TODO: check if we omp is needed beyond drive.cpp
#ifdef _OPENMP
#include <omp.h>
omp_lock_t RNGlock; // GLOBAL NAMESPACE CAUSE I NEED IT NO MATTER WHERE
#endif

#include <RcppArmadillo.h>


#endif
