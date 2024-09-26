#include <Rcpp.h>
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#define EIGEN_DONT_PARALLELIZE
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
#include "extractParams.h"
#include "latent.h"
#include "irtMod.h"
