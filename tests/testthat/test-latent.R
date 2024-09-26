test_that("internal_lat() loglikelihood", {
  skip_if_not_installed("mvtnorm")
  sig <- 2
  rho <- .6
  S <- matrix(c(1,rho*sig,rho*sig,sig^2), 2, 2)
  L <- t(chol(S))
  pars <- c(L[2,1], L[2,2])

  set.seed(123)
  for (i in 1:5) {
    lat <- rnorm(2)
    ll <- internal_lat(lat[1], lat[2], pars)$ll
    Rval <- mvtnorm::dmvnorm(lat, sigma = S, log = TRUE)
    expect_true(abs(ll-Rval)<1e-6)
  }
})

test_that("internal_lat() gradient", {
  skip_if_not_installed("numDeriv")
  sig <- 2
  rho <- .6
  S <- matrix(c(1,rho*sig,rho*sig,sig^2), 2, 2)
  L <- t(chol(S))
  pars <- c(L[2,1], L[2,2])

  rfun <- function(PAR, LAT){
    internal_lat(LAT[1], LAT[2], PAR)$ll
  }
  set.seed(123)
  for (i in 1:5) {
    lat <- rnorm(2)
    grll <- internal_lat(lat[1], lat[2], pars)$gr
    Rval <- numDeriv::grad(func = rfun, x = pars, LAT = lat)
    expect_equal(grll, Rval)
  }
})
