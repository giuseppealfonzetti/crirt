n <- 10

### gen params IRT ####
n_grades <- 4L
n_exams <- 10L
labs_exams <- paste0('ECO0',1:n_exams)
labs_grades <- c('[18,22)', '[22,25)', '[25,28)', '[29,30L]')
theta_irt <- rnorm(n_exams*(n_grades+3)+2)


### gen params CR ####
set.seed(123)
m <- 1
yb <- 3
beta_int <- cbind(-(1:yb)+rnorm(yb), -2*(1:yb)+rnorm(yb))
beta_cov <- matrix(rnorm(m*2), m, 2)
beta_lat <- rbind(rnorm(2, -2), c(rnorm(1,1), rnorm(1,-3)))
beta <- rbind(beta_int, beta_cov, beta_lat)
beta
grad_par <- 3
theta_cr <- c(grad_par, as.numeric(beta))
theta <- c(theta_irt, theta_cr)
parList <- parVec2List(
  THETA = theta,
  N_GRADES = n_grades,
  N_EXAMS = n_exams,
  N_COV = m,
  YB = yb,
  LABS_EXAMS = labs_exams,
  LABS_GRADES = labs_grades)


### gen data
covariates <- matrix(rnorm(n), n, m)
qp <- matrix(c(-.5, .5), n , 2, byrow=TRUE)
U<-chol(parList$lat_var)









Rcrll <- function(PAR, OUTCOME, EXT_COVARIATES, LAT, YB, YEAR, YLE, ROTATE=FALSE){
  out <- cpp_crmod(
    OUTCOME=OUTCOME,
    YEAR_FIRST=1,
    YEAR_LAST=YEAR,
    THETA=PAR,
    EXT_COVARIATES=EXT_COVARIATES,
    YB=YB,
    YEAR_LAST_EXAM = YLE,
    LAT_POINTS=LAT,
    ROTATE=ROTATE
  )
  out$ll
}
Rcrgrll <- function(PAR, OUTCOME, EXT_COVARIATES, LAT, YB, YEAR, YLE, ROTATE=FALSE){
  out <- cpp_crmod(
    OUTCOME=OUTCOME,
    YEAR_FIRST=1,
    YEAR_LAST=YEAR,
    THETA=PAR,
    EXT_COVARIATES=EXT_COVARIATES,
    YB=YB,
    YEAR_LAST_EXAM = YLE,
    LAT_POINTS=LAT,
    ROTATE=ROTATE
  )
  out$grll
}

for (i in 1:n) {
  for (year in 1:yb) {
    for (outcome in 0:3) {
      test_that(paste0("Outcome loglik gradient obs ", i, ", year ", year, ", outcome ",outcome),{
        skip_if_not_installed("numDeriv")
        expect_equal(
          Rcrgrll(PAR = theta, OUTCOME = outcome, EXT_COVARIATES = covariates[i,], LAT=qp[i,],YB=yb, YEAR=year, YLE=yb-1, ROTATE=TRUE),
          numDeriv::grad(x=theta,
                         func = Rcrll,
                         OUTCOME = outcome, EXT_COVARIATES = covariates[i,], LAT=qp[i,],YB=yb, YEAR=year, YLE=yb-1, ROTATE=TRUE)
        )
      })
    }
  }
}
