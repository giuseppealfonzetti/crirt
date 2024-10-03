n <- 10
set.seed(123)

### gen params ####
n_grades <- 4L
n_exams <- 3L
labs_exams <- paste0('ECO0',1:n_exams)
labs_grades <- c('[18,22)', '[22,25)', '[25,28)', '[29,30L]')
theta <- rnorm(n_exams*(n_grades+3)+2)
parList <- parVec2List(
  THETA = theta,
  N_GRADES = n_grades,
  N_EXAMS = n_exams,
  LABS_EXAMS = labs_exams,
  LABS_GRADES = labs_grades)
irtMat <- parList$irt


latMat <- matrix(rnorm(n*2), n, 2)

#### sim grades ####
gradesMat <- matrix(0, n, n_exams)

for (i in 1:n) {
  for (e in 1:n_exams) {
    #linear predictor exams X grades
    linp <- irtMat[e, n_grades+1] * latMat[i,1] - irtMat[e, 1:n_grades]

    # probabilities of greater grades
    pgg <- exp(linp)/(1+exp(linp))

    # probabilities of grades
    pg <- c(pgg[1:(n_grades-1)] - pgg[2:n_grades], pgg[n_grades])
    gradesMat[i,e] <- which(rmultinom(n = 1, size=1, prob = c(1-sum(pg), pg))==1)-1

  }
}



#### sim times ####
set.seed(123)
timeMat <- matrix(0, n, n_exams)
for (i in 1:n) {
  for (e in 1:n_exams) {
    timeMat[i,e] <- exp(
      rnorm(1,
            mean = irtMat[e, n_grades + 2]-latMat[i,2],
            sd = 1/irtMat[e, n_grades + 3])
    )
  }
}
timeMat[gradesMat==0] <- NA

#### mat to-do ####
todoMat <- matrix(1, n, n_exams)

#### censoring ####
max_day <- max(timeMat, na.rm = TRUE)+10
timeMat[timeMat>max_day] <- NA
obsMat <- matrix(1, n, n_exams)
obsMat[is.na(timeMat)] <- 0
#### checks ####

rfun <- function(PAR, ID, ROTATE){
  conditional_igrtcm(THETA = PAR,
                 EXAMS_GRADES = gradesMat[ID,],
                 EXAMS_DAYS = timeMat[ID,],
                 EXAMS_SET = todoMat[ID,],
                 EXAMS_OBSFLAG = obsMat[ID,],
                 MAX_DAY = max_day,
                 N_GRADES = n_grades,
                 N_EXAMS = n_exams,
                 ABILITY = latMat[ID,1],
                 SPEED = latMat[ID,2],
                 ROTATE = ROTATE)$ll
}

for (i in 1:n) {
  test_that("Gradient with no latent rotation", {
    val <- conditional_igrtcm(THETA = theta,
                              EXAMS_GRADES = gradesMat[i,],
                              EXAMS_DAYS = timeMat[i,],
                              EXAMS_SET = todoMat[i,],
                              EXAMS_OBSFLAG = obsMat[i,],
                              MAX_DAY = max_day,
                              N_GRADES = n_grades,
                              N_EXAMS = n_exams,
                              ABILITY = latMat[i,1],
                              SPEED = latMat[i,2],
                              ROTATE = FALSE)$gr
    Rval <- numDeriv::grad(func = rfun, x = theta, ID = i, ROTATE = FALSE)
    expect_equal(val, Rval)
  })

  test_that("Gradient with latent rotation", {
    val <- conditional_igrtcm(THETA = theta,
                              EXAMS_GRADES = gradesMat[i,],
                              EXAMS_DAYS = timeMat[i,],
                              EXAMS_SET = todoMat[i,],
                              EXAMS_OBSFLAG = obsMat[i,],
                              MAX_DAY = max_day,
                              N_GRADES = n_grades,
                              N_EXAMS = n_exams,
                              ABILITY = latMat[i,1],
                              SPEED = latMat[i,2],
                              ROTATE = TRUE)$gr
    Rval <- numDeriv::grad(func = rfun, x = theta, ID = i, ROTATE = TRUE)
    expect_equal(val, Rval)
  })
}

