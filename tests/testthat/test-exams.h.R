### gen params ####
n_grades <- 4L
n_exams <- 5L
labs_exams <- paste0('ECO0',1:n_exams)
labs_grades <- c('[18,22)', '[22,25)', '[25,28)', '[29,30L]')
theta <- rnorm(n_exams*(n_grades+3)+2)
parList <- parVec2List(
  THETA = theta,
  N_GRADES = n_grades,
  N_EXAMS = n_exams,
  LABS_EXAMS = labs_exams,
  LABS_GRADES = labs_grades)
mat <- parList$irt




#### test exam-specific likelihood ####
FUNEXAM <- function(x, OBSFLAG, ROTATE=FALSE){
  sp <- speed
  if(ROTATE){
    sp <- x[n_exams*(n_grades+3)+1]*ability + x[n_exams*(n_grades+3)+2]*speed
  }
  examLik(
    EXAM = exam-1,
    GRADE = grade,
    DAY = day,
    MAX_DAY = day,
    OBSFLAG = OBSFLAG,
    THETA_IRT = x,
    N_GRADES = n_grades,
    N_EXAMS = n_exams,
    ABILITY = ability,
    SPEED = sp,
    LOGFLAG = TRUE
  )
}



set.seed(333)
for (ability in rnorm(3,0,1)) {
  for (speed in rnorm(3,0,1)) {
    for (day in runif(3, 100, 1000)) {
      for (exam in 1:n_exams) {
        for (grade in 1:n_grades) {
          pG <- pGrade(GRADE = grade, EXAM = exam-1, THETA_IRT = theta, N_GRADES = n_grades, N_EXAMS = n_exams, ABILITY = ability)
          pT <- dlnorm(day,
                       mat[exam,n_grades+2]-speed,
                       1/mat[exam,n_grades+3])
          Rval <- pT*pG
          val <- examLik(
            EXAM = exam-1,
            GRADE = grade,
            DAY = day,
            MAX_DAY = day,
            OBSFLAG = T,
            THETA_IRT = theta,
            N_GRADES = n_grades,
            N_EXAMS = n_exams,
            ABILITY = ability,
            SPEED = speed,
            LOGFLAG = TRUE
          )
          test_that("examLik() log output observed exam", {
            skip_if(!is.finite(log(Rval)))
            expect_equal(val, log(Rval))
          })

          pG <- pGreaterGrades(GRADE = 1, EXAM = exam-1, THETA_IRT = theta, N_GRADES = n_grades, N_EXAMS = n_exams, ABILITY = ability)
          pT <- plnorm(day,
                       mat[exam, n_grades+2]-speed,
                       1/mat[exam, n_grades+3])
          Rval <- 1-pT*pG
          val <- examLik(
            EXAM = exam-1,
            GRADE = grade,
            DAY = day,
            MAX_DAY = day,
            OBSFLAG = F,
            THETA_IRT = theta,
            N_GRADES = n_grades,
            N_EXAMS = n_exams,
            ABILITY = ability,
            SPEED = speed,
            LOGFLAG = TRUE
          )
          test_that("examLik() log output not observed exam", {
            skip_if(!is.finite(log(Rval)))
            expect_equal(val, log(Rval))
          })

        }

      }
    }
  }
}


## gradient checks ######
n <- 10

set.seed(123)

#### gen params ####
n_grades <- 4L
n_exams <- 10L
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

theta_idx <- 1:length(theta)
parListIdx <- parVec2List(
  THETA = theta_idx,
  N_GRADES = n_grades,
  N_EXAMS = n_exams,
  LABS_EXAMS = labs_exams,
  LABS_GRADES = labs_grades)
parListIdx$irt



latMat <- matrix(c(rnorm((n-1)*2), -5,5),n,2, byrow = TRUE)

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
RFUN <- function(x, SPEED, ABILITY, ROTATE,  ...){
  internal_speed <- SPEED
  if(ROTATE){
    internal_speed <- x[n_exams*(n_grades+3)+1]*ABILITY + x[n_exams*(n_grades+3)+2]*SPEED
  }

  examLik(
    THETA_IRT = x,
    N_GRADES = n_grades,
    N_EXAMS = n_exams,
    ABILITY = ABILITY,
    SPEED = internal_speed,
    LOGFLAG = TRUE,
    ...
  )
}
# numDeriv::grad(func = RFUN, x = theta, ABILITY = 1, SPEED = 1, ROTATE = FALSE, EXAM = 0, GRADE = 1, DAY = 100, MAX_DAY = 100, OBSFLAG = TRUE)
#
# grcpp <- grl_examLik(
#   EXAM = exam-1,
#   GRADE = grade,
#   DAY = day,
#   MAX_DAY = day,
#   OBSFLAG = F,
#   THETA_IRT = theta,
#   N_GRADES = n_grades,
#   N_EXAMS = n_exams,
#   ABILITY = ability,
#   SPEED = theta[n_exams*(n_grades+3)+1]*ability + theta[n_exams*(n_grades+3)+2]*speed,
#   ROTATED = TRUE
# )



for (i in 1:n) {
  for (exam in 1:n_exams) {
    test_that(paste0("exam density not rotated lat, i:",i, ", e:",exam), {
      skip_if_not_installed("numDeriv")

      numGrad <- numDeriv::grad(func = RFUN, x = theta, MAX_DAY = max_day,
                                EXAM = exam-1, DAY = timeMat[i, exam], GRADE = gradesMat[i, exam],
                                SPEED = latMat[i,2], ABILITY = latMat[i,1], ROTATE = FALSE,
                                OBSFLAG = obsMat[i, exam])
      grcpp <- grl_examLik(
        EXAM = exam-1, DAY = timeMat[i, exam], GRADE = gradesMat[i, exam],
        MAX_DAY = max_day,
        SPEED = latMat[i,2], ABILITY = latMat[i,1],
        ROTATED = FALSE,
        THETA_IRT = theta,
        OBSFLAG = obsMat[i, exam],
        N_GRADES = n_grades,
        N_EXAMS = n_exams
      )
      expect_equal(grcpp, numGrad,  tolerance = 1e-5)
    })

    test_that(paste0("exam density rotated lat, i:",i, ", e:",exam), {
      skip_if_not_installed("numDeriv")

      numGrad <- numDeriv::grad(func = RFUN, x = theta, MAX_DAY = max_day,
                                EXAM = exam-1, DAY = timeMat[i, exam], GRADE = gradesMat[i, exam],
                                SPEED = latMat[i,2], ABILITY = latMat[i,1], ROTATE = TRUE,
                                OBSFLAG = obsMat[i, exam])
      grcpp <- grl_examLik(
        EXAM = exam-1, DAY = timeMat[i, exam], GRADE = gradesMat[i, exam],
        MAX_DAY = max_day,
        SPEED = theta[n_exams*(n_grades+3)+1]*latMat[i,1] + theta[n_exams*(n_grades+3)+2]*latMat[i,2],
        ABILITY = latMat[i,1],
        ROTATED = TRUE,
        THETA_IRT = theta,
        OBSFLAG = obsMat[i, exam],
        N_GRADES = n_grades,
        N_EXAMS = n_exams
      )
      expect_equal(grcpp, numGrad,  tolerance = 1e-5)
    })


  }

}

# ROTFLAG <- TRUE
# i <- 1; exam <- 3
# numGrad <- numDeriv::grad(func = RFUN, x = theta, MAX_DAY = max_day,
#                           EXAM = exam-1, DAY = timeMat[i, exam], GRADE = gradesMat[i, exam],
#                           SPEED = latMat[i,2], ABILITY = latMat[i,1], ROTATE = ROTFLAG,
#                           OBSFLAG = obsMat[i, exam])
# grcpp <- grl_examLik(
#   EXAM = exam-1, DAY = timeMat[i, exam], GRADE = gradesMat[i, exam],
#   MAX_DAY = max_day,
#   SPEED = if_else(ROTFLAG,
#     theta[n_exams*(n_grades+3)+1]*latMat[i,1] + theta[n_exams*(n_grades+3)+2]*latMat[i,2],
#     latMat[i,2]),
#   ABILITY = latMat[i,1],
#   ROTATED = ROTFLAG,
#   THETA_IRT = theta,
#   OBSFLAG = obsMat[i, exam],
#   N_GRADES = n_grades,
#   N_EXAMS = n_exams
# )
# tibble("num" = numGrad, "cpp" = grcpp) |> print(n=100)
