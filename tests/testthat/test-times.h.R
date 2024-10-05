n_grades <- 4L
n_exams <- 3L
labs_exams <- paste0('ECO0',1:n_exams)
labs_grades <- c('[18,22)', '[22,25)', '[25,28)', '[29,30L]')
set.seed(123)
theta <- rnorm(n_exams*(n_grades+3)+2)
parList <- parVec2List(
  THETA = theta,
  N_GRADES = n_grades,
  N_EXAMS = n_exams,
  LABS_EXAMS = labs_exams,
  LABS_GRADES = labs_grades)
mat <- parList$irt

#### test times probabilities output values #####
set.seed(111)
speeds <- sort(rnorm(3,0,2), decreasing=T)
FUNTIME <- function(x, CDFFLAG, LOGFLAG, SPEED){
  pTimeExam(
    EXAM = exam-1,
    DAY = day,
    THETA = x,
    N_GRADES = n_grades,
    N_EXAMS = n_exams,
    SPEED = SPEED,
    CDFFLAG = CDFFLAG,
    LOGFLAG = LOGFLAG
  )
}
for (speed_index in 1:length(speeds)) {
  for (day in runif(10, 100, 1000)) {
    for (exam in 1:n_exams) {



      Rval <- dlnorm(day,
                     mat[exam, n_grades+2]-speeds[speed_index],
                     1/mat[exam, n_grades+3])
      val <- FUNTIME(
        x = theta,
        CDFFLAG = FALSE,
        LOGFLAG = FALSE,
        SPEED = speeds[speed_index]
      )

      test_that(paste0("pTimeExam() val, speed:", round(speeds[speed_index],2), ", day:", round(day,2), ", exam:", exam), {
        expect_equal(val, Rval)
      })

      Rval <- plnorm(day,
                     mat[exam, n_grades+2]-speeds[speed_index],
                     1/mat[exam, n_grades+3])
      val <- FUNTIME(
        x = theta,
        CDFFLAG = TRUE,
        LOGFLAG = FALSE,
        SPEED = speeds[speed_index]
      )

      test_that(paste0("pTimeExam() cdf, speed:", round(speeds[speed_index],2), ", day:", round(day,2), ", exam:", exam),
                {
                  expect_equal(val, Rval)
                })


      if(speed_index>1){
        val <- FUNTIME(
          x = theta,
          CDFFLAG = TRUE,
          LOGFLAG = FALSE,
          SPEED = speeds[speed_index]
        )
        valprev <- FUNTIME(
          x = theta,
          CDFFLAG = TRUE,
          LOGFLAG = FALSE,
          SPEED = speeds[speed_index-1]
        )
        test_that(paste0("pTimeExam() val decreases with decresing speeds, speed:", round(speeds[speed_index],2), ", day:", round(day,2), ", exam:", exam),
                  {
                    expect_true(valprev>val)
                  })
      }




    }
  }

}

speeds <- sort(c(1000, -1000), decreasing=T)

set.seed(222)
for (speed_index in 1:length(speeds)) {
  for (day in runif(3, 100, 1000)) {
    for (exam in 1:n_exams) {


      val <- pTimeExam(
        EXAM = exam-1,
        DAY = day,
        THETA = theta,
        N_GRADES = n_grades,
        N_EXAMS = n_exams,
        SPEED = speeds[speed_index],
        CDFFLAG = F,
        LOGFLAG = TRUE
      )

      test_that(paste0("pTimeExam() log extreme, speed:", round(speeds[speed_index],2), ", day:", round(day,2), ", exam:", exam),
                {
                  expect_true(is.finite(val))
                })

      val <- pTimeExam(
        EXAM = exam-1,
        DAY = day,
        THETA = theta,
        N_GRADES = n_grades,
        N_EXAMS = n_exams,
        SPEED = speeds[speed_index],
        CDFFLAG = T,
        LOGFLAG = TRUE
      )
      test_that(paste0("pTimeExam() cdf log extreme, speed:", round(speeds[speed_index],2), ", day:", round(day,2), ", exam:", exam),
                {
                  expect_true(is.finite(val))
                })
    }
  }

}

########## gradient tests ##########

n <- 10
n_grades <- 4L
n_exams <- 10L
labs_exams <- paste0('ECO0',1:n_exams)
labs_grades <- c('[18,22)', '[22,25)', '[25,28)', '[29,30L]')
set.seed(123)
theta <- rnorm(n_exams*(n_grades+3)+2)
parList <- parVec2List(
  THETA = theta,
  N_GRADES = n_grades,
  N_EXAMS = n_exams,
  LABS_EXAMS = labs_exams,
  LABS_GRADES = labs_grades)
irtMat <- parList$irt
latMat <- matrix(c(rnorm((n-1)*2), -5,5),n,2, byrow = TRUE)
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


#### checks #####
FUNTIME <- function(x, SPEED, ABILITY, ROTATE,  ...){
  internal_speed <- SPEED
  if(ROTATE){
    internal_speed <- x[n_exams*(n_grades+3)+1]*ABILITY + x[n_exams*(n_grades+3)+2]*SPEED
    }
  pTimeExam(
    THETA = x,
    N_GRADES = n_grades,
    N_EXAMS = n_exams,
    SPEED = internal_speed,
    ...
  )
}


for (i in 1:n) {
  for (exam in 1:n_exams) {
    test_that(paste0("gradient time density non rotated lat, i:",i, ", e:",exam), {
      skip_if_not_installed("numDeriv")
      skip_if(is.na(timeMat[i, exam]))

      numGrad <- numDeriv::grad(func = FUNTIME, x = theta,
                                EXAM = exam-1, DAY = timeMat[i, exam], SPEED = latMat[i,2],
                                ABILITY = latMat[i,1], CDFFLAG = FALSE, LOGFLAG = FALSE, ROTATE = FALSE)
      grcpp <- gr_pTimeExam(
        EXAM = exam-1,
        DAY = timeMat[i, exam],
        THETA = theta,
        N_GRADES = n_grades,
        N_EXAMS = n_exams,
        SPEED = latMat[i,2],
        CDFFLAG = FALSE,
        ABILITY = latMat[i,1],
        ROTATED = FALSE,
        LOGFLAG = FALSE
      )
      expect_equal(grcpp, numGrad,  tolerance = 1e-5)
    })
    test_that(paste0("gradient time log density non rotated lat, i:",i, ", e:",exam), {
      skip_if_not_installed("numDeriv")
      skip_if(is.na(timeMat[i, exam]))

      numGrad <- numDeriv::grad(func = FUNTIME, x = theta,
                                EXAM = exam-1, DAY = timeMat[i, exam], SPEED = latMat[i,2],
                                ABILITY = latMat[i,1], CDFFLAG = FALSE, LOGFLAG = TRUE, ROTATE = FALSE)
      grcpp <- gr_pTimeExam(
        EXAM = exam-1,
        DAY = timeMat[i, exam],
        THETA = theta,
        N_GRADES = n_grades,
        N_EXAMS = n_exams,
        SPEED = latMat[i,2],
        CDFFLAG = FALSE,
        ABILITY = latMat[i,1],
        ROTATED = FALSE,
        LOGFLAG = TRUE
      )
      expect_equal(grcpp, numGrad,  tolerance = 1e-5)
    })

    test_that(paste0("gradient time density rotated lat, i:",i, ", e:",exam), {
      skip_if_not_installed("numDeriv")
      skip_if(is.na(timeMat[i, exam]))

      numGrad <- numDeriv::grad(func = FUNTIME, x = theta,
                                EXAM = exam-1, DAY = timeMat[i, exam], SPEED = latMat[i,2],
                                ABILITY = latMat[i,1], CDFFLAG = FALSE, LOGFLAG = FALSE, ROTATE = TRUE)
      grcpp <- gr_pTimeExam(
        EXAM = exam-1,
        DAY = timeMat[i, exam],
        THETA = theta,
        N_GRADES = n_grades,
        N_EXAMS = n_exams,
        SPEED = theta[n_exams*(n_grades+3)+1]*latMat[i,1] + theta[n_exams*(n_grades+3)+2]*latMat[i,2],
        CDFFLAG = FALSE,
        ABILITY = latMat[i,1],
        ROTATED = TRUE,
        LOGFLAG = FALSE
      )
      expect_equal(grcpp, numGrad,  tolerance = 1e-5)

    })
    test_that(paste0("gradient time log density rotated lat, i:",i, ", e:",exam), {
      skip_if_not_installed("numDeriv")
      skip_if(is.na(timeMat[i, exam]))

      numGrad <- numDeriv::grad(func = FUNTIME, x = theta,
                                EXAM = exam-1, DAY = timeMat[i, exam], SPEED = latMat[i,2],
                                ABILITY = latMat[i,1], CDFFLAG = FALSE, LOGFLAG = TRUE, ROTATE = TRUE)
      grcpp <- gr_pTimeExam(
        EXAM = exam-1,
        DAY = timeMat[i, exam],
        THETA = theta,
        N_GRADES = n_grades,
        N_EXAMS = n_exams,
        SPEED = theta[n_exams*(n_grades+3)+1]*latMat[i,1] + theta[n_exams*(n_grades+3)+2]*latMat[i,2],
        CDFFLAG = FALSE,
        ABILITY = latMat[i,1],
        ROTATED = TRUE,
        LOGFLAG = TRUE
      )
      expect_equal(grcpp, numGrad,  tolerance = 1e-5)

    })

    test_that(paste0("gradient time cdf non rotated lat, i:",i, ", e:",exam), {
      skip_if_not_installed("numDeriv")
      skip_if(is.na(timeMat[i, exam]))

      numGrad <- numDeriv::grad(func = FUNTIME, x = theta,
                                EXAM = exam-1, DAY = timeMat[i, exam], SPEED = latMat[i,2],
                                ABILITY = latMat[i,1], CDFFLAG = TRUE, LOGFLAG = FALSE, ROTATE = FALSE)
      grcpp <- gr_pTimeExam(
        EXAM = exam-1,
        DAY = timeMat[i, exam],
        THETA = theta,
        N_GRADES = n_grades,
        N_EXAMS = n_exams,
        SPEED = latMat[i,2],
        CDFFLAG = TRUE,
        ABILITY = latMat[i,1],
        ROTATED = FALSE,
        LOGFLAG = FALSE
      )
      expect_equal(grcpp, numGrad,  tolerance = 1e-5)
    })

    test_that(paste0("gradient time cdf rotated lat, i:",i, ", e:",exam), {
      skip_if_not_installed("numDeriv")
      skip_if(is.na(timeMat[i, exam]))

      numGrad <- numDeriv::grad(func = FUNTIME, x = theta,
                                EXAM = exam-1, DAY = timeMat[i, exam], SPEED = latMat[i,2],
                                ABILITY = latMat[i,1], CDFFLAG = TRUE, LOGFLAG = FALSE, ROTATE = TRUE)
      grcpp <- gr_pTimeExam(
        EXAM = exam-1,
        DAY = timeMat[i, exam],
        THETA = theta,
        N_GRADES = n_grades,
        N_EXAMS = n_exams,
        SPEED = theta[n_exams*(n_grades+3)+1]*latMat[i,1] + theta[n_exams*(n_grades+3)+2]*latMat[i,2],
        CDFFLAG = TRUE,
        ABILITY = latMat[i,1],
        ROTATED = TRUE,
        LOGFLAG = FALSE
      )
      expect_equal(grcpp, numGrad,  tolerance = 1e-5)

    })
  }
}

# i <- 3; exam <- 1
# FUNTIME(x = theta,
#         EXAM = exam-1, DAY = timeMat[i, exam], SPEED = latMat[i,2],
#         ABILITY = latMat[i,1], CDFFLAG = TRUE, LOGFLAG = FALSE, ROTATE = TRUE)
# numGrad <- numDeriv::grad(func = FUNTIME, x = theta,
#                           EXAM = exam-1, DAY = timeMat[i, exam], SPEED = latMat[i,2],
#                           ABILITY = latMat[i,1], CDFFLAG = TRUE, LOGFLAG = FALSE, ROTATE = TRUE)
# grcpp <- gr_pTimeExam(
#   EXAM = exam-1,
#   DAY = timeMat[i, exam],
#   THETA = theta,
#   N_GRADES = n_grades,
#   N_EXAMS = n_exams,
#   SPEED = theta[n_exams*(n_grades+3)+1]*latMat[i,1] + theta[n_exams*(n_grades+3)+2]*latMat[i,2],
#   CDFFLAG = TRUE,
#   ABILITY = latMat[i,1],
#   ROTATED = TRUE
# )
# expect_equal(numGrad, grcpp, tolerance = 1e-5)

