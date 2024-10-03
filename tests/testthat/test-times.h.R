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

#### test times probabilities output value and gradient #####
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

      test_that("time val () gradient",{
        skip_if_not_installed("numDeriv")
        grcpp <- gr_pTimeExam(
          EXAM = exam-1,
          DAY = day,
          THETA = theta,
          N_GRADES = n_grades,
          N_EXAMS = n_exams,
          SPEED = speeds[speed_index],
          ABILITY = 0,
          CDFFLAG = FALSE,
          ROTATED = FALSE
        )
        Rval <- numDeriv::grad(func = FUNTIME, x = theta,
                               SPEED = speeds[speed_index],
                               CDFFLAG = FALSE,
                               LOGFLAG = FALSE)

        expect_equal(Rval, grcpp)
      })

      test_that(paste0("pTimeExam() cdf gradient, speed:", round(speeds[speed_index],2), ", day:", round(day,2), ", exam:", exam),
                {
                  skip_if_not_installed("numDeriv")
                  grcpp <- gr_pTimeExam(
                    EXAM = exam-1,
                    DAY = day,
                    THETA = theta,
                    N_GRADES = n_grades,
                    N_EXAMS = n_exams,
                    SPEED = speeds[speed_index],
                    CDFFLAG = TRUE,
                    ABILITY = 0,
                    ROTATED = FALSE
                  )
                  Rval <- numDeriv::grad(func = FUNTIME, x = theta,
                                         SPEED = speeds[speed_index],
                                         CDFFLAG = TRUE,
                                         LOGFLAG = FALSE)

                  expect_true(sum(abs(Rval- grcpp)>1e-4)==0)
                })

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
