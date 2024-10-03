n_grades <- 4L
n_exams <- 3L
labs_exams <- paste0('ECO0',1:n_exams)
labs_grades <- c('[18,22)', '[22,25)', '[25,28)', '[29,30L]')
set.seed(123)
theta_irt <- rnorm(n_exams*(n_grades+3))
mat <- irtVec2Mat(
  THETA_IRT = theta_irt,
  N_GRADES = n_grades,
  N_EXAMS = n_exams,
  LABS_EXAMS = labs_exams,
  LABS_GRADES = labs_grades
)




#### test exam-specific likelihood ####
FUNEXAM <- function(x, OBSFLAG){
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
    SPEED = speed,
    LOGFLAG = TRUE
  )
}



set.seed(333)
for (ability in rnorm(2,0,1)) {
  for (speed in rnorm(2,0,1)) {
    for (day in runif(2, 100, 1000)) {
      for (exam in 1:n_exams) {
        for (grade in 1:n_grades) {
          pG <- pGrade(GRADE = grade, EXAM = exam-1, THETA_IRT = theta_irt, N_GRADES = n_grades, N_EXAMS = n_exams, ABILITY = ability)
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
            THETA_IRT = theta_irt,
            N_GRADES = n_grades,
            N_EXAMS = n_exams,
            ABILITY = ability,
            SPEED = speed,
            LOGFLAG = TRUE
          )
          test_that("examLik() log output observed exam", {
            expect_equal(val, log(Rval))
          })

          pG <- pGreaterGrades(GRADE = 1, EXAM = exam-1, THETA_IRT = theta_irt, N_GRADES = n_grades, N_EXAMS = n_exams, ABILITY = ability)
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
            THETA_IRT = theta_irt,
            N_GRADES = n_grades,
            N_EXAMS = n_exams,
            ABILITY = ability,
            SPEED = speed,
            LOGFLAG = TRUE
          )
          test_that("examLik() log output observed exam", {
            expect_equal(val, log(Rval))
          })

          grcpp <- grl_examLik(
            EXAM = exam-1,
            GRADE = grade,
            DAY = day,
            MAX_DAY = day,
            OBSFLAG = F,
            THETA_IRT = theta_irt,
            N_GRADES = n_grades,
            N_EXAMS = n_exams,
            ABILITY = ability,
            SPEED = speed,
            ROTATED = FALSE
          )
          test_that("grl_examLik() ",{
            skip_if_not_installed("numDeriv")
            Rval <- numDeriv::grad(func = FUNEXAM, x = theta_irt, OBSFLAG = F)

            expect_true(sum(abs(Rval- grcpp)>1e-4)==0)
          })

        }

      }
    }
  }
}
