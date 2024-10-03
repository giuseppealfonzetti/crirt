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

#### test grades probabilities output value and gradient ####
set.seed(123)
abilities <- sort(rnorm(3, 0, 2), decreasing = T)
for (abi_index in 1:length(abilities)) {
  for (grade in 1:n_grades) {
    for (exam in 1:n_exams) {
      probGG <- pGreaterGrades(
        GRADE = grade,
        EXAM = exam-1,
        THETA_IRT = theta,
        N_GRADES = n_grades,
        N_EXAMS = n_exams,
        ABILITY = abilities[abi_index]
      )
      lprobG <- pGrade(
        GRADE = grade,
        EXAM = exam-1,
        THETA_IRT = theta,
        N_GRADES = n_grades,
        N_EXAMS = n_exams,
        ABILITY = abilities[abi_index],
        LOGFLAG = TRUE
      )

      lprobGG <- pGreaterGrades(
        GRADE = grade,
        EXAM = exam-1,
        THETA_IRT = theta,
        N_GRADES = n_grades,
        N_EXAMS = n_exams,
        ABILITY = abilities[abi_index],
        LOGFLAG = TRUE
      )

      test_that("pGreaterGrades() val", {
        expect_equal(
          probGG,
          exp(mat[exam,n_grades+1]*abilities[abi_index]-mat[exam,grade])/(1+exp(mat[exam,n_grades+1]*abilities[abi_index]-mat[exam,grade]))
        )
      })

      test_that("pGreaterGrades() log val", {
        expect_equal(
          lprobGG,
          log(exp(mat[exam,n_grades+1]*abilities[abi_index]-mat[exam,grade])/(1+exp(mat[exam,n_grades+1]*abilities[abi_index]-mat[exam,grade])))
        )
      })

      if(grade>1){
        probprevGG <- pGreaterGrades(
          GRADE = grade-1,
          EXAM = exam-1,
          THETA_IRT = theta,
          N_GRADES = n_grades,
          N_EXAMS = n_exams,
          ABILITY = abilities[abi_index]
        )

        test_that("pGreaterGrades() decreases with higher grades",{
          expect_true(probprevGG>probGG)
        })



      }

      if(grade < n_grades){
        probnextGG <- pGreaterGrades(
          GRADE = grade+1,
          EXAM = exam-1,
          THETA_IRT = theta,
          N_GRADES = n_grades,
          N_EXAMS = n_exams,
          ABILITY = abilities[abi_index]
        )
        tmp <- log(probGG-probnextGG)
        test_that("pGrade() log",{
          expect_equal(lprobG, tmp)
        })
      }else{
        test_that("pGrade() log",{
          expect_equal(lprobG, lprobGG)
        })
      }

      if(abi_index>1){
        probprevGG <- pGreaterGrades(
          GRADE = grade,
          EXAM = exam-1,
          THETA_IRT = theta,
          N_GRADES = n_grades,
          N_EXAMS = n_exams,
          ABILITY = abilities[abi_index-1]
        )

        # check probs are decresing with higher grades
        test_that("pGreaterGrades() decreases with lower ability",{
          expect_true(probprevGG>probGG)
        })
      }






      FUNprobGG <-function(PAR){
        pGreaterGrades(
          GRADE = grade,
          EXAM = exam-1,
          THETA_IRT = PAR,
          N_GRADES = n_grades,
          N_EXAMS = n_exams,
          ABILITY = abilities[abi_index],
          LOGFLAG = FALSE
        )
      }
      grGG <- gr_pGreaterGrades(
        GRADE = grade,
        EXAM = exam-1,
        THETA_IRT = theta,
        N_GRADES = n_grades,
        N_EXAMS = n_exams,
        ABILITY = abilities[abi_index]
      )
      test_that("pGreaterGrades() gradient",{
        skip_if_not_installed("numDeriv")
        Rval <- numDeriv::grad(func = FUNprobGG, x = theta)

        expect_equal(Rval, grGG)
      })

      FUNprobG <-function(PAR){
        pGrade(
          GRADE = grade,
          EXAM = exam-1,
          THETA_IRT = PAR,
          N_GRADES = n_grades,
          N_EXAMS = n_exams,
          ABILITY = abilities[abi_index],
          LOGFLAG = FALSE
        )
      }
      grG <- gr_pGrade(
        GRADE = grade,
        EXAM = exam-1,
        THETA_IRT = theta,
        N_GRADES = n_grades,
        N_EXAMS = n_exams,
        ABILITY = abilities[abi_index]
      )
      test_that("pGrade() gradient",{
        skip_if_not_installed("numDeriv")
        Rval <- numDeriv::grad(func = FUNprobG, x = theta)

        expect_equal(Rval, grG)
      })
    }
  }
}

test_that("pGreaterGrades() log extreme ability", {

  abilities <- c(1000, -1000)
  for (abi_index in 1:length(abilities)) {
    for (grade in 1:n_grades) {
      for (exam in 1:n_exams) {
        prob <- pGreaterGrades(
          GRADE = grade,
          EXAM = exam-1,
          THETA_IRT = theta,
          N_GRADES = n_grades,
          N_EXAMS = n_exams,
          ABILITY = abilities[abi_index],
          LOGFLAG = TRUE
        )

        # check for finite values
        expect_true(is.finite(prob))

      }
    }
  }


})


test_that("check pGrade() probability space", {

  set.seed(321)
  abilities <- sort(rnorm(3, 0, 2), decreasing = T)
  for (abi_index in 1:length(abilities)) {
    for (exam in 1:n_exams) {
      probs <- rep(NA, n_grades+1)
      for (grade in 0:n_grades) {


        probs[grade+1] <- pGrade(
          GRADE = grade,
          EXAM = exam-1,
          THETA_IRT = theta,
          N_GRADES = n_grades,
          N_EXAMS = n_exams,
          ABILITY = abilities[abi_index]
        )




      }
      expect_equal(sum(probs),1)
    }

  }

})
