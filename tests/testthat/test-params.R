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
test_that("irtVec2Mat() exams labs",{
  expect_equal(
    rownames(mat),
    labs_exams
  )
})
test_that("irtVec2Mat() grades labs",{
  expect_equal(
    colnames(mat)[1:n_grades],
    labs_grades
  )
})
for (exam in 1:nrow(mat)) {

  test_that("irtVec2Mat() intercepts",{
    expect_equal(
      extract_params_irt(theta_irt, n_grades, n_exams, 2, exam-1),
      as.numeric(mat[exam, 1:n_grades])
    )
  })

  test_that("irtVec2Mat() slope",{
    expect_equal(
      extract_params_irt(theta_irt, n_grades, n_exams, 1, exam-1),
      as.numeric(mat[exam, n_grades+1])
    )
  })

  test_that("irtVec2Mat() time loc",{
    expect_equal(
      extract_params_irt(theta_irt, n_grades, n_exams, 3, exam-1),
      as.numeric(mat[exam, n_grades+2])
    )})

  test_that("irtVec2Mat() time scale",{
    expect_equal(
      extract_params_irt(theta_irt, n_grades, n_exams, 4, exam-1),
      as.numeric(mat[exam, n_grades+3])
    )})
}
test_that("irtMat2Vec()", {
  expect_equal(theta_irt, irtMat2Vec(mat))
})
