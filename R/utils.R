#' @export
getDataInfo <- function(DATA_IRT, GRADE_COL = "xcat", EXAM_COL = "itemID", STUDID_COL = "subjectID"){
  grades_labs <- sort(unique(DATA_IRT[,GRADE_COL]), decreasing = FALSE)
  exams_labs <- unique(DATA_IRT[, EXAM_COL])
  ids <- unique(DATA_IRT[, STUDID_COL])

  n <- nrow(DATA_IRT)
  n_exams <- length(exams_labs)
  n_grades <- length(grades_labs) #observed

  exam_grades_tab <-  table(DATA_IRT[, GRADE_COL], DATA_IRT[, EXAM_COL],  useNA = 'ifany')

  out <- list(
    "exams_tab" = exam_grades_tab,
    "grades_labs"= grades_labs,
    "exams_labs" = exams_labs,
    "n_observations" = n,
    "n_exams" = n_exams,
    "n_grades" = n_grades
  )

  return(out)
}
