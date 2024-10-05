#' Get Info about the dataset
#'
#' @description Helper function to extract information needed from the dataset
#'
#' @param DATA_IRT Data for the IRT model
#' @param GRADE_COL Label of the column where to find the grades
#' @param EXAM_COL Label of the column where to find exams ids
#' @param STUDID_COL Label of the column where to find students ids.
#' @param TIME_COL Label of the column where to find exams times.
#' @param MAXTIME_COL Label of the column where to find the maximum time observed for each student.
#'
#' @export
#' @importFrom rlang .data
#' @importFrom dplyr select arrange mutate group_by starts_with ungroup across everything
#' @importFrom tidyr pivot_wider nest unnest replace_na
dataRestruct <- function(
    DATA_IRT,
    GRADE_COL = "xcat",
    EXAM_COL = "itemID",
    STUDID_COL = "subjectID",
    TIME_COL = "time",
    MAXTIME_COL = "maxTime"
    ){
  grades_labs <- sort(unique(DATA_IRT[,GRADE_COL]), decreasing = FALSE)
  exams_labs <- sort(unique(DATA_IRT[, EXAM_COL]))
  ids <- sort(unique(DATA_IRT[, STUDID_COL]))

  n <- nrow(DATA_IRT)
  n_exams <- length(exams_labs)
  n_grades <- length(grades_labs) #observed

  exam_grades_tab <-  table(DATA_IRT[, GRADE_COL], DATA_IRT[, EXAM_COL],  useNA = 'ifany')

  nestedIRT <- DATA_IRT |>
    select('subject_id' = STUDID_COL,
           'max_time' = MAXTIME_COL,
           'exam_id' = EXAM_COL,
           'grade' = GRADE_COL,
           'time' = TIME_COL) |>
    arrange(.data$exam_id) |>
    mutate(todo = TRUE, obsflag=!is.na(.data$exam_id)) |>
    pivot_wider(id_cols = c('subject_id', 'max_time'), names_from = 'exam_id', values_from = c('todo', 'grade', 'time', 'todo')) |>
    group_by(.data$subject_id) |>
    nest(todo = starts_with('todo'),
                grades = starts_with('grade'),
                times = starts_with('time')) |>
    arrange(.data$subject_id)

  todoMat <- nestedIRT |> unnest('todo') |> select(.data$subject_id, starts_with('todo')) |> arrange(.data$subject_id) |>
    mutate(across(.cols = everything(), ~replace_na(., FALSE))) |>
    ungroup(.data$subject_id) |> select(-.data$subject_id) |> as.matrix()

  gradesMat <- nestedIRT |> unnest('grades') |> select(.data$subject_id, starts_with('grade')) |> arrange(.data$subject_id) |>
    ungroup(.data$subject_id) |> select(-.data$subject_id) |> as.matrix()

  gradesMat_idx <- matrix(sapply(gradesMat, function(x) if(!is.na(x)){which(grades_labs==as.character(x))}else{x}), nrow(gradesMat), ncol(gradesMat))
  gradesMat <- gradesMat_idx

  timeMat <- nestedIRT |> unnest('times') |> select(.data$subject_id, starts_with('time')) |> arrange(.data$subject_id) |>
    ungroup(.data$subject_id) |> select(-.data$subject_id) |> as.matrix()

  colnames(todoMat) <- colnames(gradesMat) <- colnames(timeMat) <- exams_labs
  rownames(todoMat) <- rownames(gradesMat) <- rownames(timeMat) <- ids

  out <- list(
    "obs_ids" = ids,
    "exams_tab" = exam_grades_tab,
    "grades_labs"= grades_labs,
    "exams_labs" = exams_labs,
    "n_observations" = n,
    "n_exams" = n_exams,
    "n_grades" = n_grades,
    "todoMat" = todoMat,
    "gradesMat" = gradesMat,
    "timeMat" = timeMat,
    "nested_IRT" = nestedIRT
  )

  return(out)
}
