#' Fit GRTCM
#'
#' Fit Graded Response Time Censored Model
#'
#' @param DATA Dataset in long form
#' @param GRADE_COL Label of the column where to find the grades
#' @param EXAM_COL Label of the column where to find exams ids
#' @param STUDID_COL Label of the column where to find students ids.
#' @param TIME_COL Label of the column where to find exams times.
#' @param MAXTIME_COL Label of the column where to find the maximum time observed for each student.
#' @param GRID Grid of quadrature points to be used.
#' @param WEIGHTS Weights for quadrature points
#' @param THETA_START Optional starting parameter vector
#'
#' @importFrom stats sd
#' @export
fit_GRTCM <- function(DATA, GRADE_COL = "xcat",
                      EXAM_COL = "itemID",
                      STUDID_COL = "subjectID",
                      TIME_COL = "time",
                      MAXTIME_COL = "maxTime",
                      GRID, WEIGHTS, THETA_START = NULL){

  dataStruct <- dataRestruct(DATA, GRADE_COL = GRADE_COL,
                             EXAM_COL = EXAM_COL,
                             STUDID_COL = STUDID_COL,
                             TIME_COL = TIME_COL,
                             MAXTIME_COL = MAXTIME_COL)

  if(is.null(THETA_START)){
    startIRTMat <- matrix(NA, dataStruct$n_exams, dataStruct$n_grades+3)
    startIRTMat[,1:dataStruct$n_grades] <- matrix(rep(seq(-2,4, length.out = dataStruct$n_grades), dataStruct$n_exams), dataStruct$n_exams, dataStruct$n_grades, byrow = TRUE)
    startIRTMat[,dataStruct$n_grades+1] <- 1
    startIRTMat[,dataStruct$n_grades+2] <- apply(dataStruct$timeMat, MARGIN = 2, FUN = function(x) mean(log(x), na.rm = TRUE))
    startIRTMat[,dataStruct$n_grades+3] <- apply(dataStruct$timeMat, MARGIN = 2, FUN = function(x) 1/sd(x, na.rm = TRUE))
    startLatMat <- diag(1, 2, 2)
    start_par <- parList2Vec(list("irt"=startIRTMat, 'lat_var'=startLatMat))
  }else{
    start_par <- THETA_START
  }

  NLL <- function(x){
    -GRTCM_GH(
      THETA = x,
      EXAMS_GRADES = dataStruct$gradesMat,
      EXAMS_DAYS = dataStruct$timeMat,
      EXAMS_SET = dataStruct$todoMat,
      EXAMS_OBSFLAG = !is.na(dataStruct$timeMat),
      MAX_DAY = dataStruct$nested_IRT$max_time,
      GRID = GRID,
      WEIGHTS = WEIGHTS,
      N_GRADES = dataStruct$n_grades,
      N_EXAMS = dataStruct$n_exams,
      GRFLAG = FALSE,
      ROTGRID = TRUE
    )$ll
  }

  NGR <- function(x){
    -GRTCM_GH(
      THETA = x,
      EXAMS_GRADES = dataStruct$gradesMat,
      EXAMS_DAYS = dataStruct$timeMat,
      EXAMS_SET = dataStruct$todoMat,
      EXAMS_OBSFLAG = !is.na(dataStruct$timeMat),
      MAX_DAY = dataStruct$nested_IRT$max_time,
      GRID = GRID,
      WEIGHTS = WEIGHTS,
      N_GRADES = dataStruct$n_grades,
      N_EXAMS = dataStruct$n_exams,
      GRFLAG = TRUE,
      ROTGRID = TRUE
    )$gr
  }

  fit <- ucminf::ucminf(par = start_par, fn = NLL, gr = NGR, hessian = 2)

  return(
    list(
      "data" = dataStruct,
      "start_par" = start_par,
      "fit" = fit,
      "parList" = parVec2List(
        THETA = fit$par,
        N_GRADES = dataStruct$n_grades,
        N_EXAMS = dataStruct$n_exams,
        LABS_EXAMS = dataStruct$exams_labs,
        LABS_GRADES = dataStruct$grades_labs)
    )
  )
}

#' MAP GRTCM
#'
#' Compute MAP scores for Graded Response Time Censored Model
#'
#' @param FIT Output from [fit_GRTCM]
#' @param TIDY Return tidy parameter table
#'
#' @importFrom dplyr as_tibble
#' @export
map_GRTCM <- function(FIT, TIDY = TRUE){



  ID_COMPLETENLL <- function(LAT, ID){
    -complete_igrtcm(
      THETA = FIT$fit$par,
      EXAMS_GRADES = FIT$data$gradesMat[ID,],
      EXAMS_DAYS = FIT$data$timeMat[ID,],
      EXAMS_SET = FIT$data$todoMat[ID,],
      EXAMS_OBSFLAG = !is.na(FIT$data$timeMat[ID,]),
      MAX_DAY = FIT$data$nested_IRT$max_time[ID],
      N_GRADES = FIT$data$n_grades,
      N_EXAMS = FIT$data$n_exams,
      ABILITY = LAT[1],
      SPEED = LAT[2])
  }

  mat <- Reduce(rbind,
                lapply(1:nrow(FIT$data$gradesMat),
                       function(ID){
                         fit <- ucminf::ucminf(par = c(0,0), fn = ID_COMPLETENLL, ID = ID)
                         as.numeric(fit$par)}))
  rownames(mat) <- FIT$data$obs_ids
  colnames(mat) <- c("ability", "speed")
  if(TIDY){
    mat <- as_tibble(mat, rownames = "subject_id")
  }

  return(mat)
}



#' Compute standard errors from GRTCM fit
#'
#' @param FIT Output from [fit_GRTCM]
#' @param TIDY Return tidy parameter table
#'
#' @importFrom rlang .data
#' @importFrom dplyr mutate
#'
#' @export
se_GRTCM <- function(FIT, TIDY = TRUE){

  reparJacob <- numDeriv::jacobian(func = parVec2Repar, x = FIT$fit$par,
                                   N_GRADES = FIT$data$n_grades, N_EXAMS = FIT$data$n_exams,
                                   LABS_EXAMS = FIT$data$exams_labs, LABS_GRADES = FIT$data$grades_labs)

  seVec <- sqrt(diag(t(reparJacob) %*% FIT$fit$invhessian %*% reparJacob))

  if(TIDY){
    out <- parVec2Repar(FIT$fit$par,
                           N_GRADES = FIT$data$n_grades, N_EXAMS = FIT$data$n_exams,
                           LABS_EXAMS = FIT$data$exams_labs, LABS_GRADES = FIT$data$grades_labs, TIDY = TRUE) |>
      mutate(
        se  = seVec,
        lb = .data$par-1.96*.data$se,
        ub = .data$par+1.96*.data$se,
        sig = !(.data$lb<0&.data$ub>0)
      )
  }else{
    out <- seVec
  }

  return(out)

}
