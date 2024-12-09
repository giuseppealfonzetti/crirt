
#' Compute Marginal Loglik
#'
#' @param FIT Output from [fit_EM] or [fit_BFGS]
#' @param GRID Grid of quadrature points to be used.
#' @param WEIGHTS Weights for quadrature points.
#' @param INIT If TRUE, it evaluates the marginal log likelihood at the initial parameter values
#'
#' @export
get_loglik <- function(FIT, GRID=NULL, WEIGHTS=NULL, INIT=FALSE){
  if(is.null(GRID)) GRID <- FIT$grid
  if(is.null(WEIGHTS)) WEIGHTS <- FIT$weights

  if(!(FIT$mod%in%c("full", "grtc"))){
    stop("Model not available. Provide fit object for `full` or `grtc` models.")
  }

  pars <- FIT$fit$par
  if(INIT){
    pars <- FIT$start_par
  }

  out <- NA
  if(FIT$mod=="full"){

    out <- CRGRTCM_GH(
      THETA=pars,
      EXAMS_GRADES = FIT$data$gradesMat,
      EXAMS_DAYS = FIT$data$timeMat,
      EXAMS_SET = FIT$data$todoMat,
      EXAMS_OBSFLAG = !is.na(FIT$data$timeMat),
      MAX_DAY=FIT$data$fulldata$max_time,
      OUTCOME=FIT$data$outcome,
      EXT_COVARIATES=FIT$data$X,
      YEAR_FIRST=FIT$data$first_year,
      YEAR_LAST=FIT$data$last_year,
      YEAR_LAST_EXAM=FIT$data$yle,
      GRID=GRID,
      WEIGHTS=WEIGHTS,
      YB=FIT$data$yb,
      N_GRADES=FIT$data$n_grades,
      N_EXAMS=FIT$data$n_exams,
      GRFLAG = FALSE,
      ROTGRID = TRUE
    )$ll
  }else if(FIT$mod=="grtc"){
    out <- GRTCM_GH(
      THETA = pars,
      EXAMS_GRADES = FIT$data$gradesMat,
      EXAMS_DAYS = FIT$data$timeMat,
      EXAMS_SET = FIT$data$todoMat,
      EXAMS_OBSFLAG = !is.na(FIT$data$timeMat),
      MAX_DAY=FIT$data$fulldata$max_time,
      GRID=GRID,
      WEIGHTS=WEIGHTS,
      N_GRADES = FIT$data$n_grades,
      N_EXAMS = FIT$data$n_exams,
      GRFLAG = FALSE,
      ROTGRID = TRUE
    )$ll
  }

  return(out)
}
