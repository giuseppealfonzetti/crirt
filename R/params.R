

#' IRT parameters to matrix
#'
#' Reparametrise IRT parameters and rearrange them in a matrix
#'
#' @param THETA_IRT Parameter vector.
#' @param N_GRADES number of grades modeled.
#' @param N_EXAMS number of exams.
#' @param LABS_EXAMS optional label for exams
#' @param LABS_GRADES optional label for grades#'
#'
#' @export
irtVec2Mat <- function(THETA_IRT, N_GRADES, N_EXAMS, LABS_EXAMS=NULL, LABS_GRADES=NULL){

  # Check input type
  if(!is.numeric(THETA_IRT)) stop('`THETA` not accepted.')
  if(!is.integer(N_GRADES)) stop('`N_GRADES` not accepted.')
  if(!is.integer(N_EXAMS)) stop('`N_EXAMS` not accepted.')

  # Check dimension
  dim <- N_EXAMS * (N_GRADES+3)
  if(dim!=length(THETA_IRT)) stop('THETA_IRT must contain N_EXAMS * (NGRADES+3) parameters')

  out <- matrix(THETA_IRT, nrow = N_EXAMS, byrow = TRUE)

  out[,1:N_GRADES] <- t(apply(out[,1:N_GRADES], 1, reparThr, CON2UN = FALSE))
  out[, N_GRADES+3] <- exp(out[, N_GRADES+3])

  if(is.null(LABS_EXAMS)) LABS_EXAMS <- paste0('exam',1:N_EXAMS)
  rownames(out) <- LABS_EXAMS
  if(is.null(LABS_GRADES)) LABS_GRADES <- paste0('grade',1:N_GRADES)
  colnames(out) <- c(LABS_GRADES, 'slope', 'time_loc', 'time_sc')

  return(out)
}

#' IRT parameters to matrix
#'
#' From constrained parameter matrix to unconstrained parameter vector
#'
#' @param MAT Parameter matrix (as returned by [irtVec2Mat]).
#'
#' @returns It reverts the transformation applied by [irtVec2Mat]
#'
#' @export
irtMat2Vec <- function(MAT){
  n_grades <- ncol(MAT)-3
  out <- MAT
  out[,1:n_grades] <- t(apply(out[,1:n_grades], 1, reparThr, CON2UN = TRUE))
  out[, n_grades+3] <- log(out[, n_grades+3])

  return(as.vector(t(out)))
}
