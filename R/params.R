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
  colnames(out) <- c(LABS_GRADES, 'slope', 'time_loc', 'time_invsd')

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

#' Unconstrained parameter vector to list
#'
#' Reparametrise unconstrained parameters and rearrange them into a list
#'
#' @param THETA Parameter vector.
#' @param N_GRADES number of grades modeled.
#' @param N_EXAMS number of exams.
#' @param N_COV number of covariates in the competing risk model
#' @param LABS_EXAMS optional label for exams
#' @param LABS_GRADES optional label for grades
#' @param LABS_COV optional labels for covariates in the competing risk model
#' @param YB Maximum number of years before graduation
#'
#' @export
parVec2List <- function(THETA, N_GRADES, N_EXAMS, N_COV, YB=5, LABS_EXAMS=NULL, LABS_GRADES=NULL, LABS_COV=NULL){
  dim_irt <- N_EXAMS * (N_GRADES+3)

  if(is.null(LABS_COV)) LABS_COV <- paste0("X", 1:N_COV)
  out <- list()

  # IRT parameters
  out[["irt"]] <- irtVec2Mat(THETA_IRT = THETA[1:dim_irt],
                           N_GRADES = N_GRADES,
                           N_EXAMS = N_EXAMS,
                           LABS_EXAMS = LABS_EXAMS,
                           LABS_GRADES = LABS_GRADES)

  # Latent parameters
  L <- matrix(c(1,THETA[dim_irt+1], 0, THETA[dim_irt+2]),2,2)
  out[["lat_var"]] <- L %*% t(L)

  # Competing Risk parameters
  parCR <- THETA[(dim_irt+3):length(THETA)]
  out[["cr"]][["beta"]] <- matrix(parCR[-1], ncol=2)
  rownames(out$cr$beta) <- c(paste0("year", 1:YB), LABS_COV, "ability", "speed")
  colnames(out$cr$beta) <- c("dropout", "transfer")

  out[["cr"]][["grad"]] <- parCR[[1]]

  return(out)
}

#' Parameters list to unconstrained vector
#'
#' @param LIST List of parameters as returned by [parVec2List]
#'
#' @export
parList2Vec <- function(LIST){

  theta_irt <- irtMat2Vec(LIST[["irt"]])
  L <- t(chol(LIST[["lat_var"]]))
  theta_lat <- c(L[2,1], L[2,2])
  theta_cr <- c(LIST$cr$grad,as.numeric(LIST$cr$beta))

  theta <- c(theta_irt, theta_lat, theta_cr)

  return(theta)
}

#' Unconstrained parameter vector to reparametrisation
#'
#' Reparametrise unconstrained parameters
#'
#' @param THETA Parameter vector.
#' @param N_GRADES number of grades modeled.
#' @param N_EXAMS number of exams.
#' @param N_COV number of covariates in the competing risk model
#' @param LABS_EXAMS optional label for exams
#' @param LABS_GRADES optional label for grades
#' @param LABS_COV optional labels for covariates in the competing risk model
#' @param YB Maximum number of years before graduation
#' @param TIDY TRUE to get tidy parameters table
#'
#' @importFrom tidyr pivot_longer tribble as_tibble
#' @importFrom dplyr all_of bind_rows
#' @export
parVec2Repar <- function(THETA, N_GRADES, N_EXAMS, N_COV, YB, LABS_EXAMS=NULL, LABS_GRADES=NULL, LABS_COV=NULL, TIDY=FALSE){
  dim_irt <- N_EXAMS * (N_GRADES+3)
  if(is.null(LABS_COV)) LABS_COV <- paste0("X", 1:N_COV)
  out <- list()

  # IRT parameters
  irtMat <- irtVec2Mat(THETA_IRT = THETA[1:dim_irt],
                             N_GRADES = N_GRADES,
                             N_EXAMS = N_EXAMS,
                             LABS_EXAMS = LABS_EXAMS,
                             LABS_GRADES = LABS_GRADES)
  irtVec <- as.numeric(t(irtMat))

  # latent params
  L <- matrix(c(1,THETA[dim_irt+1], 0, THETA[dim_irt+2]),2,2)
  S <- L %*% t(L)
  speed_sd <- sqrt(S[2,2])
  lat_cor <- S[2,1]/speed_sd

  latVec <- c(speed_sd, lat_cor)


  if(TIDY){
    out <- as_tibble(irtMat, rownames = "group") |>
      pivot_longer(cols = -all_of("group"), names_to = "type", values_to = "par") |>
      bind_rows(
        tribble(
          ~group, ~type, ~par,
          "latent", "speed_sd", speed_sd,
          "latent", "correlation", lat_cor,
          "graduation", "intercept", THETA[(dim_irt+3)]
        )
      ) |>
      bind_rows(
        tibble(group=c(rep("dropout", YB+N_COV+2), rep("transfer", YB+N_COV+2)),
               type =c(paste0("year",1:YB), LABS_COV, "ability", "speed", paste0("year",1:YB), LABS_COV, "ability", "speed"),
               par=THETA[-(1:(dim_irt+3))])
      )
  }else{
    out <- c(irtVec, latVec, THETA[-(1:(dim_irt+2))])
  }

  return(out)
}
