n <- 100
set.seed(123)

### gen params ####
n_grades <- 4L
n_exams <- 10L
labs_exams <- paste0('ECO0',1:n_exams)
labs_grades <- c('[18,22)', '[22,25)', '[25,28)', '[29,30L]')
theta <- rnorm(n_exams*(n_grades+3)+2)
parList <- parVec2List(
  THETA = theta,
  N_GRADES = n_grades,
  N_EXAMS = n_exams,
  LABS_EXAMS = labs_exams,
  LABS_GRADES = labs_grades)
irtMat <- parList$irt

theta_idx <- 1:length(theta)
parListIdx <- parVec2List(
  THETA = theta_idx,
  N_GRADES = n_grades,
  N_EXAMS = n_exams,
  LABS_EXAMS = labs_exams,
  LABS_GRADES = labs_grades)
parListIdx$irt



latMat <- matrix(rnorm(n*2), n, 2)

#### sim grades ####
gradesMat <- matrix(0, n, n_exams)

for (i in 1:n) {
  for (e in 1:n_exams) {
    #linear predictor exams X grades
    linp <- irtMat[e, n_grades+1] * latMat[i,1] - irtMat[e, 1:n_grades]

    # probabilities of greater grades
    pgg <- exp(linp)/(1+exp(linp))

    # probabilities of grades
    pg <- c(pgg[1:(n_grades-1)] - pgg[2:n_grades], pgg[n_grades])
    gradesMat[i,e] <- which(rmultinom(n = 1, size=1, prob = c(1-sum(pg), pg))==1)-1

  }
}



#### sim times ####
set.seed(123)
timeMat <- matrix(0, n, n_exams)
for (i in 1:n) {
  for (e in 1:n_exams) {
    timeMat[i,e] <- exp(
      rnorm(1,
            mean = irtMat[e, n_grades + 2]-latMat[i,2],
            sd = 1/irtMat[e, n_grades + 3])
    )
  }
}
timeMat[gradesMat==0] <- NA

#### mat to-do ####
todoMat <- matrix(1, n, n_exams)

#### censoring ####
max_day <- max(timeMat, na.rm = TRUE)+10
timeMat[timeMat>max_day] <- NA
obsMat <- matrix(1, n, n_exams)
obsMat[is.na(timeMat)] <- 0
#### checks ####

rfun <- function(PAR, ID, ROTATE){
  conditional_igrtcm(THETA = PAR,
                 EXAMS_GRADES = gradesMat[ID,],
                 EXAMS_DAYS = timeMat[ID,],
                 EXAMS_SET = todoMat[ID,],
                 EXAMS_OBSFLAG = obsMat[ID,],
                 MAX_DAY = max_day,
                 N_GRADES = n_grades,
                 N_EXAMS = n_exams,
                 ABILITY = latMat[ID,1],
                 SPEED = latMat[ID,2],
                 ROTATE = ROTATE)$ll
}

for (i in 1:n) {
  test_that(paste0("Gradient with no latent rotation, i:",i), {
    val <- conditional_igrtcm(THETA = theta,
                              EXAMS_GRADES = gradesMat[i,],
                              EXAMS_DAYS = timeMat[i,],
                              EXAMS_SET = todoMat[i,],
                              EXAMS_OBSFLAG = obsMat[i,],
                              MAX_DAY = max_day,
                              N_GRADES = n_grades,
                              N_EXAMS = n_exams,
                              ABILITY = latMat[i,1],
                              SPEED = latMat[i,2],
                              ROTATE = FALSE)$gr
    Rval <- numDeriv::grad(func = rfun, x = theta, ID = i, ROTATE = FALSE)
    expect_equal(val, Rval, tolerance = 1e-5)
  })

  test_that(paste0("Gradient with latent rotation, i:",i), {
    val <- conditional_igrtcm(THETA = theta,
                              EXAMS_GRADES = gradesMat[i,],
                              EXAMS_DAYS = timeMat[i,],
                              EXAMS_SET = todoMat[i,],
                              EXAMS_OBSFLAG = obsMat[i,],
                              MAX_DAY = max_day,
                              N_GRADES = n_grades,
                              N_EXAMS = n_exams,
                              ABILITY = latMat[i,1],
                              SPEED = latMat[i,2],
                              ROTATE = TRUE)$gr
    Rval <- numDeriv::grad(func = rfun, x = theta, ID = i, ROTATE = TRUE)
    expect_equal(val, Rval, tolerance = 1e-5)
  })
}


# ##### manual ####
# i <- 2; ROTFLAG <- TRUE
# val <- conditional_igrtcm(THETA = theta,
#                           EXAMS_GRADES = gradesMat[i,],
#                           EXAMS_DAYS = timeMat[i,],
#                           EXAMS_SET = todoMat[i,],
#                           EXAMS_OBSFLAG = obsMat[i,],
#                           MAX_DAY = max_day,
#                           N_GRADES = n_grades,
#                           N_EXAMS = n_exams,
#                           ABILITY = latMat[i,1],
#                           SPEED = latMat[i,2],
#                           ROTATE = ROTFLAG)$gr
# Rval <- numDeriv::grad(func = rfun, x = theta, ID = i, ROTATE = ROTFLAG)
# # expect_equal(val, Rval, tolerance = 1e-5)
# #val-Rval
#
#
#
#
# gru <- conditional_igrtcm(THETA = theta,
#                    EXAMS_GRADES = gradesMat[i,],
#                    EXAMS_DAYS = timeMat[i,],
#                    EXAMS_SET = todoMat[i,],
#                    EXAMS_OBSFLAG = obsMat[i,],
#                    MAX_DAY = max_day,
#                    N_GRADES = n_grades,
#                    N_EXAMS = n_exams,
#                    ABILITY = latMat[i,1],
#                    SPEED = if_else(ROTFLAG,
#                                    theta[n_exams*(n_grades+3)+1]*latMat[i,1] + theta[n_exams*(n_grades+3)+2]*latMat[i,2],
#                                    latMat[i,2]),
#                    ROTATE = FALSE)$gr
# grr <- conditional_igrtcm(THETA = theta,
#                           EXAMS_GRADES = gradesMat[i,],
#                           EXAMS_DAYS = timeMat[i,],
#                           EXAMS_SET = todoMat[i,],
#                           EXAMS_OBSFLAG = obsMat[i,],
#                           MAX_DAY = max_day,
#                           N_GRADES = n_grades,
#                           N_EXAMS = n_exams,
#                           ABILITY = latMat[i,1],
#                           SPEED = latMat[i,2],
#                           ROTATE = ROTFLAG)$gr
# tibble("no_rot" = gru, "rot" = grr, 'Rval' = Rval) |> print(n=100)
#
#
#
# #### single exams ####
# i <- 2
# FUNEXAM <- function(x, SPEED_RAW, ABILITY_RAW, ROTATE=ROTFLAG, ...){
#   sp <- SPEED_RAW
#   if(ROTATE){
#     sp <- x[n_exams*(n_grades+3)+1]*ABILITY_RAW + x[n_exams*(n_grades+3)+2]*SPEED_RAW
#   }
#   examLik(
#     ...,
#     THETA_IRT = x,
#     N_GRADES = n_grades,
#     N_EXAMS = n_exams,
#     ABILITY = ABILITY_RAW,
#     SPEED = sp,
#     LOGFLAG = TRUE
#   )
# }
#
# cppsum <- rep(0, length(theta))
# for (exm in 1:n_exams) {
#   Rexm <- numDeriv::grad(func = FUNEXAM, x = theta,
#                          SPEED_RAW = latMat[i,2],
#                          ABILITY_RAW = latMat[i,1],
#                          EXAM = exm-1,
#                          GRADE = gradesMat[i, exm],
#                          DAY = timeMat[i, exm],
#                          MAX_DAY = max_day,
#                          OBSFLAG = obsMat[i, exm], ROTATE = ROTFLAG)
#   cppexm <- grl_examLik(
#     EXAM = exm-1,
#     GRADE = gradesMat[i, exm],
#     DAY = timeMat[i, exm],
#     MAX_DAY = max_day,
#     OBSFLAG = obsMat[i, exm],
#     THETA_IRT = theta,
#     N_GRADES = n_grades,
#     N_EXAMS = n_exams,
#     ABILITY = latMat[i,1],
#     SPEED = if_else(ROTFLAG,
#                     theta[n_exams*(n_grades+3)+1]*latMat[i,1] + theta[n_exams*(n_grades+3)+2]*latMat[i,2],
#                     latMat[i,2]),#,
#     ROTATED = ROTFLAG
#   )
#   test_that(paste0("gradient exam  i:",i, ", e:",exm), {
#
#
#   expect_equal(cppexm, Rexm)
#   })
#   cppsum <- cppsum + cppexm
# }
#
# tibble("no_rot" = gru, "rot" = grr, 'Rval' = Rval, "cppsum"=cppsum, "Rexm" = Rexm, "cppexm" = cppexm) |> print(n=100)
#
# exm <- 10
# FUNEXAM(x = theta,
# SPEED_RAW = latMat[i,2],
# ABILITY_RAW = latMat[i,1],
# EXAM = exm-1,
# GRADE = gradesMat[i, exm],
# DAY = timeMat[i, exm],
# MAX_DAY = max_day,
# OBSFLAG = obsMat[i, exm], ROTATE = ROTFLAG)
#
# Rexm <- numDeriv::grad(func = FUNEXAM, x = theta,
#                        SPEED_RAW = latMat[i,2],
#                        ABILITY_RAW = latMat[i,1],
#                        EXAM = exm-1,
#                        GRADE = gradesMat[i, exm],
#                        DAY = timeMat[i, exm],
#                        MAX_DAY = max_day,
#                        OBSFLAG = obsMat[i, exm], ROTATE = ROTFLAG)
# cppexm <- grl_examLik(
#   EXAM = exm-1,
#   GRADE = gradesMat[i, exm],
#   DAY = timeMat[i, exm],
#   MAX_DAY = max_day,
#   OBSFLAG = obsMat[i, exm],
#   THETA_IRT = theta,
#   N_GRADES = n_grades,
#   N_EXAMS = n_exams,
#   ABILITY = latMat[i,1],
#   SPEED = if_else(ROTFLAG,
#                   theta[n_exams*(n_grades+3)+1]*latMat[i,1] + theta[n_exams*(n_grades+3)+2]*latMat[i,2],
#                   latMat[i,2]),#,
#   ROTATED = ROTFLAG
# )
# expect_equal(cppexm, Rexm)
#
# cppexm
# tibble("no_rot" = gru, "rot" = grr, 'Rval' = Rval,"Rexm" = Rexm, "cppexm" = cppexm) |> print(n=100)
#
#
#
#
# FUNTIME <- function(x, SPEED, ABILITY, ROTATE,  ...){
#   internal_speed <- SPEED
#   if(ROTATE){
#     internal_speed <- x[n_exams*(n_grades+3)+1]*ABILITY + x[n_exams*(n_grades+3)+2]*SPEED
#   }
#   pTimeExam(
#     THETA = x,
#     N_GRADES = n_grades,
#     N_EXAMS = n_exams,
#     SPEED = internal_speed,
#     ...
#   )
# }
# FUNTIME(x = theta,
#         EXAM = exm-1,
#         DAY = timeMat[i, exm],
#         SPEED = latMat[i,2],
#         ABILITY = latMat[i,1],
#         CDFFLAG = FALSE,
#         LOGFLAG = FALSE,
#         ROTATE = TRUE
#         )
# pTimeExam(
#   EXAM = exm-1,
#   DAY = timeMat[i, exm],
#   THETA = theta,
#   N_GRADES = n_grades,
#   N_EXAMS = n_exams,
#   SPEED = if_else(ROTFLAG,
#                   theta[n_exams*(n_grades+3)+1]*latMat[i,1] + theta[n_exams*(n_grades+3)+2]*latMat[i,2],
#                   latMat[i,2]),
#   CDFFLAG = FALSE,
#   LOGFLAG = TRUE
# )
# gr_pTimeExam(
#   EXAM = exm-1,
#   DAY = timeMat[i, exm],
#   THETA = theta,
#   N_GRADES = n_grades,
#   N_EXAMS = n_exams,
#   SPEED = if_else(ROTFLAG,
#                   theta[n_exams*(n_grades+3)+1]*latMat[i,1] + theta[n_exams*(n_grades+3)+2]*latMat[i,2],
#                   latMat[i,2]),
#   CDFFLAG = FALSE,
#   ABILITY = latMat[i,1],
#   ROTATED = ROTFLAG
# )
