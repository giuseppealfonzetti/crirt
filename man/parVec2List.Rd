% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/params.R
\name{parVec2List}
\alias{parVec2List}
\title{Unconstrained parameter vector to list}
\usage{
parVec2List(
  THETA,
  N_GRADES,
  N_EXAMS,
  N_COV,
  YB = 5,
  LABS_EXAMS = NULL,
  LABS_GRADES = NULL,
  LABS_COV = NULL
)
}
\arguments{
\item{THETA}{Parameter vector.}

\item{N_GRADES}{number of grades modeled.}

\item{N_EXAMS}{number of exams.}

\item{N_COV}{number of covariates in the competing risk model}

\item{YB}{Maximum number of years before graduation}

\item{LABS_EXAMS}{optional label for exams}

\item{LABS_GRADES}{optional label for grades}

\item{LABS_COV}{optional labels for covariates in the competing risk model}
}
\description{
Reparametrise unconstrained parameters and rearrange them into a list
}
