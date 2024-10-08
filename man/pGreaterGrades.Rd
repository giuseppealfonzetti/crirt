% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{pGreaterGrades}
\alias{pGreaterGrades}
\title{Evaluate the probability of grades greater or equal than the reference one}
\usage{
pGreaterGrades(
  GRADE,
  EXAM,
  THETA_IRT,
  N_GRADES,
  N_EXAMS,
  ABILITY,
  LOGFLAG = FALSE
)
}
\arguments{
\item{GRADE}{Grade used as reference. Integers from 1 to N_GRADES.}

\item{EXAM}{Exam of interest. Integers from 0 to N_EXAMS -1.}

\item{THETA_IRT}{Portion of the parameter vector related to the IRT model}

\item{N_GRADES}{Number of grades modelled.}

\item{N_EXAMS}{Number of exams.}

\item{ABILITY}{Ability value.}

\item{LOGFLAG}{Set TRUE to return log value.}
}
\value{
It returns the probability of obtaining grades higher than \code{GRADE} on exam \code{EXAM}.
}
\description{
Evaluate the probability of grades greater or equal than the reference one
}
