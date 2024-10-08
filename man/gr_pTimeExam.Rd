% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{gr_pTimeExam}
\alias{gr_pTimeExam}
\title{Evaluate the gradient of the c.d.f or p.d.f of the last attempt to an exam}
\usage{
gr_pTimeExam(
  EXAM,
  DAY,
  THETA,
  N_GRADES,
  N_EXAMS,
  SPEED,
  ABILITY,
  CDFFLAG,
  ROTATED,
  LOGFLAG = FALSE
)
}
\arguments{
\item{EXAM}{Exam of interest. Integers from 0 to N_EXAMS -1.}

\item{DAY}{Day of interest.}

\item{THETA}{Portion of the parameter vector related to the IRT model}

\item{N_GRADES}{Number of grades modelled.}

\item{N_EXAMS}{Number of exams.}

\item{SPEED}{speed value.}

\item{ABILITY}{ability value.}

\item{CDFFLAG}{\code{TRUE} for c.d.f. of time. \code{FALSE} for p.d.f.}

\item{ROTATED}{Have latent scores been rotated using their variance?}

\item{LOGFLAG}{TRUE to compute the gradient of the log density.(not available for cdf)}
}
\description{
Evaluate the gradient of the c.d.f or p.d.f of the last attempt to an exam
}
