% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{gr_pGrade}
\alias{gr_pGrade}
\title{Evaluate the gradient of the probability of getting a specific grade}
\usage{
gr_pGrade(GRADE, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY)
}
\arguments{
\item{GRADE}{Grade used as reference}

\item{EXAM}{Exam of interest. Integers from 0 to N_EXAMS -1.}

\item{THETA_IRT}{Portion of the parameter vector related to the IRT model}

\item{N_GRADES}{Number of grades modelled.}

\item{N_EXAMS}{Number of exams.}

\item{ABILITY}{Ability value.}
}
\value{
It returns the gradient of the probability of obtaining the grade \code{GRADE} on exam \code{EXAM}.
}
\description{
Evaluate the gradient of the probability of getting a specific grade
}
