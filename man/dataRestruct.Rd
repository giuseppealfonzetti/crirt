% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{dataRestruct}
\alias{dataRestruct}
\title{Get Info about the dataset}
\usage{
dataRestruct(
  DATA_IRT,
  GRADE_COL = "xcat",
  EXAM_COL = "itemID",
  STUDID_COL = "subjectID",
  TIME_COL = "time",
  MAXTIME_COL = "maxTime"
)
}
\arguments{
\item{DATA_IRT}{Data for the IRT model}

\item{GRADE_COL}{Label of the column where to find the grades}

\item{EXAM_COL}{Label of the column where to find exams ids}

\item{STUDID_COL}{Label of the column where to find students ids.}

\item{TIME_COL}{Label of the column where to find exams times.}

\item{MAXTIME_COL}{Label of the column where to find the maximum time observed for each student.}
}
\description{
Helper function to extract information needed from the dataset
}
