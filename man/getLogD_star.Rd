% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helperFuns.R
\name{getLogD_star}
\alias{getLogD_star}
\title{Calculate logD_star}
\usage{
getLogD_star(type, logD, pKa, pH)
}
\arguments{
\item{type}{Type of the molecule; 1=neutral, 2=monoprotic acid, 3=monoprotic base, 4=diprotic acid, 5=diprotic base, 6=monoprotic acid monoprotic base (acid comes first), 7=triprotic acid, 8=triprotic base, 9=diprotic acid monoprotic base (first two are acid), 10=diprotic base monoprotic acid (first one is acid)}

\item{logD}{Olive oil:buffer partition coefficient of nonionized species}

\item{pKa}{Negative log of the acid dissociation constant; measurement of the acidic strength of the molecule}

\item{pH}{pH}
}
\value{
A named list with required components for partition coefficient calculation
}
\description{
Takes in the molecule type and returns the value for logD_star used by Poulin and Theil and Berezhkovskiy methods
}
\keyword{internal}
