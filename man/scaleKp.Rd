% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helperFuns.R
\name{scaleKp}
\alias{scaleKp}
\title{Scale Kp values based on Vss}
\usage{
scaleKp(Kp, Vss, BP, Vt, Kpot, Ht = 0.45)
}
\arguments{
\item{Kp}{Named list of tissue:plasma partition coefficients}

\item{Vss}{Steady state volume in L}

\item{BP}{Blood:plasma concentration ratio}

\item{Vt}{Named list of tissue volumes generated by genInd}

\item{Kpot}{Partition coefficient for the "other" compartment}

\item{Ht}{Hematocrit value; default is 0.45}
}
\value{
A named list of the scaled tissue:plasma partition coefficients
}
\description{
Takes in Kp, Vss, Vt, and Kpot and returns Kp values that are scaled based on the input Vss value
}
\keyword{internal}
