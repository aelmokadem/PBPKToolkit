% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcPcoeff.R
\name{calcVss}
\alias{calcVss}
\title{Calculate steady state volume of distribution}
\usage{
calcVss(Kp, BP, Vt, Kpot = NULL, Ht = 0.45)
}
\arguments{
\item{Kp}{Named list of partition coefficients calculated by calcKp}

\item{BP}{Blood:plasma concentration ratio}

\item{Vt}{Named list of tissue volumes generated by genInd}

\item{Kpot}{Partition coefficient for the "other" compartment to be used when Vss value is provided. If Kpot is left NULL while Vss value is provided, Kpot will be calculated as the average of the non adipose partition coefficients}

\item{Ht}{Hematocrit value; default is 0.45}
}
\value{
Steady state volume of distribution
}
\description{
Takes in partition coefficients, BP, tissue volumes, Kpot, and hematocrit and returns steady state volume of distribution
}
\details{
Source: https://jcheminf.biomedcentral.com/articles/10.1186/s13321-015-0054-x
}
