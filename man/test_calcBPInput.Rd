% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helperFuns.R
\name{test_calcBPInput}
\alias{test_calcBPInput}
\title{Test if inputs to calcBP are compliant}
\usage{
test_calcBPInput(logP, fup, method)
}
\arguments{
\item{logP}{Partition coefficient of a molecule between an aqueous and lipophilic phases, usually octanol and water; measurement of lipophilicity}

\item{method}{BP calculation method; 1=logP-dependent method, 2=fup-dependent method}
}
\value{
An error if inputs to calcBP are not compliant
}
\description{
Takes in the logP and method for calcBP and returns an error if inputs to calcBP are not compliant
}
\keyword{internal}
