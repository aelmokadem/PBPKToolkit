% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcPcoeff.R
\name{calcBP}
\alias{calcBP}
\title{Calculate blood to plasma concentration ratio}
\usage{
calcBP(logP = NULL, fup, type = "total", method = 1, Ht = 0.45)
}
\arguments{
\item{logP}{Partition coefficient of a molecule between an aqueous and lipophilic phases, usually octanol and water; measurement of lipophilicity}

\item{fup}{Unbound fraction of the molecule in plasma}

\item{type}{Type of molecule; can be total, acid, base, or neutral}

\item{method}{BP calculation method; 1=fup-dependent method, 2=logP-dependent method}

\item{Ht}{Hematocrit value; default is 0.45}
}
\value{
Blood to plasma concentration ratio
}
\description{
Takes in logP, fup, and method and returns the blood to plasma concentration ratio
}
\details{
Source: https://pubmed.ncbi.nlm.nih.gov/20549836/
}
