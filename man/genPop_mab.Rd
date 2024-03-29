% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genPhys.R
\name{genPop_mab}
\alias{genPop_mab}
\title{Generate population physiological parameters for the mAb model}
\usage{
genPop_mab(
  nSubj,
  minAge,
  maxAge,
  femPerc,
  minBW = NULL,
  maxBW = NULL,
  minHT = NULL,
  maxHT = NULL,
  minBMI = NULL,
  maxBMI = NULL
)
}
\arguments{
\item{nSubj}{Number of subjects}

\item{minAge}{Minimum age in years}

\item{maxAge}{Maximum age in years}

\item{femPerc}{Percentage of females}

\item{minBW}{Minimum body weight in kg}

\item{maxBW}{Maximum body weight in kg}

\item{minHT}{Minimum height in m}

\item{maxHT}{Maximum height in m}

\item{minBMI}{Minimum body mass index in kg/m^2}

\item{maxBMI}{Maximum body mass index in kg/m^2}
}
\value{
List of named lists with physiological parameters for the mAb model for each individual in the population
}
\description{
Takes in desired population demographics and generates the population physiological parameters for the mAb model
}
