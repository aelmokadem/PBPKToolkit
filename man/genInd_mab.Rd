% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genPhys.R
\name{genInd_mab}
\alias{genInd_mab}
\title{Generate individual physiological parameters for the mAb model}
\usage{
genInd_mab(age, is.male, bw_targ = NULL, ht_targ = NULL, bmi_targ = NULL)
}
\arguments{
\item{age}{Age of the individual in years}

\item{is.male}{if TRUE, individual is male}

\item{bw_targ}{Target body weight in kg}

\item{ht_targ}{Target height in m}

\item{bmi_targ}{Target body mass index in kg/m^2}
}
\value{
Named list with physiological parameters for the mAb model for the desired individual
}
\description{
Takes in desired individual demographics and generates the individual physiological parameters for the mAb model
}
