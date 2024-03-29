% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcPcoeff.R
\name{calcKp_Schmitt}
\alias{calcKp_Schmitt}
\title{Calculate partition coefficients for a molecule based on the Schmitt method}
\usage{
calcKp_Schmitt(logP, pKa, fup, type = 1, dat)
}
\arguments{
\item{logP}{Partition coefficient of a molecule between an aqueous and lipophilic phases, usually octanol and water; measurement of lipophilicity}

\item{pKa}{Negative log of the acid dissociation constant; measurement of the acidic strength of the molecule}

\item{fup}{Unbound fraction of the molecule in plasma}

\item{type}{Type of the molecule; 1=neutral, 2=monoprotic acid, 3=monoprotic base, 4=diprotic acid, 5=diprotic base, 6=monoprotic acid monoprotic base (acid comes first), 7=triprotic acid, 8=triprotic base, 9=diprotic acid monoprotic base (first two are acid), 10=diprotic base monoprotic acid (first one is acid)}

\item{dat}{Dataframe containing tissue composition data; columns are: tissue, f_water=water fraction, f_lipids=lipids fraction, f_proteins=protein fraction, f_pl=phospholipids fraction, f_n_l=neutral lipids fraction, f_n_pl=neutral phospholipids fraction, f_a_pl=acidic phospholipids fraction, pH, f_ew=extracellular fraction, f_iw=intracellular fraction, AR=albumin ratio, LR=lipoprotein ratio,  Prediction method; PT=Poulin and Theil, Berez=Berezhkovskiy, RR=Rodgers and Rowland, Schmitt=Schmitt, pksim=PK-Sim standard}
}
\value{
A named list with tissue:plasma partition coefficients
}
\description{
Takes in the molecule's physicochemical properties and returns the tissue:plasma partition coefficients based on the Schmitt method
}
\details{
Schmitt: https://pubmed.ncbi.nlm.nih.gov/17981004/
}
\keyword{internal}
