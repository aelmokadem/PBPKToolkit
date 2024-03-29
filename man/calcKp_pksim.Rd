% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcPcoeff.R
\name{calcKp_pksim}
\alias{calcKp_pksim}
\title{Calculate partition coefficients for a molecule based on the PK_Sim method}
\usage{
calcKp_pksim(logP, fup, dat)
}
\arguments{
\item{logP}{Partition coefficient of a molecule between an aqueous and lipophilic phases, usually octanol and water; measurement of lipophilicity}

\item{fup}{Unbound fraction of the molecule in plasma}

\item{dat}{Dataframe containing tissue composition data; columns are: tissue, f_water=water fraction, f_lipids=lipids fraction, f_proteins=protein fraction, f_pl=phospholipids fraction, f_n_l=neutral lipids fraction, f_n_pl=neutral phospholipids fraction, f_a_pl=acidic phospholipids fraction, pH, f_ew=extracellular fraction, f_iw=intracellular fraction, AR=albumin ratio, LR=lipoprotein ratio,  Prediction method; PT=Poulin and Theil, Berez=Berezhkovskiy, RR=Rodgers and Rowland, Schmitt=Schmitt, pksim=PK-Sim standard}
}
\value{
A named list with tissue:plasma partition coefficients
}
\description{
Takes in the molecule's physicochemical properties and returns the tissue:plasma partition coefficients based on the PK-Sim method
}
\details{
pksim: https://www.tandfonline.com/doi/abs/10.1517/17425255.1.1.159
}
\keyword{internal}
