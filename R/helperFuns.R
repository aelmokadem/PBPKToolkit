#####################################
## Helper functions for calcPcoeff ##
#####################################

#' Extract information from tissue composition data
#'
#' Takes in the tissue composition data and returns a named list of the required components for the partition coefficient calculation
#'
#' @param dat Dataframe containing tissue composition data; columns are: tissue, f_water=water fraction, f_lipids=lipids fraction, f_proteins=protein fraction, f_pl=phospholipids fraction, f_n_l=neutral lipids fraction, f_n_pl=neutral phospholipids fraction, f_a_pl=acidic phospholipids fraction, pH, f_ew=extracellular fraction, f_iw=intracellular fraction, AR=albumin ratio, LR=lipoprotein ratio,  Prediction method; PT=Poulin and Theil, Berez=Berezhkovskiy, RR=Rodgers and Rowland, Schmitt=Schmitt, pksim=PK-Sim standard
#' @return A named list with required components for partition coefficient calculation
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @keywords internal
getTissueData <- function(dat){

  Vwp <- dat$f_water[dat$tissue == "Plasma"]
  Vnlp <- dat$f_n_l[dat$tissue == "Plasma"]
  Vphp <- dat$f_pl[dat$tissue == "Plasma"]

  dat2 <- dat %>% filter(!tissue %in% c("Plasma","RBCs"))

  Vwt <- dat2$f_water[dat2$tissue != "Adipose"]
  Vwad <- dat2$f_water[dat2$tissue == "Adipose"]
  Vnlt <- dat2$f_n_l[dat2$tissue != "Adipose"]
  Vnlad <- dat2$f_n_l[dat2$tissue == "Adipose"]
  Vpht <- dat2$f_pl[dat2$tissue != "Adipose"]
  Vphad <- dat2$f_pl[dat2$tissue == "Adipose"]

  TC <- list(Vwp=Vwp, Vnlp=Vnlp, Vphp=Vphp, Vwt=Vwt, Vwad=Vwad, Vnlt=Vnlt, Vnlad=Vnlad, Vpht=Vpht, Vphad=Vphad)

  return(TC)
}

################

#' Calculate logD_star
#'
#' Takes in the molecule type and returns the value for logD_star used by Poulin and Theil and Berezhkovskiy methods
#'
#' @param type Type of the molecule; 1=neutral, 2=monoprotic acid, 3=monoprotic base, 4=diprotic acid, 5=diprotic base, 6=monoprotic acid monoprotic base (acid comes first), 7=triprotic acid, 8=triprotic base, 9=diprotic acid monoprotic base (first two are acid), 10=diprotic base monoprotic acid (first one is acid)
#' @param logD Olive oil:buffer partition coefficient of nonionized species
#' @param pKa Negative log of the acid dissociation constant; measurement of the acidic strength of the molecule
#' @param pH pH
#' @return A named list with required components for partition coefficient calculation
#' @keywords internal
getLogD_star <- function(type, logD, pKa, pH){
  logD_star <- switch(type,
                      #1-neutral
                      logD,
                      #2-monoprotic acid
                      logD-log10(1+10^(pH-pKa)),
                      #3-monoprotic base
                      logD-log10(1+10^(pKa-pH)),
                      #4-diprotic acid
                      logD-log10(1+10^(2*pH-pKa[1]-pKa[2])),
                      #5-diprotic base
                      logD-log10(1+10^(pKa[1]+pKa[2]-2*pH)),
                      #6-monoprotic acid monoprotic base (acid comes first)
                      logD-log10(1+10^(pKa[2]-pKa[1])),
                      #7-triprotic acid
                      logD-log10(1+10^(3*pH-pKa[1]-pKa[2]-pKa[3])),
                      #8-triprotic base
                      logD-log10(1+10^(pKa[1]+pKa[2]+pKa[3]-3*pH)),
                      #9-diprotic acid monoprotic base (first two are acid)
                      logD-log10(1+10^(pH-pKa[1]-pKa[2]+pKa[3])),
                      #10-diprotic base monoprotic acid (first one is acid)
                      logD-log10(1+10^(pKa[2]+pKa[3]-pKa[1]-pH)))

  return(logD_star)
}

################

#' Calculate X, Y, and Z
#'
#' Takes in the molecule type, pKa, and pH values for intracellular water, plasma, and RBCs and returns X, Y, and Z values required by the Rodgers and Rowland method
#'
#' @param type Type of the molecule; 1=neutral, 2=monoprotic acid, 3=monoprotic base, 4=diprotic acid, 5=diprotic base, 6=monoprotic acid monoprotic base (acid comes first), 7=triprotic acid, 8=triprotic base, 9=diprotic acid monoprotic base (first two are acid), 10=diprotic base monoprotic acid (first one is acid)
#' @param pKa Negative log of the acid dissociation constant; measurement of the acidic strength of the molecule
#' @param pH_IW pH of intracellular tissue water
#' @param pH_P pH of plasma
#' @param pH_RBC pH of red blood cells
#' @return A named list with X, Y, and Z values required by the Rodgers and Rowland method
#' @keywords internal
getXYZ <- function(type, pKa, pH_IW, pH_P, pH_RBC){
  X <- switch(type,
              #1-neutral
              0,
              #2-monoprotic acid
              10^(pH_IW-pKa),
              #3-monoprotic base
              10^(pKa-pH_IW),
              #4-diprotic acid
              10^(pH_IW-pKa[1])+10^(2*pH_IW-pKa[1]-pKa[2]),
              #5-diprotic base
              10^(pKa[2]-pH_IW)+10^(pKa[1]+pKa[2]-2*pH_IW),
              #6-monoprotic acid monoprotic base (acid comes first)
              10^(pKa[2]-pH_IW)+10^(pH_IW-pKa[1]),
              #7-triprotic acid
              10^(pH_IW-pKa[1])+10^(2*pH_IW-pKa[1]-pKa[2])+10^(3*pH_IW-pKa[1]-pKa[2]-pKa[3]),
              #8-triprotic base
              10^(pKa[3]-pH_IW)+10^(pKa[3]+pKa[2]-2*pH_IW)+10^(pKa[1]+pKa[2]+pKa[3]-3*pH_IW),
              #9-diprotic acid monoprotic base (first two are acid)
              10^(pKa[3]-pH_IW)+10^(pH_IW-pKa[1])+10^(2*pH_IW-pKa[1]-pKa[2]),
              #10-diprotic base monoprotic acid (first one is acid)
              10^(pH_IW-pKa[1])+10^(pKa[3]-pH_IW)+10^(pKa[2]+pKa[3]-2*pH_IW))

  Y <- switch(type,
              #1-neutral
              0,
              #2-monoprotic acid
              10^(pH_P-pKa),
              #3-monoprotic base
              10^(pKa-pH_P),
              #4-diprotic acid
              10^(pH_P-pKa[1])+10^(2*pH_P-pKa[1]-pKa[2]),
              #5-diprotic base
              10^(pKa[2]-pH_P)+10^(pKa[1]+pKa[2]-2*pH_P),
              #6-monoprotic acid monoprotic base (acid comes first)
              10^(pKa[2]-pH_P)+10^(pH_P-pKa[1]),
              #7-triprotic acid
              10^(pH_P-pKa[1])+10^(2*pH_P-pKa[1]-pKa[2])+10^(3*pH_P-pKa[1]-pKa[2]-pKa[3]),
              #8-triprotic base
              10^(pKa[3]-pH_P)+10^(pKa[3]+pKa[2]-2*pH_P)+10^(pKa[1]+pKa[2]+pKa[3]-3*pH_P),
              #9-diprotic acid monoprotic base (first two are acid)
              10^(pKa[3]-pH_P)+10^(pH_P-pKa[1])+10^(2*pH_P-pKa[1]-pKa[2]),
              #10-diprotic base monoprotic acid (first one is acid)
              10^(pH_P-pKa[1])+10^(pKa[3]-pH_P)+10^(pKa[2]+pKa[3]-2*pH_P))

  Z <- switch(type,
              #1-neutral
              1,
              #2-monoprotic acid
              1,
              #3-monoprotic base
              10^(pKa-pH_RBC),
              #4-diprotic acid
              1,
              #5-diprotic base
              10^(pKa[2]-pH_RBC)+10^(pKa[1]+pKa[2]-2*pH_RBC),
              #6-monoprotic acid monoprotic base (acid comes first)
              10^(pKa[2]-pH_RBC)+10^(pH_RBC-pKa[1]),
              #7-triprotic acid
              1,
              #8-triprotic base
              10^(pKa[3]-pH_RBC)+10^(pKa[3]+pKa[2]-2*pH_RBC)+10^(pKa[1]+pKa[2]+pKa[3]-3*pH_RBC),
              #9-diprotic acid monoprotic base (first two are acid)
              10^(pKa[3]-pH_RBC)+10^(pH_RBC-pKa[1])+10^(2*pH_RBC-pKa[1]-pKa[2]),
              #10-diprotic base monoprotic acid (first one is acid)
              10^(pH_RBC-pKa[1])+10^(pKa[3]-pH_RBC)+10^(pKa[2]+pKa[3]-2*pH_RBC))

  XYZ <- list(X=X, Y=Y, Z=Z)

  return(XYZ)
}

###############

#' Calculate W
#'
#' Takes in the molecule type, pKa, and pH values and returns W value required by the Schmitt method
#'
#' @param type Type of the molecule; 1=neutral, 2=monoprotic acid, 3=monoprotic base, 4=diprotic acid, 5=diprotic base, 6=monoprotic acid monoprotic base (acid comes first), 7=triprotic acid, 8=triprotic base, 9=diprotic acid monoprotic base (first two are acid), 10=diprotic base monoprotic acid (first one is acid)
#' @param pKa Negative log of the acid dissociation constant; measurement of the acidic strength of the molecule
#' @param pH_IW pH of intracellular tissue water
#' @param pH_P pH of plasma
#' @param pH_RBC pH of red blood cells
#' @return W value required by the Schmitt method
#' @keywords internal
getW <- function(type, pKa, pH){
  W <- switch(type,
              #1-neutral
              0,
              #2-monoprotic acid
              10^(pH-pKa),
              #3-monoprotic base
              10^(pKa-pH),
              #4-diprotic acid
              10^(pH-pKa[1])+10^(2*pH-pKa[1]-pKa[2]),
              #5-diprotic base
              10^(pKa[2]-pH)+10^(pKa[1]+pKa[2]-2*pH),
              #6-monoprotic acid monoprotic base (acid comes first)
              10^(pKa[2]-pH)+10^(pH-pKa[1]),
              #7-triprotic acid
              10^(pH-pKa[1])+10^(2*pH-pKa[1]-pKa[2])+10^(3*pH-pKa[1]-pKa[2]-pKa[3]))

  return(W)
}

####################

#' Calculate K_n_l and K_a_pl
#'
#' Takes in the molecule type, pKa, and pH values and returns W value required by the Schmitt method
#'
#' @param type Type of the molecule; 1=neutral, 2=monoprotic acid, 3=monoprotic base, 4=diprotic acid, 5=diprotic base, 6=monoprotic acid monoprotic base (acid comes first), 7=triprotic acid, 8=triprotic base, 9=diprotic acid monoprotic base (first two are acid), 10=diprotic base monoprotic acid (first one is acid)
#' @param pKa Negative log of the acid dissociation constant; measurement of the acidic strength of the molecule
#' @param pH pH
#' @param alpha Ratio between ditribution coefficient at given pH (D) and that in neutral form (D0)
#' @param W W value calculated by the getW function
#' @param K_n_pl Neutral phospholipids:water partition coefficient
#' @return A named list with K_n_l and K_a_pl (neutral lipids:water and acidic phospholipids:water partition coefficients, respectively)
#' @keywords internal
getKs <- function(type, pKa, pH, alpha, W, K_n_pl){
  if(type==1 | type==2 | type==4 | type==7){ # neutral, monoprotic acid, diprotic acid, triprotic acid
    K_n_l <- K_n_pl*(((1-alpha)/(1+W))+alpha)
    K_a_pl <- K_n_pl*((1/(1+W))+0.05*(1-(1/(1+W))))
  }
  else if(type==3){ # monoprotic base
    K_n_l <- K_n_pl*(((1-alpha)/(1+W))+alpha)
    K_a_pl <- K_n_pl*((1/(1+W))+20*(1-(1/(1+W))))
  }
  else if(type==5){ # diprotic base
    F1 <- (1/(1+10^(pKa[1]-pH)))
    F2 <- (1/(1+10^(pKa[2]-pH)))
    K_n_l <- K_n_pl*(F1*F2 + alpha*((1-F1)*F2 + F1*(1-F2)) + (1-F1)*(1-F2))
    K_a_pl <- K_n_pl*(F1*F2 + 20*((1-F1)*F2 + F1*(1-F2)) + (1-F1)*(1-F2))
  }
  else if(type==6){ # monoprotic acid monoprotic base (acid comes first)
    F1 <- (1/(1+10^(pH-pKa[1])))
    F2 <- (1/(1+10^(pKa[2]-pH)))
    K_n_l <- K_n_pl*(F1*F2 + alpha*((1-F1)*F2 + F1*(1-F2)) + (1-F1)*(1-F2))
    K_a_pl <- K_n_pl*(F1*F2 + 0.05*(1-F1)*F2 + 20*(1-(F2))*F1 + (1-F1)*(1-F2))
  }
  else if(type==8){ # triprotic base
    F1 <- (1/(1+10^(pKa[1]-pH)))
    F2 <- (1/(1+10^(pKa[2]-pH)))
    F3 <- (1/(1+10^(pKa[3]-pH)))
    K_n_l <- K_n_pl*(F1*F2*F3 + alpha*((1-F1)*F2*F3 + F1*(1-F2)*F3 + F1*F2*(1-F3) + (1-F1)*(1-F2)*F3 + (1-F1)*F2*(1-F3) + F1*(1-F2)*(1-F3) + (1-F1)*(1-F2)*(1-F3)))
    K_a_pl <- K_n_pl*(F1*F2*F3 + 20*((1-F1)*F2*F3 + F1*(1-F2)*F3 + F1*F2*(1-F3)) + (1-F1)*(1-F2)*F3 + (1-F1)*F2*(1-F3) + F1*(1-F2)*(1-F3) + (1-F1)*(1-F2)*(1-F3))
  }
  else if(type==9){ # diprotic acid monoprotic base (first two are acid)
    F1 <- (1/(1+10^(pH-pKa[1])))
    F2 <- (1/(1+10^(pH-pKa[2])))
    F3 <- (1/(1+10^(pKa[3]-pH)))
    K_n_l <- K_n_pl*(F1*F2*F3 + alpha*((1-F1)*F2*F3 + F1*(1-F2)*F3 + F1*F2*(1-F3) + (1-F1)*(1-F2)*F3 + (1-F1)*(F2)*(1-F3) + F1*(1-F2)*(1-F3) + (1-F1)*(1-F2)*(1-F3)))
    K_a_pl <- K_n_pl*(F1*F2*F3 + 0.05*((1-F1)*F2*F3 + F1*(1-F2)*F3 + (1-F1)*(1-F2)*F3) + 20*F1*F2*(1-F3) + (1-F1)*F2*(1-F3) + F1*(1-F2)*(1-F3) + (1-F1)*(1-F2)*(1-F3))
  }
  else if(type==10){ # diprotic base monoprotic acid (first one is acid)
    F1 <- (1/(1+10^(pH-pKa[1])))
    F2 <- (1/(1+10^(pKa[2]-pH)))
    F3 <- (1/(1+10^(pKa[3]-pH)))
    K_n_l <- K_n_pl*(F1*F2*F3 + alpha*((1-F1)*F2*F3 + F1*(1-F2)*F3 + F1*F2*(1-F3) + (1-F1)*(1-F2)*F3 + (1-F1)*(F2)*(1-F3) + F1*(1-F2)*(1-F3) + (1-F1)*(1-F2)*(1-F3)))
    K_a_pl <- K_n_pl*(F1*F2*F3 + 0.05*(1-F1)*F2*F3 + 20*(F1*(1-F2)*F3 + F1*F2*(1-F3) + F1*(1-F2)*(1-F3)) + (1-F1)*F2*(1-F3) + (1-F1)*(1-F2)*F3 + (1-F1)*(1-F2)*(1-F3))
  }

  Ks <- list(K_n_l=K_n_l, K_a_pl=K_a_pl)

  return(Ks)
}

##########################

#' Test matching pKa length and type
#'
#' Takes in the molecule type and pKa and throws an error if they don't match
#'
#' @param type Type of the molecule; 1=neutral, 2=monoprotic acid, 3=monoprotic base, 4=diprotic acid, 5=diprotic base, 6=monoprotic acid monoprotic base (acid comes first), 7=triprotic acid, 8=triprotic base, 9=diprotic acid monoprotic base (first two are acid), 10=diprotic base monoprotic acid (first one is acid)
#' @param pKa Negative log of the acid dissociation constant; measurement of the acidic strength of the molecule
#' @return An error if pKa length and type don't match
#' @keywords internal
test_pKaTypeMatch <- function(type, pKa){
  if(type %in% c(1) & !is.null(pKa)) warning("Molecule type 1 does not require pKa so it will be ignored")
  if(type %in% c(2,3) & length(pKa) != 1) stop("Molecule types 2 and 3 require one pKa value")
  if(type %in% c(4,5,6) & length(pKa) != 2) stop("Molecule types 4, 5, and 6 require two pKa values")
  if(type >= 7 & length(pKa) != 3) stop("Molecule types 7, 8, 9, and 10 require three pKa values")
}

####################

#' Test for negative Kp values
#'
#' Takes in a named list of Kp values and returns an error if a negative value occurs
#'
#' @param Kps Named list with Kp values
#' @return An error if one or more Kp values is < 0
#' @keywords internal
test_negativeKps <- function(Kps){
  if(any(as.numeric(Kps) < 0)) stop("One or more partition coefficients have values < 0. Make sure the input parameters are correct or use a different method")
}

####################

#' Test for missing Vt
#'
#' Takes in Vt and returns an error if Vt is missing
#'
#' @param Vt Named list of tissue volumes generated by genInd
#' @return An error if Vss is provided without Vt
#' @keywords internal
test_missingVt <- function(Vt){
  if(is.null(Vt)) stop("Vt needs to be provided when Vss is used to calculate Kps")
}

####################

#' Scale Kp values based on Vss
#'
#' Takes in Kp, Vss, and Vt and returns an error if Vss is provided without Vt
#'
#' @param Kp Named list of tissue:plasma partition coefficients
#' @param Vss Steady state volume
#' @param BP Blood:plasma concentration ratio
#' @param Vt Named list of tissue volumes generated by genInd
#' @return A named list of the scaled tissue:plasma partition coefficients
#' @keywords internal
scaleKp <- function(Kp, Vss, BP, Vt){
  nms <- names(Vt[grepl("V", names(Vt))])
  vols <- Vt[nms]
  guNms <- c("Vla_int","Vsm_int","Vst")
  otNms <- c("Vot","Vpa","Vln","Vth","Vgo")
  Vgu <- sum(as.numeric(vols[guNms]))
  Vot <- sum(as.numeric(vols[otNms]))

  Vbl <- vols[["Vve"]] + vols[["Var"]]
  f_plasma <- ifelse(Vt[["SEX"]] == 1, 3000/(2300 + 3000), 2400/(1500 + 2400))  #ICRP publication 89
  f_rbcs <- 1 - f_plasma
  Vp <- f_plasma * Vbl
  Ve <- f_rbcs * Vbl

  # calculate EP according to https://pubmed.ncbi.nlm.nih.gov/11782904/
  Ht <- 0.45
  EP <- (BP - (1 - Ht))/Ht

  vols2 <- vols[setdiff(names(vols), c(guNms, otNms, "Vve", "Var"))]
  vols2$Vgu <- Vgu

  names(vols2) <- gsub("V","",names(vols2))
  Kp2 <- Kp
  names(Kp2) <- gsub("Kp","",names(Kp2))
  nms_Kp <- names(Kp2)

  vols3 <- vols2[nms_Kp]

  # calculate a scale for Kp by reorganizing the equation in https://jcheminf.biomedcentral.com/track/pdf/10.1186/s13321-015-0054-x.pdf
  scale <- (Vss - (Vp + EP*Ve)) / (sum(as.numeric(Kp) * as.numeric(vols3)))

  Kp2 <- lapply(Kp, function(x) scale*x)

  return(Kp2)
}

#####################################
#####################################

######################################
#### Helper functions for genPhys ####
######################################

#' Test if the target covariates are out of range
#'
#' Takes in the target covariates and corresponding ranges for the chosen age and sex and returns an error if one or more are out of range
#'
#' @param bw_targ Target body weight in kg
#' @param ht_targ Target height in m
#' @param bmi_targ Target BMI in kg/m^2
#' @param rangeBW Range of body weights in kg for the chosen age and sex
#' @param rangeHT Range of heights in m for the chosen age and sex
#' @param rangeBW Range of BMIs in kg/m^2 for the chosen age and sex
#' @return An error message if one or more of the covariates are out of range
#' @keywords internal
test_covRange <- function(bw_targ, ht_targ, bmi_targ, rangeBW, rangeHT, rangeBMI){
  if((bw_targ < rangeBW[1] || bw_targ > rangeBW[2]) & (ht_targ < rangeHT[1] || ht_targ > rangeHT[2]) & (bmi_targ < rangeBMI[1] || bmi_targ > rangeBMI[2])) stop("Target body weight, height, and BMI are out of range for the chosen age and sex")
  if((bw_targ < rangeBW[1] || bw_targ > rangeBW[2]) & (ht_targ < rangeHT[1] || ht_targ > rangeHT[2])) stop("Target body weight and height are out of range for the chosen age and sex")
  if((bw_targ < rangeBW[1] || bw_targ > rangeBW[2]) & (bmi_targ < rangeBMI[1] || bmi_targ > rangeBMI[2])) stop("Target body weight and BMI are out of range for the chosen age and sex")
  if((ht_targ < rangeHT[1] || ht_targ > rangeHT[2]) & (bmi_targ < rangeBMI[1] || bmi_targ > rangeBMI[2])) stop("Target height and BMI are out of range for the chosen age and sex")
  if(bw_targ < rangeBW[1] || bw_targ > rangeBW[2]) stop("Target body weight is out of range for the chosen age and sex")
  if(ht_targ < rangeHT[1] || ht_targ > rangeHT[2]) stop("Target height is out of range for the chosen age and sex")
  if(bmi_targ < rangeBMI[1] || bmi_targ > rangeBMI[2]) stop("Target BMI is out of range for the chosen age and sex")
}

######################################

#' Test if inputs to genInd are compliant
#'
#' Takes in the target covariates and returns an error or a warning if inputs to genInd are not compliant
#'
#' @param bw_targ Target body weight in kg
#' @param ht_targ Target height in m
#' @param bmi_targ Target BMI in kg/m^2
#' @return An error or a warning if inputs to genInd are not compliant
#' @keywords internal
test_genIndInput <- function(bw_targ, ht_targ, bmi_targ){
  if((is.null(bw_targ) && is.null(ht_targ)) | (is.null(bw_targ) && is.null(bmi_targ)) | (is.null(ht_targ) && is.null(bmi_targ))) stop("At least two inputs of bw_targ, ht_targ, or bmi_targ are required")
  if(!is.null(bw_targ) && !is.null(ht_targ) && !is.null(bmi_targ)) warning("bmi_targ will be calculated from bw_targ and bmi_targ; input value will be ignored")
}

######################################

#' Refine datasets based on age and sex
#'
#' Takes in the target age and sex and returns a list of datasets refined based on these covariates
#'
#' @param age Age in years
#' @param sex Sex
#' @return A named list with the refined datasets based on age and sex
#' @keywords internal
filterDatasets <- function(age, sex){
  if(sex == 1){
    dat <- nhanesData %>% filter(AGE_YR==age, SEX==sex)  #filter male nhanes data
    df <- icrpData %>% select(grep("_m", names(icrpData)))  #filter male physiological data
    names(df) <- gsub("_m", "", names(df))  #remove "_m" from df names
    normSD[,"cv_f"] <- NULL  #remove female CV column
    names(normSD) <- gsub("_m", "", names(normSD))  #remove "_m" from df names
    flow[,"flowPerc_f"] <- NULL  #remove female flows column
    names(flow) <- gsub("_m", "", names(flow))  #remove "_m" from df names
    BC[,"bloodPerc_f"] <- NULL  #remove female blood content column
    names(BC) <- gsub("_m", "", names(BC))  #remove "_m" from df names
  }else{
    dat <- nhanesData %>% filter(AGE_YR==age, SEX==sex)  #filter female nhanes data
    df <- icrpData %>% select(grep("_f", names(icrpData)))  #filter female physiological data
    names(df) <- gsub("_f", "", names(df))  #remove "_f" from df names
    normSD[,"cv_m"] <- NULL  #remove female CV column
    names(normSD) <- gsub("_f", "", names(normSD))  #remove "_m" from df names
    flow[,"flowPerc_m"] <- NULL  #remove male flows column
    names(flow) <- gsub("_f", "", names(flow))  #remove "_m" from df names
    BC[,"bloodPerc_m"] <- NULL  #remove male blood content column
    names(BC) <- gsub("_f", "", names(BC))  #remove "_f" from df names
  }

  l <- list(dat=dat, df=df, normSD=normSD, flow=flow, BC=BC)
  return(l)
}

######################################

#' Optimize organ volumes
#'
#' Takes in the a vector of non-optimized organ volumes and returns a probability to be minimized by an optimizer
#'
#' @param optVols Vector of non-optimized organ volumes
#' @param df_opt Dataframe with volumes of organs to be optimized except for adipose
#' @param df_ad Dataframe of adipose volume to be optimized
#' @param bw_targ Target body weight in kg
#' @param bw_mean Mean body weight in kg
#' @param ht_rel Relative height ratio
#' @param normOrgan Vector of normally distributed organs
#' @param optVols Vector of non-optimized organ volumes
#' @return A probability to be minimized by an optimizer
#' @importFrom magrittr %>%
#' @importFrom stats dnorm dlnorm
#' @importFrom dplyr mutate
#' @keywords internal
optimVolsFunc <- function(optVols, df_opt, df_ad, bw_targ, bw_mean, ht_rel, normOrgan){
  df_opt$pertMeans <- optVols
  ad <- bw_targ - sum(df_opt$pertMeans)
  p_ad <- dlnorm(ad, meanlog=log(df_ad$means_scaled), sdlog=log(df_ad$std_scaled))
  skFactor <- sqrt(bw_targ/(bw_mean*ht_rel^2))
  df_opt <- df_opt %>% mutate(means_scaled = ifelse(organ %in% c("sk"),
                                                    means_scaled + skFactor,
                                                    means_scaled)) #get a new distribution for skin
  df_opt <- suppressWarnings(df_opt %>% mutate(probs = ifelse(organ %in% normOrgan, dnorm(optVols, mean=means_scaled, sd=std_scaled),
                                                              dlnorm(optVols, meanlog=log(means_scaled), sdlog=log(std_scaled))))) #needs to change
  p <- -sum(c(log(df_opt$probs), log(p_ad)))  #optimize for the negative log probability because we want to maximize not minimize

  return(p)
}

######################################

#' Test if inputs to genPop are compliant
#'
#' Takes in the target covariates and returns an error or a warning if inputs to genPop are not compliant
#'
#' @param minBW Minimum body weight in kg
#' @param maxBW Maximum body weight in kg
#' @param minHT Minimum height in m
#' @param maxHT Maximum height in m
#' @param minBMI Minimum body mass index in kg/m^2
#' @param maxBMI Maximum body mass index in kg/m^2
#' @return An error or a warning if inputs to genPop are not compliant
#' @keywords internal
test_genPopInput <- function(minBW, maxBW, minHT, maxHT, minBMI, maxBMI){
  if((!is.null(minBW) && is.null(maxBW)) | (is.null(minBW) && !is.null(maxBW))) stop("Minimum and maximum values for body weight are required")
  if((!is.null(minHT) && is.null(maxHT)) | (is.null(minHT) && !is.null(maxHT))) stop("Minimum and maximum values for height are required")
  if((!is.null(minBMI) && is.null(maxBMI)) | (is.null(minBMI) && !is.null(maxBMI))) stop("Minimum and maximum values for BMI are required")
  if(!is.null(minBW) && !is.null(minHT) && !is.null(minBMI)) warning("Body weight and height ranges will dictate BMI range; input BMI ranges will be ignored")
}

######################################

#' Sample covariates
#'
#' Takes in the a dataframe of filtered NHANES data and minimum BMI and height and returns a named list of target body weight, height and BMI
#'
#' @param dat Dataframe of NHANES data filtered on a specific age and sex
#' @param minBW Minimum body weight in kg
#' @param maxBW Maximum body weight in kg
#' @param minHT Minimum height in m
#' @param maxHT Maximum height in m
#' @param minBMI Minimum body mass index in kg/m^2
#' @param maxBMI Maximum body mass index in kg/m^2
#' @return A named list of target body weight, height and BMI
#' @importFrom magrittr %>%
#' @importFrom stats sd
#' @importFrom truncnorm rtruncnorm
#' @keywords internal
sampleCov <- function(dat, minBW, maxBW, minHT, maxHT, minBMI, maxBMI){
  if(is.null(minBMI)){
    bw_targ <- exp(rtruncnorm(1, a=log(minBW), b=log(maxBW), mean=mean(log(dat$BW)), sd=sd(log(dat$BW))))
    ht_targ <- exp(rtruncnorm(1, a=log(minHT), b=log(maxHT), mean=mean(log(dat$HT/100)), sd=sd(log(dat$HT/100))))
    bmi_targ <- bw_targ/ht_targ^2
  }else if(is.null(minHT)){
    bw_targ <- exp(rtruncnorm(1, a=log(minBW), b=log(maxBW), mean=mean(log(dat$BW)), sd=sd(log(dat$BW))))
    bmi_targ <- exp(rtruncnorm(1, a=log(minBMI), b=log(maxBMI), mean=mean(log(dat$BMI)), sd=sd(log(dat$BMI))))
    ht_targ <- sqrt(bw_targ/bmi_targ)
  }else{
    ## the final option is same as is.null(minBW) since we can't independently sample the 3 covs but 1 has to be derived from the others and I chose weight
    bmi_targ <- exp(rtruncnorm(1, a=log(minBMI), b=log(maxBMI), mean=mean(log(dat$BMI)), sd=sd(log(dat$BMI))))
    ht_targ <- exp(rtruncnorm(1, a=log(minHT), b=log(maxHT), mean=mean(log(dat$HT/100)), sd=sd(log(dat$HT/100))))
    bw_targ <- bmi_targ*ht_targ^2
  }

  l <- list(bw_targ=bw_targ, ht_targ=ht_targ, bmi_targ=bmi_targ)
  return(l)
}

######################################

#' Sample individual parameters
#'
#' Takes in number of subjects, minimum and maximum age, sex, minimum BMI and height and the a dataframe of filtered NHANES data and minimum BMI and height and returns a named list of target body weight, height and BMI
#'
#' @param nSubj Number of subjects
#' @param minAge Minimum age in years
#' @param maxAge Maximum age in years
#' @param is.male if TRUE, individual is male
#' @param minBW Minimum body weight in kg
#' @param maxBW Maximum body weight in kg
#' @param minHT Minimum height in m
#' @param maxHT Maximum height in m
#' @param minBMI Minimum body mass index in kg/m^2
#' @param maxBMI Maximum body mass index in kg/m^2
#' @param optimize if TRUE, an optimization step is done for each individual
#' @param addBC if TRUE, blood content will be added to each organ
#' @param mab if TRUE, mAb model physiological parameters are returned
#' @return A probability to be minimized by an optimizer
#' @importFrom magrittr %>%
#' @importFrom stats runif
#' @keywords internal
sampleIndPars <- function(nSubj, minAge, maxAge, is.male, minBW, maxBW, minHT, maxHT, minBMI, maxBMI, optimize=FALSE, addBC=TRUE, mab=FALSE){
  pars <- rep(list(), nSubj)

  vol_test <- NA
  ind <- 1
  counter <- 1
  while(ind <= nSubj & counter <= 100){
    age <- round(runif(1, minAge, maxAge))
    sex <- ifelse(is.male, 1, 2)
    dat <- nhanesData %>% filter(AGE_YR == age & SEX == sex)

    covs <- sampleCov(dat, minBW=minBW, maxBW=maxBW, minHT=minHT, maxHT=maxHT, minBMI=minBMI, maxBMI=maxBMI)

    bw_targ <- covs$bw_targ
    ht_targ <- covs$ht_targ
    bmi_targ <- covs$bmi_targ

    if(mab){
      t <- try(suppressWarnings(pars[[ind]] <- genInd_mab(age=age, is.male=is.male, bw_targ=bw_targ, ht_targ=ht_targ, bmi_targ=bmi_targ)), silent = TRUE)  #get the individual parameters
      if(!"try-error" %in% class(t)) vol_test <- pars[[ind]]$V_Liver
    }else{
      t <- try(suppressWarnings(pars[[ind]] <- genInd(age=age, is.male=is.male, bw_targ=bw_targ, ht_targ=ht_targ, bmi_targ=bmi_targ, optimize=optimize, addBC=addBC)), silent = TRUE)  #get the individual parameters
      if(!"try-error" %in% class(t)) vol_test <- pars[[ind]]$Vli
    }

    if("try-error" %in% class(t) | is.na(vol_test)){
      ind <- ind
      counter <- counter + 1
    }else{
      ind <- ind + 1
      counter <- 1
    }
    if(counter == 100) stop("No generated physiological parameters after 100 iterations. Consider modifying input parameters")
  }
  return(pars)
}

######################################

#' Test if inputs to calcBP are compliant
#'
#' Takes in the logP and method for calcBP and returns an error if inputs to calcBP are not compliant
#'
#' @param logP Partition coefficient of a molecule between an aqueous and lipophilic phases, usually octanol and water; measurement of lipophilicity
#' @param method BP calculation method; 1=logP-dependent method, 2=fup-dependent method
#' @return An error if inputs to calcBP are not compliant
#' @keywords internal
test_calcBPInput <- function(logP, fup, method){
  if(method == 1 && fup == 1) stop("fup must be >= 1 for method=1")
  if(method == 2 && is.null(logP)) stop("logP is required for method=2")
}

######################################

#' Get slope and intercept for calcBP or calcFup
#'
#' Takes in the method, generic type of molecule, and the function of interest (calcBP or calcFub) and returns the intercept and slope
#'
#' @param type Type of molecule; can be total, acid, base, or neutral
#' @param method Calculation method; 1=fup-dependent method, 2=logP-dependent method
#' @details Source: https://pubmed.ncbi.nlm.nih.gov/20549836/
#' @return A named list with the intercept and slope values
#' @keywords internal
getInterceptSlope <- function(type="total", method=NULL, func){
  if(func=="calcBP"){
    if(type=="total" & method==1){
      intercept <- 0.208
      slope <- 0.617
    }else if(type=="total" & method==2){
      intercept <- 0.404
      slope <- 0.181
    }else if(type=="acid" & method==1){
      intercept <- 0.131
      slope <- 0.328
    }else if(type=="acid" & method==2){
      intercept <- 0.338
      slope <- 0.049
    }else if(type=="base" & method==1){
      intercept <- 0.232
      slope <- 0.622
    }else if(type=="base" & method==2){
      intercept <- 0.247
      slope <- 0.252
    }else if(type=="neutral" & method==1){
      intercept <- 0.259
      slope <- 0.686
    }else if(type=="neutral" & method==2){
      intercept <- 0.385
      slope <- 0.202
    }
  }else if(func=="calcFup"){
    if(type=="total"){
      intercept <- 0.376
      slope <- 0.252
    }else if(type=="acid"){
      intercept <- 0.53
      slope <- 0.048
    }else if(type=="base"){
      intercept <- 0.003
      slope <- 0.416
    }else if(type=="neutral"){
      intercept <- 0.229
      slope <- 0.26
    }
  }

  l <- list(intercept=intercept, slope=slope)
  return(l)
}

######################################
######################################
