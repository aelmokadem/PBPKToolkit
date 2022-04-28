#' Calculate partition coefficients for a molecule based on the Poulin and Theil method
#'
#' Takes in the molecule's physicochemical properties and returns the tissue:plasma partition coefficients based on the Poulin and Theil method
#'
#' @param logP Partition coefficient of a molecule between an aqueous and lipophilic phases, usually octanol and water; measurement of lipophilicity
#' @param pKa Negative log of the acid dissociation constant; measurement of the acidic strength of the molecule
#' @param fup Unbound fraction of the molecule in plasma
#' @param BP Blood:plasma concentration ratio
#' @param type Type of the molecule; 1=neutral, 2=monoprotic acid, 3=monoprotic base, 4=diprotic acid, 5=diprotic base, 6=monoprotic acid monoprotic base (acid comes first), 7=triprotic acid, 8=triprotic base, 9=diprotic acid monoprotic base (first two are acid), 10=diprotic base monoprotic acid (first one is acid)
#' @param dat Dataframe containing tissue composition data; columns are: tissue, f_water=water fraction, f_lipids=lipids fraction, f_proteins=protein fraction, f_pl=phospholipids fraction, f_n_l=neutral lipids fraction, f_n_pl=neutral phospholipids fraction, f_a_pl=acidic phospholipids fraction, pH, f_ew=extracellular fraction, f_iw=intracellular fraction, AR=albumin ratio, LR=lipoprotein ratio,  Prediction method; PT=Poulin and Theil, Berez=Berezhkovskiy, RR=Rodgers and Rowland, Schmitt=Schmitt, pksim=PK-Sim standard
#' @return A named list with tissue:plasma partition coefficients
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @details PT: https://pubmed.ncbi.nlm.nih.gov/11782904/
#' @keywords internal
calcKp_PT <- function(logP, pKa, fup, BP=1, type=1, dat){

  dat_all <- dat %>% filter(!tissue %in% c("Plasma","Adipose","RBCs"))
  n <- length(dat$tissue)
  pH <- dat$pH[dat$tissue == "Adipose"]
  Kp_all <- vector(mode = "numeric", length = n)

  TC <- getTissueData(dat)
  Vwp <- TC$Vwp
  Vnlp <- TC$Vnlp
  Vphp <- TC$Vphp
  Vwt <- TC$Vwt
  Vwad <- TC$Vwad
  Vnlt <- TC$Vnlt
  Vnlad <- TC$Vnlad
  Vpht <- TC$Vpht
  Vphad <- TC$Vphad

  logD <- 1.115*logP-1.35  #logD is the olive oil:buffer(water) partition coefficient of nonionized species

  logD_star <- getLogD_star(type, logD, pKa, pH)

  D_star <- 10^logD_star
  Kpad <- ((D_star*(Vnlad+0.3*Vphad)+(1*(Vwad+0.7*Vphad)))/(D_star*(Vnlp+0.3*Vphp)+(1*(Vwp+0.7*Vphp)))) * fup

  P <- 10^logP
  fut <- 1/(1+((1-fup)/fup)*0.5)
  Kpt <- ((P*(Vnlt+0.3*Vpht)+(1*(Vwt+0.7*Vpht)))/(P*(Vnlp+0.3*Vphp)+(1*(Vwp+0.7*Vphp)))) * (fup/fut)

  nms_all <- dat_all$tissue %>% substr(1,2) %>% tolower()
  nms_all <- paste("Kp", nms_all, sep="")
  nms <- c("Kpad",nms_all)
  Kp <- as.list(c(Kpad,Kpt))
  names(Kp) <- nms

  return(Kp)

}

######################

#' Calculate partition coefficients for a molecule based on the Rodgers and Rowland method
#'
#' Takes in the molecule's physicochemical properties and returns the tissue:plasma partition coefficients based on the Rodgers and Rowland method
#'
#' @param logP Partition coefficient of a molecule between an aqueous and lipophilic phases, usually octanol and water; measurement of lipophilicity
#' @param pKa Negative log of the acid dissociation constant; measurement of the acidic strength of the molecule
#' @param fup Unbound fraction of the molecule in plasma
#' @param BP Blood:plasma concentration ratio
#' @param type Type of the molecule; 1=neutral, 2=monoprotic acid, 3=monoprotic base, 4=diprotic acid, 5=diprotic base, 6=monoprotic acid monoprotic base (acid comes first), 7=triprotic acid, 8=triprotic base, 9=diprotic acid monoprotic base (first two are acid), 10=diprotic base monoprotic acid (first one is acid)
#' @param dat Dataframe containing tissue composition data; columns are: tissue, f_water=water fraction, f_lipids=lipids fraction, f_proteins=protein fraction, f_pl=phospholipids fraction, f_n_l=neutral lipids fraction, f_n_pl=neutral phospholipids fraction, f_a_pl=acidic phospholipids fraction, pH, f_ew=extracellular fraction, f_iw=intracellular fraction, AR=albumin ratio, LR=lipoprotein ratio,  Prediction method; PT=Poulin and Theil, Berez=Berezhkovskiy, RR=Rodgers and Rowland, Schmitt=Schmitt, pksim=PK-Sim standard
#' @param Ht Hematocrit value; default is 0.45
#' @return A named list with tissue:plasma partition coefficients
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @details RR: https://pubmed.ncbi.nlm.nih.gov/15858854/ and https://pubmed.ncbi.nlm.nih.gov/16639716/
#' @keywords internal
calcKp_RR <- function(logP, pKa, fup, BP=1, type=1, dat, Ht=0.45){

  dat_all <- dat %>% filter(!tissue %in% c("RBCs", "Adipose", "Plasma"))  #df for all tissues except for adipose, RBCs, and plasma
  dat_ad <- dat %>% filter(tissue == "Adipose")  #df for adipose
  dat_rbc <- dat %>% filter(tissue == "RBCs") #df for RBCs
  dat_plas <- dat %>% filter(tissue == "Plasma") #df for plasma

  pH_IW <- 7       #pH of intracellular tissue water
  pH_P <- 7.4      #pH of plasma
  pH_RBC <- 7.22    #pH of blood cells
  P <- 10^(logP)   # octonal:water partition coeff
  logP_OW <- 1.115*logP - 1.35 #oil:water partition coeff
  P_OW <- 10^(logP_OW)
  Ka <- 10^(-pKa)

  #Calculate Kp values
  Kpu_bc <- (Ht - 1 + BP)/(Ht)

  XYZ <- getXYZ(type, pKa, pH_IW, pH_P, pH_RBC)
  X <- XYZ$X
  Y <- XYZ$Y
  Z <- XYZ$Z

  Ka_PR <- (1/fup - 1 - (P*dat_plas$f_n_l + (0.3*P + 0.7)*dat_plas$f_n_pl)/(1+Y))
  Ka_AP <- (Kpu_bc - (1 + Z)/(1 + Y)*dat_rbc$f_iw - (P*dat_rbc$f_n_l + (0.3*P + 0.7)*dat_rbc$f_n_pl)/(1 + Y)) * (1 + Y)/dat_rbc$f_a_pl/Z

  # Assign the moderate to strong bases type_calc=1 and everything else type_calc=2
  type_calc <- ifelse((type==3 & pKa[1]>7) | (type==5 & pKa[1] > 7) | (type==6 & pKa[2] > 7) | (type==8 & pKa[1] > 7) | (type==9 & pKa[3] > 7) | (type==10 & pKa[2] > 7), 1,2)

  # Re-assign the neutrals type_calc=3
  if(type==1){type_calc=3}  #neutrals

  # Multiply by fup to get Kp rather than Kpu
  if(type_calc==1){  #moderate to strong bases
    Kp_all <- (dat_all$f_ew + ((1 + X)/(1 + Y))*dat_all$f_iw + ((P*dat_all$f_n_l + (0.3*P + 0.7)*dat_all$f_n_pl))/(1 + Y) + (Ka_AP*dat_all$f_a_pl*X)/(1 + Y))*fup  #non lipid
    Kp_ad <- (dat_ad$f_ew + ((1 + X)/(1 + Y))*dat_ad$f_iw + ((P_OW*dat_ad$f_n_l + (0.3*P_OW + 0.7)*dat_ad$f_n_pl))/(1 + Y) + (Ka_AP*dat_ad$f_a_pl*X)/(1 + Y))*fup  #lipid
  }else if(type_calc==2){   #acidic and zwitterions
    Kp_all <- (dat_all$f_ew + ((1 + X)/(1 + Y))*dat_all$f_iw + ((P*dat_all$f_n_l + (0.3*P + 0.7)*dat_all$f_n_pl))/(1 + Y) + (Ka_PR*dat_all$AR*X)/(1 + Y))*fup  #non lipid
    Kp_ad <- (dat_ad$f_ew + ((1 + X)/(1 + Y))*dat_ad$f_iw + ((P_OW*dat_ad$f_n_l + (0.3*P_OW + 0.7)*dat_ad$f_n_pl))/(1 + Y) + (Ka_PR*dat_ad$AR*X)/(1 + Y))*fup #lipid
  }else{  #neutrals
    Kp_all <- (dat_all$f_ew + ((1 + X)/(1 + Y))*dat_all$f_iw + ((P*dat_all$f_n_l + (0.3*P + 0.7)*dat_all$f_n_pl))/(1 + Y) + (Ka_PR*dat_all$LR*X)/(1 + Y))*fup  #non lipid
    Kp_ad <- (dat_ad$f_ew + ((1 + X)/(1 + Y))*dat_ad$f_iw + ((P_OW*dat_ad$f_n_l + (0.3*P_OW + 0.7)*dat_ad$f_n_pl))/(1 + Y) + (Ka_PR*dat_ad$LR*X)/(1 + Y))*fup  #lipid
  }

  nms_all <- dat_all$tissue %>% substr(1,2) %>% tolower()
  nms_all <- paste("Kp", nms_all, sep="")
  nms <- c("Kpad",nms_all)
  Kp <- as.list(c(Kp_ad,Kp_all))
  names(Kp) <- nms

  return(Kp)

}

######################

#' Calculate partition coefficients for a molecule based on the Berezhkovskiy method
#'
#' Takes in the molecule's physicochemical properties and returns the tissue:plasma partition coefficients based on the Berezhkovskiy method
#'
#' @param logP Partition coefficient of a molecule between an aqueous and lipophilic phases, usually octanol and water; measurement of lipophilicity
#' @param pKa Negative log of the acid dissociation constant; measurement of the acidic strength of the molecule
#' @param fup Unbound fraction of the molecule in plasma
#' @param BP Blood:plasma concentration ratio
#' @param type Type of the molecule; 1=neutral, 2=monoprotic acid, 3=monoprotic base, 4=diprotic acid, 5=diprotic base, 6=monoprotic acid monoprotic base (acid comes first), 7=triprotic acid, 8=triprotic base, 9=diprotic acid monoprotic base (first two are acid), 10=diprotic base monoprotic acid (first one is acid)
#' @param dat Dataframe containing tissue composition data; columns are: tissue, f_water=water fraction, f_lipids=lipids fraction, f_proteins=protein fraction, f_pl=phospholipids fraction, f_n_l=neutral lipids fraction, f_n_pl=neutral phospholipids fraction, f_a_pl=acidic phospholipids fraction, pH, f_ew=extracellular fraction, f_iw=intracellular fraction, AR=albumin ratio, LR=lipoprotein ratio,  Prediction method; PT=Poulin and Theil, Berez=Berezhkovskiy, RR=Rodgers and Rowland, Schmitt=Schmitt, pksim=PK-Sim standard
#' @return A named list with tissue:plasma partition coefficients
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @details Berez: https://pubmed.ncbi.nlm.nih.gov/15124219/
#' @keywords internal
calcKp_Berez <- function(logP, pKa, fup, BP=1, type=1, dat){

  dat_all <- dat %>% filter(!tissue %in% c("Plasma","Adipose","RBCs"))
  n <- length(dat$tissue)
  Kp_all <- vector(mode = "numeric", length = n)

  TC <- getTissueData(dat)
  Vwp <- TC$Vwp
  Vnlp <- TC$Vnlp
  Vphp <- TC$Vphp
  Vwt <- TC$Vwt
  Vwad <- TC$Vwad
  Vnlt <- TC$Vnlt
  Vnlad <- TC$Vnlad
  Vpht <- TC$Vpht
  Vphad <- TC$Vphad

  fut <- 1/(1+((1-fup)/fup)*0.5)

  pH <- dat$pH[dat$tissue == "Adipose"]
  logD <- 1.115*logP-1.35 #logD is the olive oil:buffer(water) partition coefficient of nonionized species

  logD_star <- getLogD_star(type, logD, pKa, pH)

  D_star <- 10^logD_star
  Kpad <- ((D_star*(Vnlad+0.3*Vphad)+((Vwad/fut)+0.7*Vphad))/(D_star*(Vnlp+0.3*Vphp)+((Vwp/fup)+0.7*Vphp)))

  P <- 10^logP
  Kpt <- ((P*(Vnlt+0.3*Vpht)+((Vwt/fut)+0.7*Vpht))/(P*(Vnlp+0.3*Vphp)+((Vwp/fup)+0.7*Vphp)))

  nms_all <- dat_all$tissue %>% substr(1,2) %>% tolower()
  nms_all <- paste("Kp", nms_all, sep="")
  nms <- c("Kpad",nms_all)
  # return(nms)
  Kp <- as.list(c(Kpad,Kpt))
  names(Kp) <- nms

  return(Kp)
}

######################

#' Calculate partition coefficients for a molecule based on the Schmitt method
#'
#' Takes in the molecule's physicochemical properties and returns the tissue:plasma partition coefficients based on the Schmitt method
#'
#' @param logP Partition coefficient of a molecule between an aqueous and lipophilic phases, usually octanol and water; measurement of lipophilicity
#' @param pKa Negative log of the acid dissociation constant; measurement of the acidic strength of the molecule
#' @param fup Unbound fraction of the molecule in plasma
#' @param type Type of the molecule; 1=neutral, 2=monoprotic acid, 3=monoprotic base, 4=diprotic acid, 5=diprotic base, 6=monoprotic acid monoprotic base (acid comes first), 7=triprotic acid, 8=triprotic base, 9=diprotic acid monoprotic base (first two are acid), 10=diprotic base monoprotic acid (first one is acid)
#' @param dat Dataframe containing tissue composition data; columns are: tissue, f_water=water fraction, f_lipids=lipids fraction, f_proteins=protein fraction, f_pl=phospholipids fraction, f_n_l=neutral lipids fraction, f_n_pl=neutral phospholipids fraction, f_a_pl=acidic phospholipids fraction, pH, f_ew=extracellular fraction, f_iw=intracellular fraction, AR=albumin ratio, LR=lipoprotein ratio,  Prediction method; PT=Poulin and Theil, Berez=Berezhkovskiy, RR=Rodgers and Rowland, Schmitt=Schmitt, pksim=PK-Sim standard
#' @return A named list with tissue:plasma partition coefficients
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @details Schmitt: https://pubmed.ncbi.nlm.nih.gov/17981004/
#' @keywords internal
calcKp_Schmitt <- function(logP, pKa, fup, type = 1, dat){
  #logMA is the log of membrane affinity = phosphatidylcholine:water (neutral phospholipid:water) partition coefficient;
  #we can use the available measurement of lipophilicity instead (logP or logD); from Schmitt, Walter (2008)

  dat_all <- dat %>% filter(!tissue %in% c("RBCs", "Plasma"))  #df for all tissues except for adipose and RBCs

  logMA <- logP  #in case we don't have a direct logMA
  K_n_pl <- 10^logMA    #neutral phospholipids:water partition coefficient
  K_protein <- ((0.81 + 0.11 * K_n_pl)/24.92)*5
  pH <- dat_all$pH
  alpha <- 1e-3  #ratio between distribution coefficient at given pH (D) and that in neutral form (D0)

  W <- getW(type, pKa, pH)

  Ks <- getKs(type, pKa, pH, alpha, W, K_n_pl)
  K_n_l <- Ks$K_n_l
  K_a_pl <- Ks$K_a_pl

  kp <- (dat_all$f_water+(K_n_l*dat_all$f_n_l)+(K_n_pl*dat_all$f_n_pl)+(K_a_pl*dat_all$f_a_pl)+(K_protein*dat_all$f_proteins))*fup

  dat2 <- data.frame(tissue=dat_all$tissue, Kp=kp)
  name <- dat2$tissue %>% substr(1,2) %>% tolower()
  name <- paste("Kp", name, sep="")
  Kp <- as.list(dat2$Kp)
  names(Kp) <- name

  return(Kp)

}

######################

#' Calculate partition coefficients for a molecule based on the PK_Sim method
#'
#' Takes in the molecule's physicochemical properties and returns the tissue:plasma partition coefficients based on the PK-Sim method
#'
#' @param logP Partition coefficient of a molecule between an aqueous and lipophilic phases, usually octanol and water; measurement of lipophilicity
#' @param fup Unbound fraction of the molecule in plasma
#' @param dat Dataframe containing tissue composition data; columns are: tissue, f_water=water fraction, f_lipids=lipids fraction, f_proteins=protein fraction, f_pl=phospholipids fraction, f_n_l=neutral lipids fraction, f_n_pl=neutral phospholipids fraction, f_a_pl=acidic phospholipids fraction, pH, f_ew=extracellular fraction, f_iw=intracellular fraction, AR=albumin ratio, LR=lipoprotein ratio,  Prediction method; PT=Poulin and Theil, Berez=Berezhkovskiy, RR=Rodgers and Rowland, Schmitt=Schmitt, pksim=PK-Sim standard
#' @return A named list with tissue:plasma partition coefficients
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @details pksim: https://www.tandfonline.com/doi/abs/10.1517/17425255.1.1.159
#' @keywords internal
calcKp_pksim <- function(logP, fup, dat){
  #logMA is the log of membrane affinity = phosphatidylcholin:water (neutral phospholipid:water) partition coefficient;
  #we can use the available measurement of lipophilicity instead (logP or logD); from Schmitt, Walter (2008)

  dat_all <- dat %>% filter(!tissue %in% c("RBCs", "Plasma"))  #df for all tissues except for plasma and RBCs

  logMA <- logP  #in case we don't have a direct logMA
  K_n_pl <- 10^logMA    #neutral phospholipids:water partition coefficient
  K_protein <- ((0.81 + 0.11 * K_n_pl)/24.92)*5 # From PK-Sim (very similar value to the other method)

  kp <- (dat_all$f_water + (K_n_pl*dat_all$f_lipids) + (K_protein*dat_all$f_proteins))*fup

  dat2 <- data.frame(tissue=dat_all$tissue, Kp=kp)
  name <- dat2$tissue %>% substr(1,2) %>% tolower()
  name <- paste("Kp", name, sep="")
  Kp <- as.list(dat2$Kp)
  names(Kp) <- name

  return(Kp)
}

######################

#' Calculate partition coefficients for a molecule
#'
#' Takes in the molecule's physicochemical properties and returns the tissue:plasma partition coefficients
#'
#' @param logP Partition coefficient of a molecule between an aqueous and lipophilic phases, usually octanol and water; measurement of lipophilicity
#' @param pKa Negative log of the acid dissociation constant; measurement of the acidic strength of the molecule
#' @param fup Unbound fraction of the molecule in plasma
#' @param BP Blood:plasma concentration ratio
#' @param type Type of the molecule; 1=neutral, 2=monoprotic acid, 3=monoprotic base, 4=diprotic acid, 5=diprotic base, 6=monoprotic acid monoprotic base (acid comes first), 7=triprotic acid, 8=triprotic base, 9=diprotic acid monoprotic base (first two are acid), 10=diprotic base monoprotic acid (first one is acid)
#' @param method Prediction method; PT=Poulin and Theil, Berez=Berezhkovskiy, RR=Rodgers and Rowland, Schmitt=Schmitt, pksim=PK-Sim standard
#' @param Vss Steady state volume in L
#' @param Vt Named list of tissue volumes generated by genInd
#' @param Kpot Partition coefficient for the "other" compartment to be used when Vss value is provided. If Kpot is left NULL while Vss value is provided, Kpot will be calculated as the average of the non adipose partition coefficients
#' @param Ht Hematocrit value; default is 0.45
#' @return A named list with tissue:plasma partition coefficients
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @details Sources:
#' @details PT: https://pubmed.ncbi.nlm.nih.gov/11782904/
#' @details RR: https://pubmed.ncbi.nlm.nih.gov/15858854/ and https://pubmed.ncbi.nlm.nih.gov/16639716/
#' @details Berez: https://pubmed.ncbi.nlm.nih.gov/15124219/
#' @details Schmitt: https://pubmed.ncbi.nlm.nih.gov/17981004/
#' @details pksim: https://www.tandfonline.com/doi/abs/10.1517/17425255.1.1.159
#' @export
## general function
calcKp <- function(logP, pKa=NULL, fup, BP=1, type=1, method="PT", Vss=NULL, Vt=NULL, Kpot=NULL, Ht=0.45){
  #useful error messages
  if(fup > 1) stop("fup should be <= 1")
  test_pKaTypeMatch(type, pKa)
  test_largeHt(Ht)

  if(method == "PT"){
    pcoeff <- calcKp_PT(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type, dat=TCData)
  }else if(method == "Berez"){  #Berezhkovskiy
    pcoeff <- calcKp_Berez(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type, dat=TCData)
  }else if(method == "pksim"){  #standard PK-Sim, Willmann et al. 2008
    pcoeff <- calcKp_pksim(logP=logP, fup=fup, dat=TCData)
  }else if(method == "Schmitt"){  #Schmitt, Walter 2008
    pcoeff <- calcKp_Schmitt(logP=logP, pKa=pKa, fup=fup, type=type, dat=TCData)
  }else{ #Rodgers and Rowland
    pcoeff <- calcKp_RR(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type, dat=TCData, Ht=Ht)
  }

  test_negativeKps(pcoeff)

  if(!is.null(Vss)){
    test_missingVt(Vt)
    pcoeff <- scaleKp(Kp=pcoeff, Vss=Vss, BP=BP, Vt=Vt, Kpot=Kpot, Ht=Ht)
  }

  return(pcoeff)
}

##########################

#' Calculate blood to plasma concentration ratio
#'
#' Takes in logP, fup, and method and returns the blood to plasma concentration ratio
#'
#' @param logP Partition coefficient of a molecule between an aqueous and lipophilic phases, usually octanol and water; measurement of lipophilicity
#' @param fup Unbound fraction of the molecule in plasma
#' @param type Type of molecule; can be total, acid, base, or neutral
#' @param method BP calculation method; 1=fup-dependent method, 2=logP-dependent method
#' @param Ht Hematocrit value; default is 0.45
#' @return Blood to plasma concentration ratio
#' @details Source: https://pubmed.ncbi.nlm.nih.gov/20549836/
#' @export
calcBP <- function(logP=NULL, fup, type="total", method=1, Ht=0.45){
  test_calcBPInput(logP, fup, method)

  interceptSlope <- getInterceptSlope(type=type, method=method, func="calcBP")
  intercept <- interceptSlope$intercept
  slope <- interceptSlope$slope

  if(method==1){
    logKb <- intercept + slope*log((1-fup)/fup)
  }else{
    logKb <- intercept + slope*logP
  }

  Kb <- exp(logKb)

  BP <- (Kb*fup - 1)*Ht + 1
  return(BP)
}

##########################

#' Calculate unbound fraction in plasma
#'
#' Takes in logP and type of molecule and returns the unbound fraction in plasma
#'
#' @param logP Partition coefficient of a molecule between an aqueous and lipophilic phases, usually octanol and water; measurement of lipophilicity
#' @param type Type of molecule; can be total, acid, base, or neutral
#' @return Unbound fraction in plasma
#' @details Source: https://pubmed.ncbi.nlm.nih.gov/20549836/
#' @export
calcFup <- function(logP, type="total"){
  interceptSlope <- getInterceptSlope(type=type, func="calcFup")
  intercept <- interceptSlope$intercept
  slope <- interceptSlope$slope

  logKfup <- intercept + slope*logP
  Kfup <- exp(logKfup)

  fup <- 1/(Kfup + 1)

  return(fup)
}

