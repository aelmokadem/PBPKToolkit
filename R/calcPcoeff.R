#This function generates a list of the desired partition coefficients

## load tissue composition data
TCData <- readRDS(system.file("calcPcoeffData", "unified_tissue_comp.Rds", package="mrgPBPK"))  #organ blood content

## Poulin and Theil

#' Calculate partition coefficients for a molecule based on the Poulin and Theil method
#'
#' Takes in the molecule's physicochemical properties and returns the tissue:plasma partition coefficients based on the Poulin and Theil method
#'
#' @param logP Partition coefficient of a molecule between an aqueous and lipophilic phases, usually octanol and water; measurement of lipophilicity
#' @param pKa Negative log of the acid dissociation constant; measurement of the acidic strength of the molecule
#' @param fup Unbound fraction of the molecule in plasma
#' @param BP Blood:plasma concentration ratio
#' @param type Type of the molecule; 1=neutral, 2=monoprotic acid, 3=monoprotic base, 4=diprotic acid, 5=diprotic base, 6=monoprotic acid monoprotic base (acid comes first), 7=triprotic acid, 8=triprotic base, 9=diprotic acid monoprotic base (first two are acid), 10=diprotic base monoprotic acid (first one is acid)
#' @param dat Dataframe containing tissue composition data; columns are: tissue, f_water=water fraction, f_lipids=lipids fraction, f_proteins=protein fraction, f_pl=phospholipids fraction, f_n_l=neutral lipids fraction, f_n_pl=neutral phospholipids fraction, f_a_pl=acidic phospholipids fraction, pH, f_ew=extracellular fraction, f_iw=intracellular fraction, AR, LR  Prediction method; PT=Poulin and Theil, Berez=Berezhkovskiy, RR=Rodgers and Rowland, Schmitt=Schmitt, pksim=PK-Sim standard
#' @return A named list with tissue:plasma partition coefficients
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @keywords internal
calcKp_PT <- function(logP, pKa, fup, BP=1, type=1, dat){

  dat_all <- dat %>% filter(!tissue %in% c("Plasma","Adipose","RBCs"))

  n <- length(dat$tissue)
  Kp_all <- vector(mode = "numeric", length = n)

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

  pH <- dat$pH[dat$tissue == "Adipose"]
  logD <- 1.115*logP-1.35 #logD is the olive oil:buffer(water) partition coefficient of nonionized species

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

## Rodgers and Rowland

#' Calculate partition coefficients for a molecule based on the Rodgers and Rowland method
#'
#' Takes in the molecule's physicochemical properties and returns the tissue:plasma partition coefficients based on the Rodgers and Rowland method
#'
#' @param logP Partition coefficient of a molecule between an aqueous and lipophilic phases, usually octanol and water; measurement of lipophilicity
#' @param pKa Negative log of the acid dissociation constant; measurement of the acidic strength of the molecule
#' @param fup Unbound fraction of the molecule in plasma
#' @param BP Blood:plasma concentration ratio
#' @param type Type of the molecule; 1=neutral, 2=monoprotic acid, 3=monoprotic base, 4=diprotic acid, 5=diprotic base, 6=monoprotic acid monoprotic base (acid comes first), 7=triprotic acid, 8=triprotic base, 9=diprotic acid monoprotic base (first two are acid), 10=diprotic base monoprotic acid (first one is acid)
#' @param dat Dataframe containing tissue composition data; columns are: tissue, f_water=water fraction, f_lipids=lipids fraction, f_proteins=protein fraction, f_pl=phospholipids fraction, f_n_l=neutral lipids fraction, f_n_pl=neutral phospholipids fraction, f_a_pl=acidic phospholipids fraction, pH, f_ew=extracellular fraction, f_iw=intracellular fraction, AR, LR  Prediction method; PT=Poulin and Theil, Berez=Berezhkovskiy, RR=Rodgers and Rowland, Schmitt=Schmitt, pksim=PK-Sim standard
#' @return A named list with tissue:plasma partition coefficients
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @keywords internal
calcKp_RR <- function(logP, pKa, fup, BP=1, type=1, dat){

  dat_all <- dat %>% filter(!tissue %in% c("RBCs", "Adipose", "Plasma"))  #df for all tissues except for adipose, RBCs, and plasma
  dat_ad <- dat %>% filter(tissue == "Adipose")  #df for adipose
  dat_rbc <- dat %>% filter(tissue == "RBCs") #df for RBCs
  dat_plas <- dat %>% filter(tissue == "Plasma") #df for aplasma

  pH_IW <- 7       #pH of intracellular tissue water
  pH_P <- 7.4      #pH of plasma
  pH_RBC <- 7.22    #pH of blood cells
  P <- 10^(logP)   # octonal:water partition coeff
  logP_OW <- 1.115*logP - 1.35 #oil:water partition coeff
  P_OW <- 10^(logP_OW)
  Ka <- 10^(-pKa)
  HCT <- 0.45 #hematocrit

  #Calculate Kp values
  Kpu_bc <- (HCT - 1 + BP)/(HCT)

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
              10^(pKa[3]-pH_P)+10^(pKa[3]+pka[2]-2*pH_P)+10^(pKa[1]+pKa[2]+pKa[3]-3*pH_P),
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
              10^(pKa[3]-pH_RBC)+10^(pKa[3]+pka[2]-2*pH_RBC)+10^(pKa[1]+pKa[2]+pKa[3]-3*pH_RBC),
              #9-diprotic acid monoprotic base (first two are acid)
              10^(pKa[3]-pH_RBC)+10^(pH_RBC-pKa[1])+10^(2*pH_RBC-pKa[1]-pKa[2]),
              #10-diprotic base monoprotic acid (first one is acid)
              10^(pH_RBC-pKa[1])+10^(pKa[3]-pH_RBC)+10^(pKa[2]+pKa[3]-2*pH_RBC))

  Ka_PR <- (1/fup - 1 - (P*dat_plas$f_n_l + (0.3*P + 0.7)*dat_plas$f_n_pl)/(1+Y))
  Ka_AP <- (Kpu_bc - (1 + Z)/(1 + Y)*dat_rbc$f_iw - (P*dat_rbc$f_n_l + (0.3*P + 0.7)*dat_rbc$f_n_pl)/(1 + Y)) * (1 + Y)/dat_rbc$f_a_pl/Z

  # Assign the moderate to strong bases type_calc=1 and everything else type_calc=2
  type_calc <- ifelse((type==3 & pKa[1]>7) | (type==5 & pKa[1] >7) | (type==6 & pKa[2] > 7) | (type==8 & pKa[1] > 7) | (type==9 & pKa[3]>7) | (type==10 & pKa[2]>7), 1,2)

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

## Berezhkovskiy

#' Calculate partition coefficients for a molecule based on the Berezhkovskiy method
#'
#' Takes in the molecule's physicochemical properties and returns the tissue:plasma partition coefficients based on the Berezhkovskiy method
#'
#' @param logP Partition coefficient of a molecule between an aqueous and lipophilic phases, usually octanol and water; measurement of lipophilicity
#' @param pKa Negative log of the acid dissociation constant; measurement of the acidic strength of the molecule
#' @param fup Unbound fraction of the molecule in plasma
#' @param BP Blood:plasma concentration ratio
#' @param type Type of the molecule; 1=neutral, 2=monoprotic acid, 3=monoprotic base, 4=diprotic acid, 5=diprotic base, 6=monoprotic acid monoprotic base (acid comes first), 7=triprotic acid, 8=triprotic base, 9=diprotic acid monoprotic base (first two are acid), 10=diprotic base monoprotic acid (first one is acid)
#' @param dat Dataframe containing tissue composition data; columns are: tissue, f_water=water fraction, f_lipids=lipids fraction, f_proteins=protein fraction, f_pl=phospholipids fraction, f_n_l=neutral lipids fraction, f_n_pl=neutral phospholipids fraction, f_a_pl=acidic phospholipids fraction, pH, f_ew=extracellular fraction, f_iw=intracellular fraction, AR, LR  Prediction method; PT=Poulin and Theil, Berez=Berezhkovskiy, RR=Rodgers and Rowland, Schmitt=Schmitt, pksim=PK-Sim standard
#' @return A named list with tissue:plasma partition coefficients
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @keywords internal
calcKp_Berez <- function(logP, pKa, fup, BP=1, type=1, dat){

  dat_all <- dat %>% filter(!tissue %in% c("Plasma","Adipose","RBCs"))

  n <- length(dat$tissue)
  Kp_all <- vector(mode = "numeric", length = n)

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

  fut <- 1/(1+((1-fup)/fup)*0.5)

  pH <- dat$pH[dat$tissue == "Adipose"]
  logD <- 1.115*logP-1.35 #logD is the olive oil:buffer(water) partition coefficient of nonionized species

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

## Schmitt

#' Calculate partition coefficients for a molecule based on the Schmitt method
#'
#' Takes in the molecule's physicochemical properties and returns the tissue:plasma partition coefficients based on the Schmitt method
#'
#' @param logP Partition coefficient of a molecule between an aqueous and lipophilic phases, usually octanol and water; measurement of lipophilicity
#' @param pKa Negative log of the acid dissociation constant; measurement of the acidic strength of the molecule
#' @param fup Unbound fraction of the molecule in plasma
#' @param type Type of the molecule; 1=neutral, 2=monoprotic acid, 3=monoprotic base, 4=diprotic acid, 5=diprotic base, 6=monoprotic acid monoprotic base (acid comes first), 7=triprotic acid, 8=triprotic base, 9=diprotic acid monoprotic base (first two are acid), 10=diprotic base monoprotic acid (first one is acid)
#' @param dat Dataframe containing tissue composition data; columns are: tissue, f_water=water fraction, f_lipids=lipids fraction, f_proteins=protein fraction, f_pl=phospholipids fraction, f_n_l=neutral lipids fraction, f_n_pl=neutral phospholipids fraction, f_a_pl=acidic phospholipids fraction, pH, f_ew=extracellular fraction, f_iw=intracellular fraction, AR, LR  Prediction method; PT=Poulin and Theil, Berez=Berezhkovskiy, RR=Rodgers and Rowland, Schmitt=Schmitt, pksim=PK-Sim standard
#' @return A named list with tissue:plasma partition coefficients
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @keywords internal
calcKp_Schmitt <- function(logP, pKa, fup, type = 1, dat){
  #logMA is the log of membrane affinity = phosphatidylcholine:water (neutral phospholipid:water) partition coefficient;
  #we can use the available measurement of lipophilicity instead (logP or logD); from Schmitt, Walter (2008)

  dat_all <- dat %>% filter(!tissue %in% c("RBCs", "Plasma"))  #df for all tissues except for adipose and RBCs

  logMA <- logP  #in case we don't have a direct logMA
  K_n_pl <- 10^logMA    #neutral phospholipids:water partition coefficient
  K_protein <- ((0.81 + 0.11 * K_n_pl)/24.92)*5
  pH <- dat_all$pH
  alpha <- 1e-3  #ratio between ditribution coefficient at given pH (D) and that in neutral form (D0)

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

  kp <- (dat_all$f_water+(K_n_l*dat_all$f_n_l)+(K_n_pl*dat_all$f_n_pl)+(K_a_pl*dat_all$f_a_pl)+(K_protein*dat_all$f_proteins))*fup

  dat2 <- data.frame(tissue=dat_all$tissue, Kp=kp)
  name <- dat2$tissue %>% substr(1,2) %>% tolower()
  name <- paste("Kp", name, sep="")
  Kp <- as.list(dat2$Kp)
  names(Kp) <- name

  return(Kp)
}

######################

## PK-Sim Standard

#' Calculate partition coefficients for a molecule based on the PK_Sim method
#'
#' Takes in the molecule's physicochemical properties and returns the tissue:plasma partition coefficients based on the PK-Sim method
#'
#' @param logP Partition coefficient of a molecule between an aqueous and lipophilic phases, usually octanol and water; measurement of lipophilicity
#' @param fup Unbound fraction of the molecule in plasma
#' @param dat Dataframe containing tissue composition data; columns are: tissue, f_water=water fraction, f_lipids=lipids fraction, f_proteins=protein fraction, f_pl=phospholipids fraction, f_n_l=neutral lipids fraction, f_n_pl=neutral phospholipids fraction, f_a_pl=acidic phospholipids fraction, pH, f_ew=extracellular fraction, f_iw=intracellular fraction, AR, LR  Prediction method; PT=Poulin and Theil, Berez=Berezhkovskiy, RR=Rodgers and Rowland, Schmitt=Schmitt, pksim=PK-Sim standard
#' @return A named list with tissue:plasma partition coefficients
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @keywords internal
calcKp_pksim <- function(logP, fup, dat){
  #logMA is the log of membrane affinity = phosphatidylcholin:water (neutral phospholipid:water) partition coefficient;
  #we can use the available measurement of lipophilicity instead (logP or logD); from Schmitt, Walter (2008)

  dat_all <- dat %>% filter(!tissue %in% c("RBCs", "Plasma"))  #df for all tissues except for adipose and RBCs

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
#' @return A named list with tissue:plasma partition coefficients
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @export
## general function
calcKp <- function(logP, pKa=NULL, fup, BP=1, type=1, method="PT"){
  if(method == "PT"){
    pcoeff <- calcKp_PT(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type, dat=TCData)
  }else if(method == "Berez"){  #Berezhkovskiy
    pcoeff <- calcKp_Berez(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type, dat=TCData)
  }else if(method == "pksim"){  #standard PK-Sim, Willmann et al. 2008
    pcoeff <- calcKp_pksim(logP=logP, fup=fup, dat=TCData)
  }else if(method == "Schmitt"){  #Schmitt, Walter 2008
    pcoeff <- calcKp_Schmitt(logP=logP, pKa=pKa, fup=fup, type=type, dat=TCData)
  }else{ #Rodgers and Rowland
    pcoeff <- calcKp_RR(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type, dat=TCData)
  }

  return(pcoeff)
}
