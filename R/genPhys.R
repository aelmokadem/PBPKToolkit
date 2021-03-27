## load databases
nhanesData <- readRDS(system.file("genPhysData", "nhanesData.Rds", package="mrgPBPK"))  #load nhanes anthropometric data
icrpData <- readRDS(system.file("genPhysData", "icrpParams.Rds", package="mrgPBPK"))  #load physiology data from ICRP publication 89
SF <- readRDS(system.file("genPhysData", "allometricSF.Rds", package="mrgPBPK"))  #allometric scaling factors; source: https://github.com/Open-Systems-Pharmacology/OSPSuite.Documentation/wiki/Create-Individual-Algorithm
normSD <- readRDS(system.file("genPhysData", "normalOrgansData.Rds", package="mrgPBPK"))  #normally distributed organs' sd
lnormSD <- readRDS(system.file("genPhysData", "lnormalOrgansData.Rds", package="mrgPBPK"))  #log normally distributed organs' sd
flow <- readRDS(system.file("genPhysData", "organBloodFlow.Rds", package="mrgPBPK"))  #organ blood flows
BC <- readRDS(system.file("genPhysData", "organBloodCont.Rds", package="mrgPBPK"))  #organ blood content

#' Generate individual physiological parameters
#'
#' Takes in desired individual demographics and generates the individual physiological parameters
#'
#' @param age Age of the individual
#' @param is.male if TRUE, individual is male
#' @param bw_targ Target body weight
#' @param ht_targ Target height
#' @param bmi_targ Target body mass index
#' @param optimize if TRUE, an optimization step is done
#' @return Named list with physiological parameters for the desired individual
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom dplyr bind_rows
#' @importFrom dplyr left_join
#' @importFrom nloptr newuoa
#' @export

#This function generates the desired individual parameters; using linear interpolaton
#source is mainly the PK-Sim Github page https://github.com/Open-Systems-Pharmacology/OSPSuite.Documentation/wiki/Create-Individual-Algorithm

genInd <- function(age, is.male, bw_targ=NULL, ht_targ=NULL, bmi_targ=NULL, optimize=F){
  #nhanesData is nhanes anthropometric dataset; icrpData are physiological parameters from ICRP; SF is allomeric scaling factors
  #BC is organ relative blood content

  if(is.null(bmi_targ)) bmi_targ <- bw_targ/ht_targ^2
  if(is.null(bw_targ)) bw_targ <- bmi_targ*ht_targ^2
  if(is.null(ht_targ)) ht_targ <- sqrt(bw_targ/bmi_targ)

  ##filter nhanes database for desired age and gender
  if(is.male){
    sex <- 1
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
    sex <- 2
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

  #get ranges for body weight and height correlated with age and sex
  rangeBW <- quantile(dat$BW, c(0.025, 0.975))  #get the 95% range of body weights
  rangeHT <- quantile(dat$HT/100, c(0.025, 0.975))  #get the 95% range of heights
  rangeBMI <- quantile(dat$BMI, c(0.025, 0.975))  #get the 95% range of body mass index

  #throw error when target measurements are out of range
  if (bw_targ < rangeBW[1] || bw_targ > rangeBW[2] || ht_targ < rangeHT[1] || ht_targ > rangeHT[2] || bmi_targ < rangeBMI[1] || bmi_targ > rangeBMI[2]){
    stop("Input values are out of range")
  }

  ##get mean bw, ht and bmi
  bw_mean <- exp(mean(log(dat$BW)))  #geomteric mean for lognormally distributed weights
  ht_mean <- mean(dat$HT)/100  #arithmetic mean for normally distributed heights in m
  bmi_mean <- exp(mean(log(dat$BMI)))  #geomteric mean for lognormally distributed BMI; Jim and I think this is better

  ##get mean physiological parameters scaled by linear interpolation with height
  df2 <- df %>% select(-c(bw, ht, bsa))  #get rid of anthropometric measurements
  ht_ref <- df$ht/100  #get reference body heights from ICRP in m
  linIntFns <- apply(df2, 2, FUN=function(x) approxfun(ht_ref, x, rule=2))  #linear interpolation functions with height; if outside range use closest value
  physPars <- lapply(linIntFns, FUN=function(x) x(ht_mean))  #get list of mean physiological parameters

  ##########piling mean physiological values and standard deviations############
  df_temp <- data.frame(organ=names(physPars), means=as.numeric(physPars))  #get scaling for blood and cardiac output
  blVol <- df_temp$means[df_temp$organ == "bl"]  #extract blood volume
  df_temp <- bind_rows(df_temp, data.frame(organ=c("ve","ar"), means=c(0, 0)))  #temporarily use whole blood volume for arterial and venous blood vols
  df_temp <- df_temp %>% filter(organ != "bl")  #remove blood volume from df
  l_temp <- split(df_temp, df_temp$organ %in% c("co"))  #separate vols from co
  df_vols <- l_temp[["FALSE"]]
  df_co <- l_temp[["TRUE"]]

  ################################# getting organ volumes ##########################################
  #add blood content to vols
  df_vols <- left_join(df_vols, BC, by="organ")
  df_vols <- df_vols %>% mutate(means=means + (bloodPerc*blVol/100))
  df_temp <- bind_rows(df_vols %>% select(organ, means), df_co)  #join with co again

  normOrgan <- normSD$organ
  lnormOrgan <- lnormSD$organ
  df_norm <- df_temp %>% filter(organ %in% normOrgan)  #prepare a df for normally distributed organs
  df_norm <- left_join(df_norm, normSD, by="organ")  #integrate with stds
  df_norm <- df_norm %>% mutate(std = cv*means/100) %>% select(-cv)  #get standard deviations from CVs then drop it
  df_lnorm <- df_temp %>% filter(organ %in% lnormOrgan)  #prepare a df for lognormally distributed organs
  df_lnorm <- left_join(df_lnorm, lnormSD, by="organ")  #integrate with stds
  names(df_lnorm)[names(df_lnorm) == "geomSD"] <- "std"  #change geomSD to std
  df3 <- bind_rows(df_norm, df_lnorm)  #get combined dataframe for normally and lognormally distributed organs

  ##scaling of physiological parameters
  ht_rel <- ht_targ/ht_mean  #get relative height ratio
  SF2 <- bind_rows(SF, data.frame(organ="co", factor=0.75))  #add allometric scaling factor for cardiac output
  df3 <- left_join(df3, SF2, by=("organ"))  #integrate allometric scaling factors
  df3 <- df3 %>% mutate(means_scaled = ifelse(organ %in% normOrgan, means*ht_rel^factor, means + log(ht_rel^factor)))  #get scaled mean values
  df3 <- df3 %>% mutate(std_scaled = ifelse(organ %in% normOrgan, std*ht_rel^factor, std))  #get scaled standard deviations
  df4 <- df3 %>% filter(organ != "co")  #remove cardiac output

  ######################### optimization of organ volumes #################################
  #set initial values for volumes
  l_temp <- split(df4, df4$organ %in% c("ad"))  #separate vols from ad
  df_opt <- l_temp[["FALSE"]]
  df_ad <- l_temp[["TRUE"]]

  #df_opt <- df4 %>% filter(organ != "ad") #remove adipose as we will not optimize for it
  pert <- runif(1)  #perturbation to parameters
  df_opt <- df_opt %>% mutate(pertMeans = ifelse(organ %in% normOrgan,
                                                        means_scaled + std_scaled*pert,
                                                        means_scaled*exp(log(std_scaled)*pert))) %>%
    mutate(pertMeans = ifelse(organ == "sk", pertMeans*sqrt((bw_targ/(bw_mean*ht_rel^2))), pertMeans))#temporary fix
  ad <- bw_targ - sum(df_opt$pertMeans) #compute adipose volume by subtracting current weight from target bw_targ
  p_ad <- plnorm(ad, meanlog=log(df_ad$means_scaled), sdlog=log(df_ad$std_scaled))  #compute prob for adipose vol

  optVols <- df_opt$pertMeans  #initial parameter values

  #set the optimization function to return probability
  optimVols <- function(optVols){
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

  if(isTRUE(optimize)){
    optMod <- suppressWarnings(newuoa(optVols, optimVols))
    df_opt <- df_opt %>% mutate(optimized=optMod$par)
  }else{
    df_opt <- df_opt %>% mutate(optimized=means_scaled)
  }

  ad <- bw_targ - sum(df_opt$optimized)

  v_ov <- c(df_opt$optimized, ad)

  l_ov <- as.list(v_ov)  #get the final scaled organ volumes in a list
  names(l_ov) <- paste("V", c(as.character(df_opt$organ), "ad"), sep="")  #name the list


  ################################# getting organ blood flows #####################################
  co <- df3$means_scaled[df3$organ == "co"]
  flow2 <- flow %>% mutate(bfs = (flowPerc*co/100)*60) %>% filter(!organ %in% c("ve","ar","lu","li")) #get flow rates in L/h and remove lungs, liver, arterial and venous blood flows as they will be calculated within the model
  flow2 <- flow2 %>%
    mutate(sds = 0.05*bfs,
                  bfs2 = rnorm(nrow(flow2), mean=bfs, sd=sds))
  l_bf <- as.list(flow2$bfs2)  #list of blood flows
  names(l_bf) <- paste("Q", flow2$organ, sep="")  #name the list

  ##getting final parameter list
  pars <- c(l_ov, l_bf, BW=bw_targ, HT=ht_targ, BMI=bmi_targ, SEX=sex)
}


#####################

#' Generate population physiological parameters
#'
#' Takes in desired population demographics and generates the population physiological parameters
#'
#' @param nSubj Number of subjects
#' @param minAge Minimum age
#' @param maxAge Maximum age
#' @param femPerc Percentage of females
#' @param minBW Minimum body weight
#' @param maxBW Maximum body weight
#' @param minHT Minimum height
#' @param maxHT Maximum height
#' @param minBMI Minimum body mass index
#' @param maxBMI Maximum body mass index
#' @param optimize if TRUE, an optimization step is done for each individual
#' @return List of named lists with physiological parameters for each individual in the population
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom dplyr bind_rows
#' @importFrom truncnorm rtruncnorm
#' @export
#This function generates the desired individual parameters. It takes in the model (mod), number of subjects (n), ranges for age, height
#and weight and and the percentage of females in the population
#library(tidyverse)
#library(Runuran)

genPop <- function(nSubj, minAge, maxAge, femPerc, minBW = NULL, maxBW = NULL, minHT = NULL, maxHT = NULL, minBMI = NULL, maxBMI = NULL, optimize=F){

  ##age is in years; weight is in kg; height is in m
  nFemale <- (femPerc/100)*nSubj
  nMale <- nSubj-nFemale
  pars_m <- rep(list(), nMale)
  pars_f <- rep(list(), nFemale)

  ind <- 1
  counter <- 1
  while(ind <= nMale & counter <= 100){

    age <- round(runif(1, minAge, maxAge))
    is.male <- TRUE
    dat <- nhanesData %>% dplyr::filter(AGE_YR == age & SEX == 1)

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

    t <- try(pars_m[[ind]] <- genInd(age=age, is.male=is.male, bw_targ=bw_targ, ht_targ=ht_targ, bmi_targ=bmi_targ, optimize=optimize))  #get the individual parameters

    if("try-error" %in% class(t)){
      ind <- ind
      counter <- counter + 1
    }else{
      ind <- ind + 1
      counter <- 1
    }
  }

  ind <- 1
  counter <- 1
  while(ind <= nFemale & counter <= 100){

    # test for failed subjects to resample
    vol_test <- NA
    age <- round(runif(1, minAge, maxAge))
    is.male <- FALSE
    dat <- nhanesData %>% dplyr::filter(AGE_YR == age & SEX == 2)

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

    t <- try(pars_f[[ind]] <- genInd(age=age, is.male=is.male, bw_targ=bw_targ, ht_targ=ht_targ, bmi_targ=bmi_targ, optimize=optimize))  #get the individual parameters and get error if input out of range and try again

    if("try-error" %in% class(t)){
      ind <- ind
      counter <- counter + 1
    }else{
      ind <- ind + 1
      counter <- 1
    }
  }

  pars <- c(pars_m, pars_f)

  # add IDs
  pars2 <- lapply(1:length(pars), function(i) return(c(pars[[i]],list(ID = i))))

  return(pars2)
}


