#' Generate individual physiological parameters
#'
#' Takes in desired individual demographics and generates the individual physiological parameters
#'
#' @param age Age of the individual in years
#' @param is.male if TRUE, individual is male
#' @param bw_targ Target body weight in kg
#' @param ht_targ Target height in m
#' @param bmi_targ Target body mass index in kg/m^2
#' @param bsa_targ Target body surface area in m^2; will be used by "Huisinga" method; if missing the value will be calculated from body weight and height
#' @param lbw_targ Target lean body weight in kg; will be used by "Huisinga" method; if missing the value will be approximated by fat free mass
#' @param optimize if TRUE, an optimization step is done; for "Willmann" method
#' @param addBC if TRUE, blood content will be added to each organ
#' @param method Method to use for generating physiological parameters; choices are "Willmann" and "Huisinga"
#' @return Named list with physiological parameters for the desired individual
#' @details Sources: Willmann: https://pubmed.ncbi.nlm.nih.gov/17431751/, Huisinga: https://ascpt.onlinelibrary.wiley.com/doi/10.1038/psp.2012.3
#' @importFrom magrittr %>%
#' @importFrom dplyr filter select mutate bind_rows left_join case_when
#' @importFrom stats approxfun plnorm quantile rnorm runif sd
#' @importFrom nloptr newuoa
#' @export

#This function generates the desired individual parameters; using linear interpolaton
#source is mainly the PK-Sim Github page https://github.com/Open-Systems-Pharmacology/OSPSuite.Documentation/wiki/Create-Individual-Algorithm

genInd <- function(age, is.male, bw_targ=NULL, ht_targ=NULL, bmi_targ=NULL, bsa_targ=NULL, lbw_targ=NULL, optimize=FALSE, addBC=TRUE, method="Willmann"){
  ##filter nhanes database for desired age and sex
  sex <- ifelse(is.male, 1, 2)

  test_genIndInput(age, bw_targ, ht_targ, bmi_targ)

  # inputs for lbw ; https://ascpt.onlinelibrary.wiley.com/doi/10.1038/psp.2012.3
  M_bmi <- 9270.0/216.0
  K_bmi <-  6680.0/216.0
  if(sex == 2){
    M_bmi <- 9270.0/244.0
    K_bmi <-  6680.0/244.0
  }

  if(!is.null(bw_targ) & !is.null(ht_targ)) bmi_targ <- bw_targ/ht_targ^2
  if(is.null(bw_targ)) bw_targ <- bmi_targ*ht_targ^2
  if(is.null(ht_targ)) ht_targ <- sqrt(bw_targ/bmi_targ)
  if(is.null(bsa_targ)) bsa_targ <- sqrt(ht_targ*bw_targ/36)
  if(is.null(lbw_targ)) lbw_targ = (M_bmi/(K_bmi + bmi_targ))*bw_targ  # approximate lbw by fat free mass

  refDatasets <- filterDatasets(age, sex)

  dat <- refDatasets$dat
  df <- refDatasets$df
  normSD <- refDatasets$normSD
  flow <- refDatasets$flow
  BC <- refDatasets$BC

  #get ranges for body weight and height correlated with age and sex
  rangeBW <- quantile(dat$BW, c(0.025, 0.975))  #get the 95% range of body weights
  rangeHT <- quantile(dat$HT/100, c(0.025, 0.975))  #get the 95% range of heights
  rangeBMI <- quantile(dat$BMI, c(0.025, 0.975))  #get the 95% range of body mass index
  rangeBSA <- quantile(dat$BSA, c(0.025, 0.975))  #get the 95% range of body surface area

  #check covariate range and throw an error when target measurements are out of range
  test_covRange(bw_targ, ht_targ, bmi_targ, rangeBW, rangeHT, rangeBMI)

  ##get mean bw, ht and bmi
  dat$LBW <- (M_bmi/(K_bmi + dat$BMI))*dat$BW
  bw_mean <- exp(mean(log(dat$BW)))  #geomteric mean for lognormally distributed weights in kg
  ht_mean <- mean(dat$HT)/100  #arithmetic mean for normally distributed heights in m
  bmi_mean <- exp(mean(log(dat$BMI)))  #geomteric mean for lognormally distributed BMI
  bsa_mean <- mean(dat$BSA)  #arithmetic mean for normally distributed BSA in m^2
  lbw_mean <- (M_bmi/(K_bmi + bmi_mean))*bw_mean  #exp(mean(log(dat$LBW)))  #geomteric mean for lognormally distributed lean body weights in kg

  ##get mean physiological parameters scaled by linear interpolation with height
  df2 <- df %>% select(-c(bw, ht, bsa))  #get rid of anthropometric measurements
  ht_refs <- df$ht/100  #get reference body heights from ICRP in m
  bw_refs <- df$bw
  bsa_refs <- df$bsa  # get reference body surface area from ICRP in m^2
  bmi_refs <- (df$bw)/(df$ht)^2
  lbw_refs <- df$bw - df$ad

  if(method == "Willmann"){
    linIntFns <- apply(df2, 2, FUN=function(x) approxfun(ht_refs, x, rule=2))  #linear interpolation functions with height; if outside range use closest value
    physPars <- lapply(linIntFns, FUN=function(x) x(ht_mean))  #get list of mean physiological parameters
  }else if(method == "Huisinga"){
    ages <- icrpData$age
    linIntFns <- apply(df, 2, FUN=function(x) approxfun(ages, x, rule=2))  #linear interpolation functions with lean body weight; if outside range use closest value; used df as we need ref covariates
    physPars <- lapply(linIntFns, FUN=function(x) x(age))  #get list of mean physiological parameters
    bw_ref <- physPars$bw
    lbw_ref <- bw_ref - physPars$ad
    bsa_ref <- physPars$bsa
    physPars <- physPars[4:length(physPars)]
  }

  blVol <- physPars$bl  #extract blood volume

  df_temp <- data.frame(organ=names(physPars), means=as.numeric(physPars))  #get scaling for blood and cardiac output
  blVol <- df_temp$means[df_temp$organ == "bl"]  #extract blood volume
  df_temp <- bind_rows(df_temp, data.frame(organ=c("ve","ar"), means=c(0, 0)))  #temporarily use whole blood volume for arterial and venous blood vols
  df_temp <- df_temp %>% dplyr::filter(organ != "bl")  #remove blood volume from df
  df_vols <- df_temp %>% filter(organ != "co")
  df_co <- df_temp %>% filter(organ == "co")

  ################################# getting organ volumes ##########################################
  #add blood content to vols
  df_vols <- left_join(df_vols, BC, by="organ")
  if(addBC){
    df_vols <- df_vols %>% dplyr::mutate(means = means + (bloodPerc*blVol/100))
  }else{
    df_vols <- df_vols %>% mutate(means = case_when(organ == "ve" ~ 0.705*blVol,
                                                    organ == "ar" ~ 0.295*blVol,
                                                    TRUE ~ means))
  }
  df_temp <- bind_rows(df_vols %>% select(organ, means), df_co)  #join with co again

  #--# Willmann #--#

  if(method == "Willmann"){
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
    df_opt <- df4 %>% filter(organ != "ad")
    df_ad <- df4 %>% filter(organ == "ad")

    pert <- runif(1)  #perturbation to parameters
    df_opt <- df_opt %>% mutate(pertMeans = ifelse(organ %in% normOrgan,
                                                   means_scaled + std_scaled*pert,
                                                   means_scaled*exp(log(std_scaled)*pert))) %>%
      mutate(pertMeans = ifelse(organ == "sk", pertMeans*sqrt((bw_targ/(bw_mean*ht_rel^2))), pertMeans))
    ad <- bw_targ - sum(df_opt$pertMeans) #compute adipose volume by subtracting current weight from target bw_targ

    optVols <- df_opt$pertMeans  #initial parameter values

    if(optimize){
      optMod <- suppressWarnings(newuoa(optVols, optimVolsFunc, df_opt=df_opt, df_ad=df_ad, bw_targ=bw_targ, bw_mean=bw_mean, ht_rel=ht_rel, normOrgan=normOrgan))
      df_opt <- df_opt %>% mutate(optimized=optMod$par)
    }else{
      df_opt <- df_opt %>% mutate(optimized=means_scaled)
    }

    ad <- bw_targ - sum(df_opt$optimized)

    test_realisticInd(ad)

    v_ov <- c(df_opt$optimized, ad)

    l_ov <- as.list(v_ov)  #get the final scaled organ volumes in a list
    names(l_ov) <- paste("V", c(as.character(df_opt$organ), "ad"), sep="")  #name the list


    ################################# getting organ blood flows #####################################
    co <- df3$means_scaled[df3$organ == "co"]
    flow2 <- flow %>% dplyr::mutate(bfs = (flowPerc*co/100)*60) %>% dplyr::filter(!organ %in% c("ve","ar","lu","li")) #get flow rates in L/h and remove lungs, liver, arterial and venous blood flows as they will be calculated within the model
    flow2 <- flow2 %>%
      dplyr::mutate(sds = 0.05*bfs,
                    bfs2 = rnorm(nrow(flow2), mean=bfs, sd=sds))
    l_bf <- as.list(flow2$bfs2)  #list of blood flows
    names(l_bf) <- paste("Q", flow2$organ, sep="")  #name the list
    flow_li <- flow2 %>% filter(organ %in% c("st","sm_int","la_int","sp","pa","ha"))
    l_bf <- c(l_bf, Qli=sum(flow_li$bfs2), Qlu=sum(flow2$bfs2))

#--# Huisinga #--#

  }else if(method == "Huisinga"){
    scale_ad <- (bw_targ - lbw_targ) / (bw_ref - lbw_ref)
    scale_br <- 1
    scale_sk <- bsa_targ/bsa_ref
    w_br <- scale_br*df_temp$means[df_temp$organ == "br"]
    w_sk <- scale_sk*df_temp$means[df_temp$organ == "sk"]
    w_ad <- scale_ad*df_temp$means[df_temp$organ == "ad"]
    scale_ti <- (lbw_targ - w_br - w_sk) / (lbw_ref - df_temp$means[df_temp$organ == "br"] - df_temp$means[df_temp$organ == "sk"])

    df_temp$density <- 1
    df_temp2 <- df_temp %>% mutate(density = case_when(organ == "ad" ~ 0.92,
                                                       organ == "bo" ~ 1.3,
                                                       TRUE ~ density)) %>%
      mutate(wts = case_when(organ == "br" ~ w_br,
                              organ == "sk" ~ w_sk,
                              organ == "ad" ~ w_ad,
                              TRUE ~ scale_ti*means))

    # organ volumes #
    df_co <- df_temp2 %>% filter(organ == "co")

    df_vols <- df_temp2 %>%
      mutate(vols = wts/density) %>%
      filter(organ != "co")

    v_ov <- df_vols$vols

    l_ov <- as.list(v_ov)  #get the final scaled organ volumes in a list
    names(l_ov) <- paste("V", as.character(df_vols$organ), sep="")  #name the list

    # organ flows
    co <- df_co$wts
    flow2 <- flow %>% dplyr::mutate(bfs = (flowPerc*co/100)*60) %>% dplyr::filter(!organ %in% c("ve","ar","lu","li")) #get flow rates in L/h and remove lungs, liver, arterial and venous blood flows as they will be calculated within the model
    l_bf <- as.list(flow2$bfs)  #list of blood flows
    names(l_bf) <- paste("Q", flow2$organ, sep="")  #name the list
    flow_li <- flow2 %>% filter(organ %in% c("st","sm_int","la_int","sp","pa","ha"))
    l_bf <- c(l_bf, Qli=sum(flow_li$bfs), Qlu=sum(flow2$bfs))
  }

  ##getting final parameter list
  pars <- c(AGE=age, SEX=sex, BW=bw_targ, HT=ht_targ, BMI=bmi_targ, BSA=bsa_targ, LBW=lbw_targ, l_ov, l_bf)
}


#####################

#' Generate population physiological parameters
#'
#' Takes in desired population demographics and generates the population physiological parameters
#'
#' @param nSubj Number of subjects
#' @param minAge Minimum age in years
#' @param maxAge Maximum age in years
#' @param femPerc Percentage of females
#' @param minBW Minimum body weight in kg
#' @param maxBW Maximum body weight in kg
#' @param minHT Minimum height in m
#' @param maxHT Maximum height in m
#' @param minBMI Minimum body mass index in kg/m^2
#' @param maxBMI Maximum body mass index in kg/m^2
#' @param optimize if TRUE, an optimization step is done; for "Willmann" method
#' @param addBC if TRUE, blood content will be added to each organ
#' @param method Method to use for generating physiological parameters; choices are "Willmann" and "Huisinga"
##' @return List of named lists with physiological parameters for each individual in the population
#' @export
genPop <- function(nSubj, minAge, maxAge, femPerc, minBW = NULL, maxBW = NULL, minHT = NULL, maxHT = NULL, minBMI = NULL, maxBMI = NULL, optimize=FALSE, addBC=TRUE, method="Willmann"){
  ##age is in years; weight is in kg; height is in m

  test_genPopInput(minBW, maxBW, minHT, maxHT, minBMI, maxBMI)

  nFemale <- floor((femPerc/100)*nSubj)
  nMale <- nSubj-nFemale

  pars_m <- sampleIndPars(nSubj=nMale, minAge, maxAge, is.male=TRUE, minBW, maxBW, minHT, maxHT, minBMI, maxBMI, optimize=optimize, addBC=addBC, mab=FALSE, method=method)
  pars_f <- sampleIndPars(nSubj=nFemale, minAge, maxAge, is.male=FALSE, minBW, maxBW, minHT, maxHT, minBMI, maxBMI, optimize=optimize, addBC=addBC, mab=FALSE, method=method)

  pars <- c(pars_m, pars_f)

  # add IDs
  pars2 <- lapply(1:length(pars), function(i) return(c(list(ID = i), pars[[i]])))

  return(pars2)
}

#################################

## mA model physiological parameters

#' Generate individual physiological parameters for the mAb model
#'
#' Takes in desired individual demographics and generates the individual physiological parameters for the mAb model
#'
#' @param age Age of the individual in years
#' @param is.male if TRUE, individual is male
#' @param bw_targ Target body weight in kg
#' @param ht_targ Target height in m
#' @param bmi_targ Target body mass index in kg/m^2
#' @return Named list with physiological parameters for the mAb model for the desired individual
#' @importFrom magrittr %>%
#' @importFrom dplyr filter select mutate case_when tibble
#' @export

genInd_mab <- function(age, is.male, bw_targ=NULL, ht_targ=NULL, bmi_targ=NULL){
  indPars <- genInd(age=age, is.male=is.male, bw_targ=bw_targ, ht_targ=ht_targ, bmi_targ=bmi_targ, optimize=FALSE, addBC=TRUE)
  foo <- tibble(param = names(indPars),
                value = as.numeric(indPars))
  # add stomach and gonads volumes and flows to the "Other" compartment as they are missing in the mAb model
  # the lymph node blood flow will also be added to "Other" compartment as this flow is ignored in the mAb model
  vols_add <- foo$value[foo$param == "Vst"] + foo$value[foo$param == "Vgo"]
  flows_add <- foo$value[foo$param == "Qst"] + foo$value[foo$param == "Qgo"] + foo$value[foo$param == "Qln"]
  foo$value[foo$param == "Vot"] <- foo$value[foo$param == "Vot"] + vols_add
  foo$value[foo$param == "Qot"] <- foo$value[foo$param == "Qot"] + flows_add
  BCbl <- (BC$bloodPerc[BC$organ == "ar"] + BC$bloodPerc[BC$organ == "ve"])/100
  Vbl <- (foo$value[foo$param == "Var"] + foo$value[foo$param == "Vve"])/BCbl  #get full body volume

  ## make necessary adjustments
  ### vols
  physvols <- mabShahData %>% select(Organ, Vtot, Vplas, Vbc, Vi, Ve, Vc)
  names(physvols) <- c("Organ","total", "plas", "bc", "int", "endo", "cell")

  physvols <- physvols %>%
    mutate(param = paste0("V", tolower(substr(Organ, start = 1, stop = 2))),
           param = case_when(param == "Vbl" ~ "Vbc",
                             param == "Vsi" ~ "Vsm_int",
                             param == "Vli" & Organ == "LI" ~ "Vla_int",
                             TRUE ~ param))

  foo_vols <- left_join(physvols, foo, by="param") %>%
    mutate(value = ifelse(param == "Vpl", (3/5.3)*Vbl, value),
           value = ifelse(param == "Vbc", (2.3/5.3)*Vbl, value),
           plas2 = (plas/total) * value,
           int2 = (int/total) * value) %>%
    select(Organ, Total=value, Plas=plas2, Int=int2)
  l_vol <- as.list(foo_vols$Total)
  names(l_vol) <- paste0("V_", foo_vols$Organ)
  l_int <- as.list(foo_vols$Int)
  names(l_int) <- paste0("V_", foo_vols$Organ, "_IS")
  l_plas <- as.list(foo_vols$Plas)
  names(l_plas) <- paste0("V_", foo_vols$Organ, "_V")

  l_volAll <- c(l_vol,l_int,l_plas)

  ### flows
  physflows <- mabShahData %>% select(Organ, Qplas, Qbc)
  names(physflows) <- c("Organ","plas","bc")

  physflows <- physflows %>%
    mutate(param = paste0("Q", tolower(substr(Organ, start = 1, stop = 2))),
           param = case_when(param == "Qbl" ~ "Qbc",
                             param == "Qsi" ~ "Qsm_int",
                             param == "Qli" & Organ == "LI" ~ "Qla_int",
                             param == "Qli" & Organ == "Liver" ~ "Qha",
                             TRUE ~ param))

  foo_flows <- left_join(physflows, foo, by="param")
  foo_flows2 <- foo_flows %>%
    # get just plasma flow that is required by the model
    mutate(PLQ = value * (plas/(plas + bc) * 1000),
           # sum LN plasma flow to Other
           plas = ifelse(Organ == "Other", plas + plas[Organ == "LN"], plas),
           # reassign PLQ_SI and PLQ_Bone
           PLQ = ifelse(Organ == "SI", PLQ[Organ=="Lung"] * (plas[Organ=="SI"] / plas[Organ=="Lung"]), PLQ),
           PLQ = ifelse(Organ == "Bone", PLQ[Organ=="Lung"] * (plas[Organ=="Bone"] / plas[Organ=="Lung"]), PLQ),
           PLQ = ifelse(Organ == "Kidney", PLQ[Organ=="Lung"] * (plas[Organ=="Kidney"] / plas[Organ=="Lung"]), PLQ),
           PLQ = ifelse(Organ == "Muscle", PLQ[Organ=="Lung"] * (plas[Organ=="Muscle"] / plas[Organ=="Lung"]), PLQ),
           PLQ = ifelse(Organ == "LI", PLQ[Organ=="Lung"] * (plas[Organ=="LI"] / plas[Organ=="Lung"]), PLQ)) %>%
    filter(!Organ %in% c("LN","Plasma","Bloodcells"))

  # balance blood flow
  PLQSum <- sum(foo_flows2$PLQ) - foo_flows2$PLQ[foo_flows2$Organ == "Lung"]
  foo_flows3 <- foo_flows2 %>%
    mutate(PLQ = ifelse(Organ == "Lung", PLQSum, PLQ),
           plas = plas/1000,
           PLQ = PLQ/1000) %>%
    select(Organ, PLQ)
  l_flows <- as.list(foo_flows3$PLQ)
  names(l_flows) <- paste0("PLQ_", foo_flows3$Organ)

  physPars <- c(SEX = indPars$SEX, BW = indPars$BW, HT = indPars$HT, BMI = indPars$BMI, l_volAll, l_flows)

  # remove entries with NA that won't be used in model like IS in plasma..
  physPars2 <- physPars[!is.na(physPars)]

  return(physPars2)
}

#################################

#' Generate population physiological parameters for the mAb model
#'
#' Takes in desired population demographics and generates the population physiological parameters for the mAb model
#'
#' @param nSubj Number of subjects
#' @param minAge Minimum age in years
#' @param maxAge Maximum age in years
#' @param femPerc Percentage of females
#' @param minBW Minimum body weight in kg
#' @param maxBW Maximum body weight in kg
#' @param minHT Minimum height in m
#' @param maxHT Maximum height in m
#' @param minBMI Minimum body mass index in kg/m^2
#' @param maxBMI Maximum body mass index in kg/m^2
#' @return List of named lists with physiological parameters for the mAb model for each individual in the population
#' @export
genPop_mab <- function(nSubj, minAge, maxAge, femPerc, minBW = NULL, maxBW = NULL, minHT = NULL, maxHT = NULL, minBMI = NULL, maxBMI = NULL){
  ##age is in years; weight is in kg; height is in m

  test_genPopInput(minBW, maxBW, minHT, maxHT, minBMI, maxBMI)

  nFemale <- (femPerc/100)*nSubj
  nMale <- nSubj-nFemale

  pars_m <- sampleIndPars(nSubj=nMale, minAge, maxAge, is.male=TRUE, minBW, maxBW, minHT, maxHT, minBMI, maxBMI, optimize=FALSE, addBC=FALSE, mab=TRUE)
  pars_f <- sampleIndPars(nSubj=nFemale, minAge, maxAge, is.male=FALSE, minBW, maxBW, minHT, maxHT, minBMI, maxBMI, optimize=FALSE, addBC=FALSE, mab=TRUE)

  pars <- c(pars_m, pars_f)

  # add IDs
  pars2 <- lapply(1:length(pars), function(i) return(c(list(ID = i), pars[[i]])))

  return(pars2)
}

#################################

#' Get the built-in virtual population data
#'
#' Takes in the model type and returns the built-in virtual population data for 1000 subjects
#'
#' @param model Model type; can be "general" or "mAb"
#' @details The covariates used to generate the virtual population are:
#' @details Number of subjects: 1000
#' @details Age: 20-80 years
#' @details Female percentage: 50 percent
#' @details Body weight: 50-100 kg
#' @details Height: 1.5-1.9 m
#' @return List of named lists with physiological parameters for the virtual population
#' @export
getVirtPop <- function(model="general"){
  if(model == "general"){
    virtPop <- virtPop_gen
  }else{
    virtPop <- virtPop_mab
  }
  return(virtPop)
}
