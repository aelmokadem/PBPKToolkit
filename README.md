
# Load libraries

``` r
library(mrgPBPK)
library(dplyr)
library(ggplot2)
library(cowplot)
library(GGally)
```

# Calculate tissue:plasma partition coefficients

``` r
# define molecule's physicochemical properties
logP <- 2  #lipophilicity
pKa <- 1  #acidic strength
type <- 3  #type of molecule
BP <- 1  #blood:plasma concentration ratio
fup <- 0.5  #unbound fraction in plasma
method <- "PT"  #prediction method

# calculate partition coefficients
Kps <- calcKp(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type, method=method)
df_Kps <- bind_rows(Kps)
df_Kps2 <- tibble(Parameter=names(df_Kps), Value=as.numeric(Kps))

ggplot(data=df_Kps2, aes(Parameter, Value)) +
  geom_col() + theme_bw()
```

![](README_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

# General PBPK modeling

## Generate individual physiological parameters

``` r
# define the individual demographics
## You only need to define two of the three covariates: body weight, height, and BMI
age <- 30
ismale <- TRUE
bw <- 73
ht <- 1.76

# generate individual physiological parameters
indPars <- genInd(age=age, is.male=ismale, bw_targ=bw, ht_targ=ht, optimize = FALSE)
df_indPars <- bind_rows(indPars) %>% select(contains(c("Q","V")))
df_indPars2 <- tibble(Parameter=names(df_indPars), Value=as.numeric(df_indPars))

plot_vols <- ggplot(data=df_indPars2 %>% filter(grepl("V", Parameter)), aes(Parameter, Value)) + geom_col() + theme_bw() + labs(title = "Volumes")
plot_flows <- ggplot(data=df_indPars2 %>% filter(grepl("Q", Parameter)), aes(Parameter, Value)) + geom_col() + theme_bw() + labs(title = "Flows")

plot_grid(plot_vols, plot_flows, ncol=1)
```

![](README_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

## Generate population physiological parameters

``` r
# define the population demographics
## You only need to define the range for two of the three covariates: body weight, height, and BMI
nSubj <- 40  #number of subjects
minAge <- 20  #minimum age
maxAge <- 80  #maximum age
femPerc <- 50  #percentage of females
minBW <- 50  #minimum body weight
maxBW <- 100  #maximum body weight
minHT <- 1.5  #minimum height
maxHT <- 1.9  #maximum height

# generate population physiological parameters
popPars <- genPop(nSubj=nSubj, minAge=minAge, maxAge=maxAge, femPerc=femPerc, minBW=minBW, maxBW=maxBW, minHT=minHT, maxHT=maxHT, optimize = FALSE)
df_popPars <- bind_rows(popPars) %>% select(ID, everything())

df_popPars <- bind_rows(popPars) %>%
  mutate(SEX = ifelse(SEX == 1, "male", "female"),
         SEX = factor(SEX))
df_vols <- df_popPars %>% 
  select(AGE, HT, BW, SEX, contains("V"))

df_flows <- df_popPars %>% 
  select(AGE, HT, BW, SEX, contains("Q"))

plot_corr_vols <- ggcorr(df_vols %>% select(-SEX, -AGE, -BW, -HT, -Vot)) + labs(title = "Volumes")
plot_corr_flows <- ggcorr(df_flows %>% select(-SEX, -AGE, -BW, -HT, -Qot)) + labs(title = "Flows")
plot_grid(plot_corr_vols, plot_corr_flows, ncol=2)
```

![](README_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

# Monoclonal antibody modeling

## Generate individual physiological parameters

``` r
# define the individual demographics
## You only need to define two of the three covariates: body weight, height, and BMI
age <- 30
ismale <- T
bw <- 73
ht <- 1.76

# generate individual physiological parameters
indPars <- genInd_mab(age=age, is.male=ismale, bw_targ=bw, ht_targ=ht)
df_indPars <- bind_rows(indPars)

head(t(df_indPars), 10)
```

    ##                [,1]
    ## SEX        1.000000
    ## BW        73.000000
    ## HT         1.760000
    ## BMI       23.566632
    ## V_Heart    0.366579
    ## V_Lung     1.152143
    ## V_Muscle  28.748212
    ## V_Skin     3.249019
    ## V_Adipose 16.829509
    ## V_Bone    10.575605

## Generate population physiological parameters

``` r
# define the population demographics
## You only need to define the range for two of the three covariates: body weight, height, and BMI
nSubj <- 40  #number of subjects
minAge <- 20  #minimum age
maxAge <- 80  #maximum age
femPerc <- 50  #percentage of females
minBW <- 50  #minimum body weight
maxBW <- 100  #maximum body weight
minHT <- 1.5  #minimum height
maxHT <- 1.9  #maximum height

# generate population physiological parameters
popPars <- genPop_mab(nSubj=nSubj, minAge=minAge, maxAge=maxAge, femPerc=femPerc, minBW=minBW, maxBW=maxBW, minHT=minHT, maxHT=maxHT)
df_popPars <- bind_rows(popPars) %>% select(ID, everything())

summary(df_popPars[,1:11])
```

    ##        ID             SEX            BW              HT             BMI       
    ##  Min.   : 1.00   Min.   :1.0   Min.   :53.23   Min.   :1.508   Min.   :19.61  
    ##  1st Qu.:10.75   1st Qu.:1.0   1st Qu.:69.01   1st Qu.:1.603   1st Qu.:24.87  
    ##  Median :20.50   Median :1.5   Median :78.00   Median :1.661   Median :28.27  
    ##  Mean   :20.50   Mean   :1.5   Mean   :78.20   Mean   :1.667   Mean   :28.17  
    ##  3rd Qu.:30.25   3rd Qu.:2.0   3rd Qu.:89.49   3rd Qu.:1.735   3rd Qu.:30.64  
    ##  Max.   :40.00   Max.   :2.0   Max.   :97.25   Max.   :1.839   Max.   :38.75  
    ##     V_Heart           V_Lung          V_Muscle         V_Skin     
    ##  Min.   :0.2349   Min.   :0.6659   Min.   :16.59   Min.   :1.565  
    ##  1st Qu.:0.2496   1st Qu.:0.7218   1st Qu.:17.23   1st Qu.:1.775  
    ##  Median :0.3000   Median :0.9666   Median :21.98   Median :2.490  
    ##  Mean   :0.3090   Mean   :0.9488   Mean   :22.91   Mean   :2.524  
    ##  3rd Qu.:0.3658   3rd Qu.:1.1484   3rd Qu.:28.99   3rd Qu.:3.186  
    ##  Max.   :0.3976   Max.   :1.2297   Max.   :29.86   Max.   :3.694  
    ##    V_Adipose          V_Bone      
    ##  Min.   : 8.753   Min.   : 6.438  
    ##  1st Qu.:21.167   1st Qu.: 7.285  
    ##  Median :31.472   Median : 8.656  
    ##  Mean   :31.152   Mean   : 8.828  
    ##  3rd Qu.:37.717   3rd Qu.:10.243  
    ##  Max.   :56.791   Max.   :11.788
