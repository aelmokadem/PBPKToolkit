
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

![](README_files/figure-gfm/unnamed-chunk-3-1.png)<!-- --> \# Monoclonal
antibody modeling

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

    ##                 [,1]
    ## SEX        1.0000000
    ## BW        73.0000000
    ## HT         1.7600000
    ## BMI       23.5666322
    ## V_Heart    0.3117982
    ## V_Lung     0.4728204
    ## V_Muscle  27.9873714
    ## V_Skin     3.0831851
    ## V_Adipose 16.5495677
    ## V_Bone    10.1870098

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
    ##  Min.   : 1.00   Min.   :1.0   Min.   :50.60   Min.   :1.503   Min.   :18.36  
    ##  1st Qu.:10.75   1st Qu.:1.0   1st Qu.:69.49   1st Qu.:1.594   1st Qu.:25.47  
    ##  Median :20.50   Median :1.5   Median :80.40   Median :1.698   Median :27.78  
    ##  Mean   :20.50   Mean   :1.5   Mean   :78.14   Mean   :1.676   Mean   :27.90  
    ##  3rd Qu.:30.25   3rd Qu.:2.0   3rd Qu.:85.56   3rd Qu.:1.731   3rd Qu.:29.94  
    ##  Max.   :40.00   Max.   :2.0   Max.   :98.16   Max.   :1.852   Max.   :37.76  
    ##     V_Heart           V_Lung          V_Muscle         V_Skin     
    ##  Min.   :0.2018   Min.   :0.2420   Min.   :16.03   Min.   :1.463  
    ##  1st Qu.:0.2217   1st Qu.:0.2970   1st Qu.:16.99   1st Qu.:1.750  
    ##  Median :0.2725   Median :0.4290   Median :22.28   Median :2.585  
    ##  Mean   :0.2711   Mean   :0.4001   Mean   :22.63   Mean   :2.481  
    ##  3rd Qu.:0.3189   3rd Qu.:0.4758   3rd Qu.:28.58   3rd Qu.:3.133  
    ##  Max.   :0.3427   Max.   :0.5378   Max.   :29.10   Max.   :3.577  
    ##    V_Adipose         V_Bone      
    ##  Min.   :10.33   Min.   : 6.179  
    ##  1st Qu.:23.96   1st Qu.: 7.056  
    ##  Median :28.09   Median : 9.015  
    ##  Mean   :30.03   Mean   : 8.701  
    ##  3rd Qu.:38.80   3rd Qu.: 9.851  
    ##  Max.   :57.44   Max.   :11.615
