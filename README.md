
# mrgPBPK

mrgPBPK is an open-source R package that provides a set of functions to
be used for PBPK modeling. The functions mainly generate drug- and
system-specific parameters required to build a PBPK model.

# Load libraries

``` r
library(mrgPBPK)
library(dplyr)
library(ggplot2)
library(cowplot)
library(GGally)
```

# Generate drug-specific parameters

## Calculate tissue:plasma partition coefficients (Kp)

The function `calcKp` can calculate the molecule’s Kp values for
different organs using one of five different calculation methods:

  - PT: Poulin and Theil
    ([source](https://pubmed.ncbi.nlm.nih.gov/11782904/)).
  - RR: Rodgers and Rowland
    ([source1](https://pubmed.ncbi.nlm.nih.gov/15858854/) and
    [source2](https://pubmed.ncbi.nlm.nih.gov/16639716/)).
  - Berez: Berezhkovskiy
    ([source](https://pubmed.ncbi.nlm.nih.gov/15124219/)).
  - Schmitt: Schimtt
    ([source](https://pubmed.ncbi.nlm.nih.gov/17981004/)).
  - pksim: PK-Sim
    ([source](https://www.tandfonline.com/doi/abs/10.1517/17425255.1.1.159)).

The function uses the unified tissue composition data reported in
<https://dmd.aspetjournals.org/content/48/10/903>

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

![](man/figures/README-calcKp-1.png)<!-- -->

## Calculate blood:plasma concentration ratio (BP)

In case BP parameter was missing, the function `calcBP` can be used to
calculate BP based on the methods reported
[here](https://pubmed.ncbi.nlm.nih.gov/20549836/). There are two methods
to chose from:

  - Method 1: uses molecule’s fup to calculate BP.
  - Method 2: uses molecule’s logP (or another measurement of
    lipophilicity like logD) to calculate BP.

Note: drug type = “total” is the default type and it uses the regression
coefficients calculated by fitting different molecule types (acids,
bases, and neutrals) together.

``` r
# in case BP parameter was missing, the function calcBP
calcBP(fup = fup, method = 1)
```

    . [1] 0.827023

``` r
calcBP(logP = logP, fup = fup, method = 2)
```

    . [1] 1.034007

## Calculate unbound fraction in plasma (fup)

In case fup parameter was missing, the function `calcFup` can be used to
calculate fup using the molecule’s logP (or another measurement of
lipophilicity like logD) based on the method reported
[here](https://pubmed.ncbi.nlm.nih.gov/20549836/).

``` r
# in case BP parameter was missing, the function calcBP
calcFup(logP = logP)
```

    . [1] 0.2931778

# Generate system-specific parameters

System-related (physiologic) parameters required for PBPK modeling are
the organ volumes and blood flows. There are two versions of the
functions that generate the physiologic parameters, a version to be
applied for general PBPK modeling of small molecules and another
specific for monoclonal antibody PBPK modeling.

## General small molecule PBPK modeling

Two functions `genInd` and `genPop` can be used to generate a virtual
individual or a population, respectively, for PBPK modeling of small
molecules. The algorithm used is adapted from Willmann et al 2007
[source](https://pubmed.ncbi.nlm.nih.gov/17431751/) and it maintains a
correlation between the physiologic parameters and the sampled
individual’s covariates to generate realistic individuals. The algorithm
is guided by these databases:

  - NHANES data for anthropometric data
    ([source](https://www.cdc.gov/nchs/nhanes/index.htm)).
  - ICRP data for physiologic parameters of typical individuals
    ([source](https://journals.sagepub.com/doi/pdf/10.1177/ANIB_32_3-4)).
  - Willmann et al 2007 data for organ variability
    ([source](https://pubmed.ncbi.nlm.nih.gov/17431751/)).

The function `genInd` takes age, sex, and two inputs of body weight,
height, and BMI. It returns a named list of the covariates and the
physiologic parameters with volumes prefixed with `V` and flows prefixed
with `Q`. The logical argument `optimize` states if the user wants to
run an additional optimization step on the generated organ volumes. This
step would give better parameter predictions but will take longer to
run.

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

![](man/figures/README-genInd-1.png)<!-- -->

The function `genPop` takes ranges of ages, and two ranges of body
weights, heights, or BMIs as well as the percentage of females in the
desired population. It returns a named list of lists with the upper
level being the individuals’ IDs and the nested list is each
individual’s covariates and physiologic parameters.

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

![](man/figures/README-genPop-1.png)<!-- -->

# Monoclonal antibody PBPK modeling

The functions `genInd_mab` and `genPop_mab` take the same inputs as
`genInd` and `genPop`, respectively, but they return physiologic
parameters to be used for monoclonal antibody (mAb) PBPK modeling. The
parameter names are consistent with the names used in the mAb PBPK model
developed
[here](https://github.com/metrumresearchgroup/bioPBPK/tree/main/mAb_bamlanivimab)
and was based on the model reported by Shah and Betts 2012
([source](https://pubmed.ncbi.nlm.nih.gov/22143261/)) and Jones et al
2019
([source](https://ascpt.onlinelibrary.wiley.com/doi/full/10.1002/psp4.12461)).

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

    .                [,1]
    . SEX        1.000000
    . BW        73.000000
    . HT         1.760000
    . BMI       23.566632
    . V_Heart    0.366579
    . V_Lung     1.152143
    . V_Muscle  28.748212
    . V_Skin     3.249019
    . V_Adipose 16.829509
    . V_Bone    10.575605

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

    .        ID             SEX            BW              HT             BMI       
    .  Min.   : 1.00   Min.   :1.0   Min.   :52.89   Min.   :1.518   Min.   :19.21  
    .  1st Qu.:10.75   1st Qu.:1.0   1st Qu.:66.01   1st Qu.:1.618   1st Qu.:22.87  
    .  Median :20.50   Median :1.5   Median :79.66   Median :1.689   Median :26.23  
    .  Mean   :20.50   Mean   :1.5   Mean   :76.91   Mean   :1.690   Mean   :26.99  
    .  3rd Qu.:30.25   3rd Qu.:2.0   3rd Qu.:84.35   3rd Qu.:1.751   3rd Qu.:30.69  
    .  Max.   :40.00   Max.   :2.0   Max.   :97.84   Max.   :1.891   Max.   :34.88  
    .     V_Heart           V_Lung          V_Muscle         V_Skin     
    .  Min.   :0.2408   Min.   :0.6920   Min.   :16.61   Min.   :1.653  
    .  1st Qu.:0.2606   1st Qu.:0.7514   1st Qu.:17.47   1st Qu.:1.897  
    .  Median :0.3094   Median :0.9923   Median :22.75   Median :2.589  
    .  Mean   :0.3180   Mean   :0.9802   Mean   :23.20   Mean   :2.652  
    .  3rd Qu.:0.3722   3rd Qu.:1.1643   3rd Qu.:29.07   3rd Qu.:3.303  
    .  Max.   :0.4072   Max.   :1.2535   Max.   :29.93   Max.   :3.887  
    .    V_Adipose          V_Bone      
    .  Min.   : 8.462   Min.   : 6.891  
    .  1st Qu.:19.919   1st Qu.: 7.672  
    .  Median :25.223   Median : 8.784  
    .  Mean   :28.862   Mean   : 9.151  
    .  3rd Qu.:39.408   3rd Qu.:10.515  
    .  Max.   :55.480   Max.   :12.562
