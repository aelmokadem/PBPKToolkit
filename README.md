
# Load libraries

``` r
library(mrgPBPK)
library(dplyr)
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
df_Kps
```

    ## # A tibble: 1 x 11
    ##    Kpad  Kpbo  Kpbr  Kphe  Kpki  Kpgu  Kpli  Kplu  Kpmu  Kpsk  Kpsp
    ##   <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
    ## 1  3.24  4.56  4.06  5.89  2.81  3.54  3.01 0.785  1.31  2.67  1.62

# Generate individual physiological parameters

``` r
# define the individual demographics
## You only need to define two of the three covariates: body weight, height, and BMI
age <- 30
ismale <- T
bw <- 73
ht <- 1.76

# generate individual physiological parameters
indPars <- genInd(age=age, is.male=ismale, bw_targ=bw, ht_targ=ht, optimize = FALSE)
df_indPars <- bind_rows(indPars)
df_indPars
```

    ## # A tibble: 1 x 43
    ##     Vbo   Vbr    Vgo   Vhe   Vki Vla_int   Vli   Vpa   Vsk Vsm_int   Vst    Vth
    ##   <dbl> <dbl>  <dbl> <dbl> <dbl>   <dbl> <dbl> <dbl> <dbl>   <dbl> <dbl>  <dbl>
    ## 1  10.6  1.51 0.0335 0.367 0.410   0.479  2.26 0.168  3.25   0.836 0.200 0.0306
    ## # … with 31 more variables: Vln <dbl>, Vot <dbl>, Vve <dbl>, Var <dbl>,
    ## #   Vlu <dbl>, Vmu <dbl>, Vsp <dbl>, Vad <dbl>, Qbo <dbl>, Qbr <dbl>,
    ## #   Qgo <dbl>, Qhe <dbl>, Qki <dbl>, Qla_int <dbl>, Qpa <dbl>, Qsk <dbl>,
    ## #   Qsm_int <dbl>, Qst <dbl>, Qad <dbl>, Qmu <dbl>, Qsp <dbl>, Qha <dbl>,
    ## #   Qth <dbl>, Qln <dbl>, Qot <dbl>, Qli <dbl>, Qlu <dbl>, BW <dbl>, HT <dbl>,
    ## #   BMI <dbl>, SEX <dbl>

# Generate population physiological parameters

``` r
# define the population demographics
## You only need to define the range for two of the three covariates: body weight, height, and BMI
nSubj <- 10  #number of subjects
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
df_popPars
```

    ## # A tibble: 10 x 44
    ##       ID   Vbo   Vbr     Vgo   Vhe   Vki Vla_int   Vli   Vpa   Vsk Vsm_int   Vst
    ##    <int> <dbl> <dbl>   <dbl> <dbl> <dbl>   <dbl> <dbl> <dbl> <dbl>   <dbl> <dbl>
    ##  1     1 11.2   1.52 0.0370  0.387 0.424   0.496  2.37 0.174  3.51   0.868 0.207
    ##  2     2 10.7   1.51 0.0347  0.373 0.414   0.484  2.29 0.170  3.33   0.845 0.202
    ##  3     3  9.51  1.52 0.0347  0.363 0.399   0.467  2.23 0.164  3.08   0.816 0.195
    ##  4     4 11.4   1.51 0.0368  0.388 0.426   0.498  2.37 0.175  3.54   0.871 0.208
    ##  5     5  9.77  1.52 0.0358  0.371 0.405   0.474  2.27 0.167  3.18   0.828 0.198
    ##  6     6  7.23  1.34 0.00658 0.251 0.306   0.371  1.62 0.119  1.76   0.643 0.153
    ##  7     7  8.38  1.35 0.0103  0.285 0.349   0.435  1.80 0.140  2.32   0.738 0.176
    ##  8     8  8.38  1.35 0.00987 0.283 0.346   0.430  1.79 0.138  2.28   0.731 0.174
    ##  9     9  7.22  1.34 0.00844 0.263 0.320   0.396  1.67 0.127  1.94   0.676 0.161
    ## 10    10  8.63  1.35 0.0113  0.293 0.360   0.451  1.84 0.145  2.47   0.761 0.182
    ## # … with 32 more variables: Vth <dbl>, Vln <dbl>, Vot <dbl>, Vve <dbl>,
    ## #   Var <dbl>, Vlu <dbl>, Vmu <dbl>, Vsp <dbl>, Vad <dbl>, Qbo <dbl>,
    ## #   Qbr <dbl>, Qgo <dbl>, Qhe <dbl>, Qki <dbl>, Qla_int <dbl>, Qpa <dbl>,
    ## #   Qsk <dbl>, Qsm_int <dbl>, Qst <dbl>, Qad <dbl>, Qmu <dbl>, Qsp <dbl>,
    ## #   Qha <dbl>, Qth <dbl>, Qln <dbl>, Qot <dbl>, Qli <dbl>, Qlu <dbl>, BW <dbl>,
    ## #   HT <dbl>, BMI <dbl>, SEX <dbl>
