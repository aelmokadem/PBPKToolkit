
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
indPars <- genInd(age=age, is.male=ismale, bw_targ=bw, ht_targ=ht, optimize = TRUE)
df_indPars <- bind_rows(indPars)
df_indPars
```

    ## # A tibble: 1 x 43
    ##     Vbo   Vbr    Vgo   Vhe   Vki Vla_int   Vli   Vpa   Vsk Vsm_int   Vst    Vth
    ##   <dbl> <dbl>  <dbl> <dbl> <dbl>   <dbl> <dbl> <dbl> <dbl>   <dbl> <dbl>  <dbl>
    ## 1  10.6  1.51 0.0335 0.367 0.411   0.480  2.28 0.168  4.39   0.837 0.200 0.0306
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
```

    ## Error in test_covRange(bw_targ, ht_targ, bmi_targ, rangeBW, rangeHT, rangeBMI) : 
    ##   Target body weight and BMI are out of range for the chosen age and sex

``` r
df_popPars <- bind_rows(popPars) %>% select(ID, everything())
df_popPars
```

    ## # A tibble: 10 x 44
    ##       ID   Vbo   Vbr     Vgo   Vhe   Vki Vla_int   Vli   Vpa   Vsk Vsm_int   Vst
    ##    <int> <dbl> <dbl>   <dbl> <dbl> <dbl>   <dbl> <dbl> <dbl> <dbl>   <dbl> <dbl>
    ##  1     1 12.3   1.52 0.0390  0.404 0.442   0.516  2.47 0.182  3.83   0.903 0.216
    ##  2     2 10.9   1.52 0.0373  0.386 0.422   0.493  2.36 0.174  3.47   0.863 0.206
    ##  3     3 10.9   1.51 0.0333  0.368 0.413   0.483  2.27 0.169  3.30   0.842 0.201
    ##  4     4  9.97  1.50 0.0308  0.349 0.396   0.463  2.16 0.162  3.01   0.808 0.193
    ##  5     5 10.3   1.52 0.0364  0.378 0.413   0.483  2.31 0.170  3.31   0.844 0.202
    ##  6     6  7.74  1.35 0.0109  0.282 0.346   0.434  1.77 0.139  2.27   0.732 0.175
    ##  7     7  7.42  1.33 0.00646 0.247 0.305   0.369  1.60 0.117  1.75   0.640 0.152
    ##  8     8  8.04  1.35 0.0103  0.282 0.345   0.431  1.77 0.138  2.27   0.729 0.174
    ##  9     9  7.68  1.34 0.00844 0.268 0.326   0.402  1.71 0.129  2.02   0.689 0.164
    ## 10    10  6.52  1.34 0.00636 0.242 0.295   0.358  1.57 0.115  1.62   0.621 0.147
    ## # … with 32 more variables: Vth <dbl>, Vln <dbl>, Vot <dbl>, Vve <dbl>,
    ## #   Var <dbl>, Vlu <dbl>, Vmu <dbl>, Vsp <dbl>, Vad <dbl>, Qbo <dbl>,
    ## #   Qbr <dbl>, Qgo <dbl>, Qhe <dbl>, Qki <dbl>, Qla_int <dbl>, Qpa <dbl>,
    ## #   Qsk <dbl>, Qsm_int <dbl>, Qst <dbl>, Qad <dbl>, Qmu <dbl>, Qsp <dbl>,
    ## #   Qha <dbl>, Qth <dbl>, Qln <dbl>, Qot <dbl>, Qli <dbl>, Qlu <dbl>, BW <dbl>,
    ## #   HT <dbl>, BMI <dbl>, SEX <dbl>
