
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
ht <- 1.7

# generate individual physiological parameters
indPars <- genInd(age=age, is.male=ismale, bw_targ=bw, ht_targ=ht)
df_indPars <- bind_rows(indPars)
df_indPars
```

    ## # A tibble: 1 x 39
    ##     Vbo   Vbr    Vgo   Vhe   Vki Vla_int   Vli   Vpa   Vsk Vsm_int   Vst    Vth
    ##   <dbl> <dbl>  <dbl> <dbl> <dbl>   <dbl> <dbl> <dbl> <dbl>   <dbl> <dbl>  <dbl>
    ## 1  9.87  1.51 0.0327 0.357 0.399   0.467  2.20 0.163  3.07   0.815 0.195 0.0298
    ## # … with 27 more variables: Vln <dbl>, Vve <dbl>, Var <dbl>, Vlu <dbl>,
    ## #   Vmu <dbl>, Vsp <dbl>, Vad <dbl>, Qbo <dbl>, Qbr <dbl>, Qgo <dbl>,
    ## #   Qhe <dbl>, Qki <dbl>, Qla_int <dbl>, Qpa <dbl>, Qsk <dbl>, Qsm_int <dbl>,
    ## #   Qst <dbl>, Qad <dbl>, Qmu <dbl>, Qsp <dbl>, Qha <dbl>, Qth <dbl>,
    ## #   Qln <dbl>, BW <dbl>, HT <dbl>, BMI <dbl>, SEX <dbl>

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
popPars <- genPop(nSubj=nSubj, minAge=minAge, maxAge=maxAge, femPerc=femPerc, minBW=minBW, maxBW=maxBW, minHT=minHT, maxHT=maxHT)
df_popPars <- bind_rows(popPars) %>% select(ID, everything())
df_popPars
```

    ## # A tibble: 10 x 40
    ##       ID   Vbo   Vbr     Vgo   Vhe   Vki Vla_int   Vli   Vpa   Vsk Vsm_int   Vst
    ##    <int> <dbl> <dbl>   <dbl> <dbl> <dbl>   <dbl> <dbl> <dbl> <dbl>   <dbl> <dbl>
    ##  1     1 10.0   1.51 0.0329  0.360 0.402   0.470  2.21 0.164  3.12   0.820 0.196
    ##  2     2 10.1   1.52 0.0362  0.375 0.411   0.480  2.30 0.169  3.27   0.840 0.200
    ##  3     3  9.66  1.51 0.0334  0.359 0.398   0.466  2.21 0.163  3.06   0.813 0.194
    ##  4     4  9.87  1.50 0.0296  0.344 0.392   0.459  2.13 0.160  2.95   0.800 0.191
    ##  5     5 10.0   1.52 0.0361  0.374 0.409   0.478  2.29 0.168  3.24   0.836 0.200
    ##  6     6  7.41  1.34 0.00668 0.254 0.309   0.376  1.65 0.121  1.80   0.652 0.155
    ##  7     7  6.60  1.33 0.00617 0.237 0.292   0.353  1.53 0.112  1.59   0.612 0.145
    ##  8     8  7.09  1.33 0.00635 0.244 0.300   0.363  1.58 0.115  1.69   0.629 0.149
    ##  9     9  8.55  1.35 0.0121  0.297 0.365   0.460  1.85 0.148  2.53   0.772 0.185
    ## 10    10  6.87  1.33 0.00636 0.243 0.298   0.361  1.58 0.115  1.67   0.626 0.149
    ## # … with 28 more variables: Vth <dbl>, Vln <dbl>, Vve <dbl>, Var <dbl>,
    ## #   Vlu <dbl>, Vmu <dbl>, Vsp <dbl>, Vad <dbl>, Qbo <dbl>, Qbr <dbl>,
    ## #   Qgo <dbl>, Qhe <dbl>, Qki <dbl>, Qla_int <dbl>, Qpa <dbl>, Qsk <dbl>,
    ## #   Qsm_int <dbl>, Qst <dbl>, Qad <dbl>, Qmu <dbl>, Qsp <dbl>, Qha <dbl>,
    ## #   Qth <dbl>, Qln <dbl>, BW <dbl>, HT <dbl>, BMI <dbl>, SEX <dbl>
