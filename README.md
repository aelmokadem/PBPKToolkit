
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
indPars <- genInd(age=age, is.male=ismale, bw_targ=bw, ht_targ=ht, optimize=TRUE)
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
popPars <- genPop(nSubj=nSubj, minAge=minAge, maxAge=maxAge, femPerc=femPerc, minBW=minBW, maxBW=maxBW, minHT=minHT, maxHT=maxHT)
```

    ## Error in test_covRange(bw_targ, ht_targ, bmi_targ, rangeBW, rangeHT, rangeBMI) : 
    ##   Target height is out of range for the chosen age and sex
    ## Error in test_covRange(bw_targ, ht_targ, bmi_targ, rangeBW, rangeHT, rangeBMI) : 
    ##   Target BMI is out of range for the chosen age and sex

``` r
df_popPars <- bind_rows(popPars) %>% select(ID, everything())
df_popPars
```

    ## # A tibble: 10 x 44
    ##       ID   Vbo   Vbr     Vgo   Vhe   Vki Vla_int   Vli   Vpa   Vsk Vsm_int   Vst
    ##    <int> <dbl> <dbl>   <dbl> <dbl> <dbl>   <dbl> <dbl> <dbl> <dbl>   <dbl> <dbl>
    ##  1     1 11.2   1.52 0.0369  0.387 0.425   0.496  2.37 0.174  3.51   0.868 0.207
    ##  2     2  9.65  1.50 0.0292  0.340 0.388   0.454  2.11 0.158  2.89   0.792 0.189
    ##  3     3 11.6   1.51 0.0367  0.388 0.428   0.501  2.38 0.176  3.58   0.875 0.209
    ##  4     4  9.74  1.52 0.0357  0.370 0.405   0.473  2.26 0.166  3.17   0.828 0.198
    ##  5     5 10.8   1.51 0.0347  0.374 0.416   0.486  2.30 0.170  3.35   0.849 0.203
    ##  6     6  7.00  1.34 0.00653 0.249 0.303   0.368  1.61 0.118  1.72   0.637 0.151
    ##  7     7  8.64  1.34 0.00709 0.270 0.328   0.399  1.75 0.128  2.04   0.691 0.164
    ##  8     8  7.08  1.34 0.00655 0.249 0.304   0.369  1.61 0.118  1.73   0.640 0.152
    ##  9     9  7.29  1.33 0.00640 0.246 0.303   0.366  1.59 0.116  1.72   0.635 0.151
    ## 10    10  7.47  1.34 0.00668 0.254 0.310   0.377  1.65 0.121  1.81   0.652 0.155
    ## # … with 32 more variables: Vth <dbl>, Vln <dbl>, Vot <dbl>, Vve <dbl>,
    ## #   Var <dbl>, Vlu <dbl>, Vmu <dbl>, Vsp <dbl>, Vad <dbl>, Qbo <dbl>,
    ## #   Qbr <dbl>, Qgo <dbl>, Qhe <dbl>, Qki <dbl>, Qla_int <dbl>, Qpa <dbl>,
    ## #   Qsk <dbl>, Qsm_int <dbl>, Qst <dbl>, Qad <dbl>, Qmu <dbl>, Qsp <dbl>,
    ## #   Qha <dbl>, Qth <dbl>, Qln <dbl>, Qot <dbl>, Qli <dbl>, Qlu <dbl>, BW <dbl>,
    ## #   HT <dbl>, BMI <dbl>, SEX <dbl>
