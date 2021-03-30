
# Load libraries

``` r
library(mrgPBPK)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
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
#> # A tibble: 1 x 11
#>    Kpad  Kpbo  Kpbr  Kphe  Kpki  Kpgu  Kpli  Kplu  Kpmu  Kpsk  Kpsp
#>   <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#> 1  3.24  4.56  4.06  5.89  2.81  3.54  3.01 0.785  1.31  2.67  1.62
```

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
#> # A tibble: 1 x 39
#>     Vbo   Vbr    Vgo   Vhe   Vki Vla_int   Vli   Vpa   Vsk Vsm_int   Vst    Vth
#>   <dbl> <dbl>  <dbl> <dbl> <dbl>   <dbl> <dbl> <dbl> <dbl>   <dbl> <dbl>  <dbl>
#> 1  9.87  1.51 0.0327 0.357 0.399   0.467  2.20 0.163  3.07   0.815 0.195 0.0298
#> # … with 27 more variables: Vln <dbl>, Vve <dbl>, Var <dbl>, Vlu <dbl>,
#> #   Vmu <dbl>, Vsp <dbl>, Vad <dbl>, Qbo <dbl>, Qbr <dbl>, Qgo <dbl>,
#> #   Qhe <dbl>, Qki <dbl>, Qla_int <dbl>, Qpa <dbl>, Qsk <dbl>, Qsm_int <dbl>,
#> #   Qst <dbl>, Qad <dbl>, Qmu <dbl>, Qsp <dbl>, Qha <dbl>, Qth <dbl>,
#> #   Qln <dbl>, BW <dbl>, HT <dbl>, BMI <dbl>, SEX <dbl>
```

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
#> Error in genInd(age = age, is.male = is.male, bw_targ = bw_targ, ht_targ = ht_targ,  : 
#>   Input values are out of range
#> Error in genInd(age = age, is.male = is.male, bw_targ = bw_targ, ht_targ = ht_targ,  : 
#>   Input values are out of range
df_popPars <- bind_rows(popPars) %>% select(ID, everything())
df_popPars
#> # A tibble: 10 x 40
#>       ID   Vbo   Vbr     Vgo   Vhe   Vki Vla_int   Vli   Vpa   Vsk Vsm_int   Vst
#>    <int> <dbl> <dbl>   <dbl> <dbl> <dbl>   <dbl> <dbl> <dbl> <dbl>   <dbl> <dbl>
#>  1     1  9.96  1.51 0.0344  0.366 0.404   0.473  2.25 0.166  3.16   0.826 0.197
#>  2     2 12.1   1.52 0.0387  0.401 0.439   0.513  2.45 0.180  3.77   0.897 0.214
#>  3     3  8.57  1.50 0.0282  0.326 0.372   0.435  2.02 0.152  2.64   0.759 0.181
#>  4     4  9.94  1.52 0.0360  0.373 0.408   0.477  2.28 0.168  3.22   0.834 0.199
#>  5     5 11.9   1.51 0.0363  0.389 0.431   0.504  2.39 0.177  3.63   0.881 0.210
#>  6     6  8.11  1.35 0.0110  0.287 0.351   0.441  1.79 0.142  2.35   0.743 0.178
#>  7     7  7.05  1.35 0.00953 0.267 0.326   0.406  1.68 0.131  2.01   0.690 0.165
#>  8     8  6.45  1.33 0.00614 0.235 0.290   0.350  1.52 0.111  1.57   0.608 0.144
#>  9     9  7.68  1.34 0.00722 0.261 0.317   0.387  1.68 0.124  1.90   0.669 0.159
#> 10    10  8.14  1.35 0.00985 0.280 0.343   0.426  1.77 0.137  2.24   0.724 0.173
#> # … with 28 more variables: Vth <dbl>, Vln <dbl>, Vve <dbl>, Var <dbl>,
#> #   Vlu <dbl>, Vmu <dbl>, Vsp <dbl>, Vad <dbl>, Qbo <dbl>, Qbr <dbl>,
#> #   Qgo <dbl>, Qhe <dbl>, Qki <dbl>, Qla_int <dbl>, Qpa <dbl>, Qsk <dbl>,
#> #   Qsm_int <dbl>, Qst <dbl>, Qad <dbl>, Qmu <dbl>, Qsp <dbl>, Qha <dbl>,
#> #   Qth <dbl>, Qln <dbl>, BW <dbl>, HT <dbl>, BMI <dbl>, SEX <dbl>
```
