## script used to generate the built-in virtual PBPK populations
nSubj <- 1000  #number of subjects
minAge <- 20  #minimum age
maxAge <- 80  #maximum age
femPerc <- 50  #percentage of females
minBW <- 50  #minimum body weight
maxBW <- 100  #maximum body weight
minHT <- 1.5  #minimum height
maxHT <- 1.9  #maximum height

# generate population physiological parameters for small molecule PBPK modeling
set.seed(29837461)
virtPop_gen <- genPop(nSubj=nSubj, minAge=minAge, maxAge=maxAge, femPerc=femPerc, minBW=minBW, maxBW=maxBW, minHT=minHT, maxHT=maxHT, optimize = TRUE)

## generate population physiological parameters for mAb PBPK modeling
set.seed(29837461)
virtPop_mab <- genPop_mab(nSubj=nSubj, minAge=minAge, maxAge=maxAge, femPerc=femPerc, minBW=minBW, maxBW=maxBW, minHT=minHT, maxHT=maxHT)

