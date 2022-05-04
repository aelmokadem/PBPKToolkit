out_dir <- "genPhysRefs"
out_dir <- file.path(system.file("test-refs", package = "PBPKToolkit"), out_dir)
if(!fs::dir_exists(out_dir)) fs::dir_create(out_dir)

## genInd with adding BC for Willmann method
set.seed(123)
res <- genInd(age=30, is.male=TRUE, bw_targ=73, ht_targ=1.76, optimize = TRUE, method="Willmann")
dput(res, file.path(out_dir, "genInd_BC_willmann"))

## genInd without BC for Willmann method
set.seed(123)
res <- genInd(age=30, is.male=TRUE, bw_targ=73, ht_targ=1.76, optimize = TRUE, addBC=FALSE, method="Willmann")
dput(res, file.path(out_dir, "genInd_noBC_willmann"))

## genInd with adding BC for Huisinga method
set.seed(123)
res <- genInd(age=30, is.male=TRUE, bw_targ=73, ht_targ=1.76, optimize = TRUE, method="Huisinga")
dput(res, file.path(out_dir, "genInd_BC_huisinga"))

## genInd without BC for Huisinga method
set.seed(123)
res <- genInd(age=30, is.male=TRUE, bw_targ=73, ht_targ=1.76, optimize = TRUE, addBC=FALSE, method="Huisinga")
dput(res, file.path(out_dir, "genInd_noBC_huisinga"))

## genPop for Willmann
set.seed(123)
res <- genPop(nSubj=2, minAge=20, maxAge=80, femPerc=50, minBW=50, maxBW=100, minHT=1.5, maxHT=1.9, optimize=FALSE, method="Willmann")
dput(res, file.path(out_dir, "genPop_willmann"))

## genPop for Huisinga
set.seed(123)
res <- genPop(nSubj=2, minAge=20, maxAge=80, femPerc=50, minBW=50, maxBW=100, minHT=1.5, maxHT=1.9, optimize=FALSE, method="Huisinga")
dput(res, file.path(out_dir, "genPop_huisinga"))

## genInd_mab
set.seed(123)
res <- genInd_mab(age=30, is.male=TRUE, bw_targ=73, ht_targ=1.76)
dput(res, file.path(out_dir, "genInd_mab"))

## genPop_mab
set.seed(123)
res <- genPop_mab(nSubj=2, minAge=20, maxAge=80, femPerc=50, minBW=50, maxBW=100, minHT=1.5, maxHT=1.9)
dput(res, file.path(out_dir, "genPop_mab"))
