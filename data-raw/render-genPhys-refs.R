out_dir <- "genPhysRefs"
out_dir <- file.path(system.file("test-refs", package = "mrgPBPK"), out_dir)
if(!fs::dir_exists(out_dir)) fs::dir_create(out_dir)

## genInd with adding BC
set.seed(123)
res <- genInd(age=30, is.male=TRUE, bw_targ=73, ht_targ=1.76, optimize = TRUE)

dput(res, file.path(out_dir, "genInd_BC"))

## genInd without BC
set.seed(123)
res <- genInd(age=30, is.male=TRUE, bw_targ=73, ht_targ=1.76, optimize = TRUE, addBC=FALSE)

dput(res, file.path(out_dir, "genInd_noBC"))

## genPop
set.seed(123)
res <- genPop(nSubj=2, minAge=20, maxAge=80, femPerc=50, minBW=50, maxBW=100, minHT=1.5, maxHT=1.9, optimize=FALSE)

dput(res, file.path(out_dir, "genPop"))
