out_dir <- "genPhysRefs"
out_dir <- file.path(system.file("test-refs", package = "mrgPBPK"), out_dir)

## genInd
ref <- dget(file.path(out_dir, "genInd"))

set.seed(123)
res <- genInd(age=30, is.male=TRUE, bw_targ=73, ht_targ=1.76, optimize = TRUE)

test_that("genInd", {
  expect_equal(
    res,
    ref
  )
})

## genPop
ref <- dget(file.path(out_dir, "genPop"))

set.seed(123)
res <- genPop(nSubj=2, minAge=20, maxAge=80, femPerc=50, minBW=50, maxBW=100, minHT=1.5, maxHT=1.9, optimize=FALSE)

test_that("genInd", {
  expect_equal(
    res,
    ref
  )
})
