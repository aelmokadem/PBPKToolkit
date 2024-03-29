out_dir <- "genPhysRefs"
out_dir <- file.path(system.file("test-refs", package = "PBPKToolkit"), out_dir)

## genInd_BC for Willmann
ref <- dget(file.path(out_dir, "genInd_BC_willmann"))

set.seed(123)
res <- genInd(age=30, is.male=TRUE, bw_targ=73, ht_targ=1.76, addBC=TRUE, optimize = TRUE, method="Willmann")

test_that("genInd", {
  expect_equal(
    res,
    ref
  )
})

## genInd_noBC for Willmann
ref <- dget(file.path(out_dir, "genInd_noBC_willmann"))

set.seed(123)
res <- genInd(age=30, is.male=TRUE, bw_targ=73, ht_targ=1.76, addBC=FALSE, optimize = TRUE, method="Willmann")

test_that("genInd", {
  expect_equal(
    res,
    ref
  )
})

## genInd_BC for Huisinga
ref <- dget(file.path(out_dir, "genInd_BC_huisinga"))

set.seed(123)
res <- genInd(age=30, is.male=TRUE, bw_targ=73, ht_targ=1.76, addBC=TRUE, method="Huisinga")

test_that("genInd", {
  expect_equal(
    res,
    ref
  )
})

## genInd_noBC for Huisinga
ref <- dget(file.path(out_dir, "genInd_noBC_huisinga"))

set.seed(123)
res <- genInd(age=30, is.male=TRUE, bw_targ=73, ht_targ=1.76, addBC=FALSE, method="Huisinga")

test_that("genInd", {
  expect_equal(
    res,
    ref
  )
})

test_that("genInd_error_outOfRange", {
  expect_error(
    genInd(age=30, is.male=TRUE, bw_targ=73, ht_targ=1.1),
    regexp = "Target height and BMI are out of range for the chosen age and sex"
  )
})

test_that("genInd_error_inputs", {
  expect_error(
    genInd(age=30, is.male=TRUE, bw_targ=73),
    regexp = "At least two inputs of bw_targ, ht_targ, or bmi_targ are required"
  )
})

## genPop for Willmann
ref <- dget(file.path(out_dir, "genPop_willmann"))

set.seed(123)
res <- genPop(nSubj=2, minAge=20, maxAge=80, femPerc=50, minBW=50, maxBW=100, minHT=1.5, maxHT=1.9, optimize=FALSE, method="Willmann")

test_that("genPop", {
  expect_equal(
    res,
    ref
  )
})

## genPop for Huisinga
ref <- dget(file.path(out_dir, "genPop_huisinga"))

set.seed(123)
res <- genPop(nSubj=2, minAge=20, maxAge=80, femPerc=50, minBW=50, maxBW=100, minHT=1.5, maxHT=1.9, optimize=FALSE, method="Huisinga")

test_that("genPop", {
  expect_equal(
    res,
    ref
  )
})

test_that("genPop_error_inputs", {
  expect_error(
    genPop(nSubj=2, minAge=20, maxAge=80, femPerc=50, minBW=50, maxBW=100, minHT=1.5),
    regexp = "Minimum and maximum values for height are required"
  )
})

## genInd_mab
ref <- dget(file.path(out_dir, "genInd_mab"))

set.seed(123)
res <- genInd_mab(age=30, is.male=TRUE, bw_targ=73, ht_targ=1.76)

test_that("genInd_mab", {
  expect_equal(
    res,
    ref
  )
})

## genPop_mab
ref <- dget(file.path(out_dir, "genPop_mab"))

set.seed(123)
res <- genPop_mab(nSubj=2, minAge=20, maxAge=80, femPerc=50, minBW=50, maxBW=100, minHT=1.5, maxHT=1.9)

test_that("genPop_mab", {
  expect_equal(
    res,
    ref
  )
})

