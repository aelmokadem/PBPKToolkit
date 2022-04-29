## calcKp for different types and methods
ref_df <- readRDS(file.path(REF_DIR, "calcKp-ref.Rds"))
out_dir <- "calcKpRefs"
out_dir <- file.path(system.file("test-refs", package = "PBPKToolkit"), out_dir)

purrr::pwalk(ref_df, ~ {
  # capture all the columns for this row as a list
  .row <- list(...)

  # extract reference path from list
  ref_path <- .row$ref_path
  .row$ref_path <- NULL

  test_that(paste("calcKp works for", ref_path), {
    res <- do.call(calcKp, .row)
    ref <- dget(file.path(REF_DIR, "calcKpRefs", ref_path))
    expect_equal(res, ref)
  })
})

## calcKp type 1 warning
test_that("types_1_warning", {
  expect_warning(
    calcKp(logP=2, pKa=c(1,2), fup=0.5, BP=1, type=1, method="PT"),
    regexp = "Molecule type 1 does not require pKa so it will be ignored"
  )
})

## calcKp error for types 2 and 3
test_that("types_2_3_error", {
  expect_error(
    calcKp(logP=2, pKa=c(1,2), fup=0.5, BP=1, type=2, method="PT"),
    regexp = "Molecule types 2 and 3 require one pKa value"
  )
})

## calcKp error for types 4-6
test_that("types_4_5_6_error", {
  expect_error(
    calcKp(logP=2, pKa=c(1), fup=0.5, BP=1, type=4, method="PT"),
    regexp = "Molecule types 4, 5, and 6 require two pKa values"
  )
})

## calcKp error for types 7-10
test_that("types_7_8_9_10_error", {
  expect_error(
    calcKp(logP=2, pKa=c(1), fup=0.5, BP=1, type=7, method="PT"),
    regexp = "Molecule types 7, 8, 9, and 10 require three pKa values"
  )
})

## calcBP method=1
test_that("calcBP_method1", {
  expect_equal(
    calcBP(logP=1, fup=0.5, method=1),
    0.827023,
    tolerance = 6
  )
})

## calcBP method=2
test_that("calcBP_method2", {
  expect_equal(
    calcBP(logP=1, fup=0.5, type="total", method=2),
    0.953873,
    tolerance = 6
  )
})

## calcBP correct inputs
test_that("calcBPInput_error", {
  expect_error(
    calcBP(logP=NULL, fup=0.5, type="total", method=2),
    regexp = "logP is required for method=2"
  )
})

## calcFup
test_that("calcFup", {
  expect_equal(
    calcFup(logP=1, type="total"),
    0.3479642,
    tolerance = 7
  )
})

## calcKp with Vss input
ref <- dget(file.path(out_dir, "calcKp_Vss"))

logP <- 2  #lipophilicity
pKa <- 1  #acidic strength
type <- 3  #type of molecule
BP <- 1  #blood:plasma concentration ratio
fup <- 0.5  #unbound fraction in plasma
method <- "PT"  #prediction method
Vss <- 50

# generate individual physiological parameters to pass to Vt
age <- 30
ismale <- TRUE
bw <- 73
ht <- 1.76

set.seed(123)
indPars <- genInd(age=age, is.male=ismale, bw_targ=bw, ht_targ=ht, optimize = FALSE)

# get result and save
res <- calcKp(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type, method=method, Vss=Vss, Vt=indPars)

test_that("calcKp_Vss", {
  expect_equal(
    res,
    ref
  )
})

## calcVss
Kps <- calcKp(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type, method=method)

test_that("calcVss", {
  expect_equal(
    calcVss(Kp=Kps, BP=BP, Vt=indPars),
    186.9979,
    tolerance = 6
  )
})

