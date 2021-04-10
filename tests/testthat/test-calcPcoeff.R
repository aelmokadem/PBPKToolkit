library(testthat)

ref_df <- readRDS(file.path(REF_DIR, "calcKp-ref.Rds"))

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

test_that("types_1_warning", {
  expect_warning(
    calcKp(logP=2, pKa=c(1,2), fup=0.5, BP=1, type=1, method="PT"),
    regexp = "Molecule type 1 does not require pKa so it will be ignored"
  )
})

test_that("types_2_3_error", {
  expect_error(
    calcKp(logP=2, pKa=c(1,2), fup=0.5, BP=1, type=2, method="PT"),
    regexp = "Molecule types 2 and 3 require one pKa value"
  )
})

test_that("types_4_5_6_error", {
  expect_error(
    calcKp(logP=2, pKa=c(1), fup=0.5, BP=1, type=4, method="PT"),
    regexp = "Molecule types 4, 5, and 6 require two pKa values"
  )
})

test_that("types_7_8_9_10_error", {
  expect_error(
    calcKp(logP=2, pKa=c(1), fup=0.5, BP=1, type=7, method="PT"),
    regexp = "Molecule types 7, 8, 9, and 10 require three pKa values"
  )
})

