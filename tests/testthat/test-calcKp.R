library(testthat)

#ref_df <- readr::read_csv("calcKp-ref.csv", col_types = readr::cols())
#ref_df <- readRDS(system.file("test-refs","calcKp-ref.Rds", package = "mrgPBPK"))
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


# test_that("type13_error", {
#   expect_error(
#     calcKp(2, pKa=c(1,2), 0.5, 1, 1, "PT"),
#     regexp = "Molecule types 1, 2, and 3 require one pKa value"
#   )
# })


