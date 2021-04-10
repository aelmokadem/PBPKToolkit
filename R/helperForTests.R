#' Render a ref from a test case row
#'
#' This takes a single row, runs the `.func` on it, and dput's the result to a
#' file. It also returns the same row with a `path` column added, which contains
#' the path the result was written to.
#' @param ... Dataframe with test cases in each row
#' @param .func Function to use the row inputs for
#' @param out_dir Output directory where the ref cases will be saved
#' @importFrom magrittr %>%
#' @importFrom dplyr as_tibble mutate
#' @keywords internal
#' **Meant to be called with `pmap_dfr()`**
render_ref <- function(..., .func, out_dir) {
  # capture all the columns for this row as a list
  .row <- list(...)
  # run the function
  res <- do.call(.func, .row)
  # write the result to a file
  out_path <- paste(
    paste(names(unlist(.row)), unlist(.row), sep = "="),
    collapse = "_"
  )

  out_dir <- file.path(system.file("test-refs", package = "mrgPBPK"), out_dir)
  if(!fs::dir_exists(out_dir)) fs::dir_create(out_dir)

  dput(res, file.path(out_dir, out_path))
  # return the same row, with the output path added
  .row %>%
    vec_to_list() %>%
    as_tibble() %>%
    mutate(ref_path = out_path)
}

#########################

#' Maps over a list and wraps any elements that are longer than 1 in list().
#' This is for converting lists to rows in a tibble with as_tibble()
#' @keywords internal
vec_to_list <- function(.row) {
  if ("pKa" %in% names(.row)) {
    .row$pKa <- list(.row$pKa)
  }
  return(.row)
}

