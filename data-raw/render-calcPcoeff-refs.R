devtools::load_all()
#' @importFrom dplyr bind_rows
#' @importFrom tibble tribble
#' @importFrom purrr pmap_dfr

test_cases <- tibble::tribble(
  ~logP,     ~pKa,   ~type,   ~BP,   ~fup,
  2,        1,       1,     1,    0.5,
  2,        1,       2,     1,    0.5,
  2,        1,       3,     1,    0.5,
  2,   c(1,2),       4,     1,    0.5,
  2,   c(1,2),       5,     1,    0.5,
  2,   c(1,2),       6,     1,    0.5,
  2, c(1,2,3),       7,     1,    0.5,
  2, c(1,2,3),       8,     1,    0.5,
  2, c(1,2,3),       9,     1,    0.5,
  2, c(1,2,3),      10,     1,    0.5
)
test_cases <- bind_rows(
  mutate(test_cases, method = "PT"),
  mutate(test_cases, method = "RR"),
  mutate(test_cases, method = "Berez"),
  mutate(test_cases, method = "Schmitt"),
  mutate(slice(test_cases, 1), method = "pksim")
)

res_df <- pmap_dfr(test_cases, render_ref, .func=calcKp, out_dir="calcKpRefs")
saveRDS(res_df, file.path(system.file("test-refs", package = "mrgPBPK"), "calcKp-ref.Rds"))

