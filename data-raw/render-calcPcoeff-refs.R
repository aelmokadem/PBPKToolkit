devtools::load_all()
#' @importFrom dplyr bind_rows
#' @importFrom tibble tribble
#' @importFrom purrr pmap_dfr

## calcKp
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

## scaleKp
out_dir <- "calcKpRefs"
out_dir <- file.path(system.file("test-refs", package = "mrgPBPK"), out_dir)
if(!fs::dir_exists(out_dir)) fs::dir_create(out_dir)

logP <- 2  #lipophilicity
pKa <- 1  #acidic strength
type <- 3  #type of molecule
BP <- 1  #blood:plasma concentration ratio
fup <- 0.5  #unbound fraction in plasma
method <- "PT"  #prediction method
Kp <- calcKp(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type, method=method)

age <- 30
ismale <- TRUE
bw <- 73
ht <- 1.76

# generate individual physiological parameters
set.seed(123)
indPars <- genInd(age=age, is.male=ismale, bw_targ=bw, ht_targ=ht, optimize = FALSE)

res <- scaleKp(Kp=Kp, Vss=10, BP=BP, Vt=indPars)
dput(res, file.path(out_dir, "scaleKp"))

