devtools::load_all()

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
test_cases <- dplyr::bind_rows(
  dplyr::mutate(test_cases, method = "PT"),
  dplyr::mutate(test_cases, method = "RR"),
  dplyr::mutate(test_cases, method = "Berez"),
  dplyr::mutate(test_cases, method = "Schmitt"),
  dplyr::mutate(dplyr::slice(test_cases, 1), method = "pksim")
)

res_df <- purrr::pmap_dfr(test_cases, render_ref, .func=calcKp, out_dir="calcKpRefs")
saveRDS(res_df, file.path(system.file("test-refs", package = "PBPKToolkit"), "calcKp-ref.Rds"))

## calcKp with Vss input
out_dir <- "calcKpRefs"
out_dir <- file.path(system.file("test-refs", package = "PBPKToolkit"), out_dir)
if(!fs::dir_exists(out_dir)) fs::dir_create(out_dir)

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
dput(res, file.path(out_dir, "calcKp_Vss"))

