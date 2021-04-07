methods <- c("PT","RR","Beres","Schmitt","pksim")  #prediction method
refFolder <- system.file("test-refs", "calcKpRefs", package = "mrgPBPK")
logP <- 2  #lipophilicity
pKa <- 1  #acidic strength
type <- 3  #type of molecule
BP <- 1  #blood:plasma concentration ratio
fup <- 0.5  #unbound fraction in plasma

for(method in methods){
  test_that(paste0("calcKp for ", method), {
    Kps <- calcKp(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type, method=method)
    ref <- dget(file.path(refFolder, paste0(method, ".R")))

    expect_equal(Kps, ref)
  })
}


