logP <- 2  #lipophilicity
pKa <- 1  #acidic strength
type <- 3  #type of molecule
BP <- 1  #blood:plasma concentration ratio
fup <- 0.5  #unbound fraction in plasma
methods <- c("PT","RR","Beres","Schmitt","pksim")  #prediction method
outputFolder <- system.file("test-refs", "calcKpRefs", package = "mrgPBPK")

# calculate partition coefficients
for(method in methods){
  Kps <- calcKp(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type, method=method)
  dput(Kps, file=file.path(outputFolder, paste0(method, ".R")))
}

