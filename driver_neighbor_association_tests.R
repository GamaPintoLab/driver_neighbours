#perform association tests

driver_neib_pairs <- read.csv("driver_neib_pairs.csv")
load("mutationtab3.RData")
load("exptab3.RData")
load("cancertype3.RData")
load("dfrespair2.RData")
source("driver_neighbour_functions.R")


ctdtab3=mkctdtab(cancertype3,mutationtab3)
ctetab3=mkctetab(cancertype3,exptab3)

save(ctdtab3,file="ctdtab3.RData")
save(ctetab3,file="ctetab3.RData")


unidrivers=unique(driver_neib_pairs$driver)


library(doParallel)
library(progressr)
library(doSNOW)
handlers(global = TRUE)


cl <- makeCluster(7)
registerDoSNOW(cl) #registerDoParallel(cl)
unidrivers=unique(driver_neib_pairs$driver)
dfreslmd3=foreach(i=1:length(unidrivers), .combine=rbind) %dopar% {
  s=mktestslmd(unidrivers[i],driver_neib_pairs$neighbour[driver_neib_pairs$driver==unidrivers[i]],mutationtab3,exptab3,cancertype3,ctdtab3,ctetab3,minmut=1)  
  return(s)
}
stopCluster(cl) 


dfreslmd3$rho[is.na(dfreslmd3$rho)]=0
dfreslmd3$rho_p[is.na(dfreslmd3$rho_p)]=1
dfreslmd3$coef[is.na(dfreslmd3$coef)]=0
dfreslmd3$coef_p[is.na(dfreslmd3$coef_p)]=1


save(dfreslmd3,file="dfreslmd3.RData") #current version (made with mktestslmd, filter for mutations, ate least 3 ctypes with mutated individuals)

#perform association tests for randomly permutated datasets
for (k in 1:101){
  
  rdata3=mkrandexp(exptab3,cancertype3)
  
  tic()
  cl <- makeCluster(7)
  registerDoSNOW(cl) #registerDoParallel(cl)
  unidrivers=unique(driver_neib_pairs$driver)
  dfreslmd3_r=foreach(i=1:length(unidrivers), .combine=rbind) %dopar% {
    s=mktestslmd(unidrivers[i],driver_neib_pairs$neighbour[driver_neib_pairs$driver==unidrivers[i]],mutationtab3,rdata3$rexptab,cancertype3,ctdtab3,rdata3$rctetab,minmut=1)  
    return(s)
  }
  stopCluster(cl) 
  toc()
  
  save(dfreslmd3_r,file=paste0("dfreslmd3_r",k,".RData"))
  
}


