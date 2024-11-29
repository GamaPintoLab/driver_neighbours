#perform paired comparison tests

driver_neib_pairs <- read.csv("driver_neib_pairs.csv")
load("mutationtab1.RData")
load("exptab1_0.RData")
load("cancertype1.RData")

source("driver_neighbour_functions.R")
ctdtab1=mkctdtab(cancertype1,mutationtab1)


library(doParallel)
library(progressr)
library(doSNOW)
handlers(global = TRUE)


cl <- makeCluster(9)
registerDoSNOW(cl) #registerDoParallel(cl)
pb = txtProgressBar(min = 0, max = 362563, initial = 0) 
dfrespair=foreach(i=1:362563, .combine=rbind) %dopar% {
  s=mkpairedtest(driver_neib_pairs$driver[i],driver_neib_pairs$neighbour[i],mutationtab1,exptab1_0,cancertype1,ctdtab1,minmut=1)
  setTxtProgressBar(pb,i)
  return(s)
}
close(pb)
stopCluster(cl) 


dfrespair2=as.data.frame(dfrespair)

dfrespair2$paircoef_p[which(is.na(dfrespair2$paircoef_p))]=1


save(dfrespair2,file="dfrespair2.RData")
