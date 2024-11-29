
library(survival)
library(survminer)

unict=unique(clintab2$cancer.type.abbreviation)
load("clintab2.RData")
load("exptab2.RData")
load("neibrho.RData")
source("driver_neighbour_functions.R")


neibsurv=list()
neibsurv$coef=matrix(0,nrow=nrow(neibrho),ncol=length(unict))
neibsurv$p=matrix(1,nrow=nrow(neibrho),ncol=length(unict))
for (i in 1:nrow(neibrho)){
  for (j in 1:length(unict)){
    out=survtest(neibrho$Group.1[i],unict[j],clintab2,exptab2)
    neibsurv$coef[i,j]=out$survcoef
    neibsurv$p[i,j]=out$survp
  }
  print(i)
}

save(neibsurv, file="neibsurv.RData")