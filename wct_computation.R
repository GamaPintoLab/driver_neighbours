#libraries
library(arrow)
library(tictoc)
library(doParallel)
library(doSNOW)
nclust=7 #adjust according to computer processor

#load data
mutdf=read_feather(file="./data/mutation.feather")
expdf=read_feather(file="./data/expression.feather")
maingraph=read.csv("./data/main_graph.csv",header=T,stringsAsFactors = F)
mtb=read.csv(".data/mutation_burden.csv",header=T,stringsAsFactors = F)

pairs_to_exclude_wct=read.csv("./data/pairs_to_exclude_wct.csv", header=T, stringsAsFactors = F)

pairid=paste(maingraph$driver,maingraph$neighbour,sep="_")
maingraph$pairid=pairid
wctmg=maingraph[which(!is.element(maingraph$pairid,pairs_to_exclude_wct$pairid)),]
wcttab=mkctdtab(mtb,mutdf[,1:3081])

unidrivers=unique(wctmg$driver)


#compute wct
tic()
cl <- makeCluster(nclust)
registerDoSNOW(cl)

wct_obs=foreach(i=1:length(unidrivers), .combine=rbind) %dopar% {
  s=mkwcttests(unidrivers[i],wctmg$neighbour[wctmg$driver==unidrivers[i]],mutdf[,1:3081],expdf[,1:15464],mtb,wcttab,minmut=1)  
  return(s)
}
stopCluster(cl) 
toc()

write_feather(wct_obs,"wct_obs.feather")

#shuffle

for (k in 1:100){
  
  rdata=mkrandexpwct(expdf[,1:15464])
  
  tic()
  cl <- makeCluster(nclust)
  registerDoSNOW(cl)
  
  wct_r=foreach(i=1:length(unidrivers), .combine=rbind) %dopar% {
    s=mkwcttests(unidrivers[i],wctmg$neighbour[wctmg$driver==unidrivers[i]],mutdf[,1:3081],rdata,mtb,wcttab,minmut=1)  
    return(s)
  }
  stopCluster(cl) 
  toc()
  
  save(wct_r,file=paste0("./data/shuffled_wct/wct_r",k,".RData"))
  
}



