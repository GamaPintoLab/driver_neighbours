library(arrow)
library(tictoc)
library(doParallel)
library(doSNOW)
nclust=7 #adjust to computer processor

#paired analysis

mutpairdf=read_feather(file="./data/processed/mutation_paired.feather")
exppairdf=read_feather(file="./data/processed/expression_paired.feather")
ctpair=mutpairdf[,1328:1329]

mgpaired=read.csv("./data/processed/main_graph_paired.csv",header=T, stringsAsFactors = F)

ctdtab1=mkctdtab(ctpair,mutpairdf[,1:1327])


cl <- makeCluster(nclust)
registerDoSNOW(cl)
tic()
dfrespair=foreach(i=1:nrow(mgpaired), .combine=rbind) %dopar% { #nrow(mgpaired)
  s=mkpairedtest(mgpaired$driver[i],mgpaired$neighbour[i],mutpairdf[,1:1329],exppairdf[,1:14245],ctpair,ctdtab1,minmut=1)
  return(s)
}
toc()
stopCluster(cl) 

pairedlm=cbind(mgpaired,dfrespair)
pairedlm$paircoef_p[is.na(pairedlm$paircoef_p)]=1


#compare paired with pathway analysis

fullpathanalysis=read.csv("./data/processed/full_pathway_analysis.csv",header=T,stringsAsFactors = F)

pairid=paste(fullpathanalysis$driver,fullpathanalysis$neighbour,sep="_")
fullpathanalysis$pairid=pairid

pairid=paste(pairedlm$driver,pairedlm$neighbour,sep="_")
pairedlm$pairid=pairid

pairedlm_comp=merge(pairedlm,fullpathanalysis,by="pairid",all.x=F,all.y=F)
pairedlm_comp$sig=(pairedlm_comp$paircoef_p<=0.05)*1

dist_sig=(table(pairedlm_comp[,c("distance","sig")]))
sigfrac=vector()
for (i in 1:nrow(dist_sig)){
  sigfrac[i]=dist_sig[i,2]/(dist_sig[i,1]+dist_sig[i,2])
}

plot(sigfrac)
sigfrac

psigfrac=vector()
for (i in 1:8){
  psigfrac[i]=phyper(dist_sig[i,2],nsig,ntot-nsig,dist_sig[i,2]+dist_sig[i,1],lower.tail=F)
}

sigdf=data.frame(SPL=rownames(dist_sig),NonSig=dist_sig[,1],Sig=dist_sig[,2],Freq=sigfrac,pval=psigfrac)

write.csv(sigdf,file="./data/processed/sigdf.csv",row.names=F)

pairs_to_exclude1=pairedlm[pairedlm$paircoef_p<=0.05,1:2]

pairs_to_exclude2=fullpathanalysis[which(fullpathanalysis$distance<=2),1:2]

pairstat=read_feather("./data/processed/pairstat.feather")

pairs_to_exclude3=pairstat[pairstat$nct<10,1:2]

pairs_to_exclude_wct=rbind(pairs_to_exclude1,pairs_to_exclude2)
pairs_to_exclude_wct=pairs_to_exclude_wct[!duplicated(pairs_to_exclude_wct),]

pairs_to_exclude_bct=rbind(pairs_to_exclude_wct,pairs_to_exclude3)
pairs_to_exclude_bct=pairs_to_exclude_bct[!duplicated(pairs_to_exclude_bct),]

write.csv(pairs_to_exclude_wct,file="./data/processed/pairs_to_exclude_wct.csv",row.names = F)
write.csv(pairs_to_exclude_bct,file="./data/processed/pairs_to_exclude_bct.csv",row.names = F)

