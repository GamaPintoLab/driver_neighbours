# libraries
library("arrow")
source("bct_wct_functions")
#read wct results

wct_obs=read_feather("./data/processed/wct_obs.feather")

#read shuffled p values

wctrpmat=matrix(data=1,nrow=nrow(wct_o),ncol=100)
wctrcmat=matrix(data=1,nrow=nrow(wct_o),ncol=100)
for (i in 1:100){
  load(paste0("./data/processed/shuffled_wct/wct_r",i,".RData"))
  wctrpmat[,i]=wct_r$coef_p
  wctrcmat[,i]=wct_r$coef
}
wctrsigmat=1*(wctrpmat<=0.05)
wctrsignmat=1*(wctrcmat>0)
wctrsigmat[is.na(wctrsigmat)]=0
wctrsignmat[is.na(wctrsignmat)]=0

write_feather(as.data.frame(wctrpmat),"./data/processed/wctrpmat.feather")

wctbyneib=xct_by_dn(wct_obs,wctrsigmat,wctrsignmat,"n","coef",0.1)

wctbydriv=xct_by_dn(wct_obs,wctrsigmat,wctrsignmat,"d","coef",0.1)


# read bct results

bct_obs=read_feather("./data/processed/bct_obs.feather")

bct_obs$pairid=paste(bct_obs$driver,bct_obs$neighbour,sep="_")

# filter out interactions were driver may regulate neighbour expression ou that have less than 10
# cancer types 
pairs_to_exclude_bct=read.csv("./data/processed/pairs_to_exclude_bct.csv",header=T,stringsAsFactors = F)

pairs_to_exclude_bct$pairid=paste(pairs_to_exclude_bct$driver,pairs_to_exclude_bct$neighbour,sep="_")

bct_obs=bct_obs[!is.element(bct_obs$pairid,pairs_to_exclude_bct$pairid),]

#read shuffle p-values
bctrpmat=matrix(data=1,nrow=nrow(bct_obs),ncol=100)
bctrcmat=matrix(data=0,nrow=nrow(bct_obs),ncol=100)
bct_o_match=bct_obs[,c(1,4)]
bct_o_match$pairid=paste(bct_o_match$driver,bct_o_match$neighbour,sep="_")
bct_o_match$order=(1:nrow(bct_o_match))

for (i in 1:100){
  filename=paste0("./data/processed/shuffled_bct/shuffle",i,".feather")
  bct_r=read_feather(filename)
  bct_r$pairid=paste(bct_r$driver,bct_r$neighbour,sep="_")
  bct_r_all=merge(bct_o_match,bct_r,by="pairid",all.x=T)
  bct_r_all$rho_pval[is.na(bct_r_all$rho_pval)]=1
  bct_r_all=bct_r_all[order(bct_r_all$order),]
  bctrpmat[,i]=bct_r_all$rho_pval
  bctrcmat[,i]=bct_r_all$rho
}

write_feather(as.data.frame(bctrpmat),"./data/processed/bctrpmat.feather")

bctrsigmat=1*(bctrpmat<=0.05)
bctrsigmat[is.na(bctrsigmat)]=0
bctrsignmat=1*(bctrcmat>0)
bctrsignmat[is.na(bctrsignmat)]=0


bctbyneib=xct_by_dn(bct_obs,bctrsigmat,bctrsignmat,"n","rho",0.1)

bctbydriv=xct_by_dn(bct_obs,bctrsigmat,bctrsignmat,"d","rho",0.1)

save(wctbyneib,file="./data/processed/wctbyneib.RData")
save(wctbydriv,file="./data/processed/wctbydriv.RData")
save(bctbyneib,file="./data/processed/bctbyneib.RData")
save(bctbydriv,file="./data/processed/bctbydriv.RData")

#non random distribution of significant associations

wneibn=length(unique(wct_o$neighbour))
wdrivn=length(unique(wct_o$driver))

bneibn=length(unique(bct_obs$neighbour))
bdrivn=length(unique(bct_obs$driver))

rdistwct=sigdist(wct_o,1000,2)
rdistbct=sigdist(bct_obs,1000,3)

mean(rdistwct$rnr)/wneibn
rdistwct$znr

mean(rdistbct$rnr)/bneibn
rdistbct$znr

mean(rdistwct$rndr)/wdrivn
rdistwct$zndr

mean(rdistbct$rndr)/bdrivn
rdistbct$zndr

rdf1=data.frame(analysis="WCT",type="Neighbours",fraction=rdistwct$rnr/wneibn,int=rdistwct$rmeann)
rdf2=data.frame(analysis="WCT",type="Drivers",fraction=rdistwct$rndr/wdrivn,int=rdistwct$rmeandn)

rdf3=data.frame(analysis="BCT",type="Neighbours",fraction=rdistbct$rnr/bneibn,int=rdistbct$rmeann)
rdf4=data.frame(analysis="BCT",type="Drivers",fraction=rdistbct$rndr/bdrivn, int=rdistbct$rmeandn)

rdfall=rbind(rdf1,rdf2)
rdfall=rbind(rdfall,rdf3)
rdfall=rbind(rdfall,rdf4)

obsdf=data.frame(analysis=c("WCT","WCT","BCT","BCT"),type=c("Neighbours","Drivers","Neighbours","Drivers"),fraction=0, int=0)
obsdf$fraction[1]=rdistwct$nr/wneibn
obsdf$fraction[2]=rdistwct$ndr/wdrivn
obsdf$fraction[3]=rdistbct$nr/bneibn
obsdf$fraction[4]=rdistbct$ndr/bdrivn
obsdf$int[1]=rdistwct$meann
obsdf$int[2]=rdistwct$meandn
obsdf$int[3]=rdistbct$meann
obsdf$int[4]=rdistbct$meandn


# include survival associations from human protein atlas

hpa_surv=read.delim("./data/raw/cancer_prognostic_data.tsv",header=T,stringsAsFactors = F)
hpa_surv_tcga=hpa_surv[grep("(TCGA)",hpa_surv$Cancer,fixed=T),]

survgenes=sort(unique(hpa_surv_tcga$Gene.name))

survtable=data.frame(gene=survgenes,posstrict=0,pos=0, negstrict=0, neg=0, allstrict=0, all=0, stringsAsFactors = F)

for (i in 1:length(survgenes)){
  lines=which(hpa_surv_tcga$Gene.name==survgenes[i])
  survtable$negstrict[i]=sum(!is.na(hpa_surv_tcga$potential.prognostic...favorable[lines[1:21]]))+sum(!is.na(hpa_surv_tcga$validated.prognostic...favorable[lines[1:21]]))
  survtable$neg[i]=survtable$posstrict[i]+sum(hpa_surv_tcga$unprognostic...favorable[lines[1:21]]<=0.05,na.rm=T)
  survtable$posstrict[i]=sum(!is.na(hpa_surv_tcga$potential.prognostic...unfavorable[lines[1:21]]))+sum(!is.na(hpa_surv_tcga$validated.prognostic...unfavorable[lines[1:21]]))
  survtable$pos[i]=survtable$posstrict[i]+sum(hpa_surv_tcga$unprognostic...unfavorable[lines[1:21]]<=0.05,na.rm=T)
  survtable$allstrict[i]=survtable$posstrict[i]+survtable$negstrict[i]
  survtable$all[i]=survtable$pos[i]+survtable$neg[i]
}

wctneibtab=data.frame(wctbyneib$stat,stringsAsFactors = F)
bctneibtab=data.frame(bctbyneib$stat,stringsAsFactors = F)

survneibtable=survtable[is.element(survtable$gene,wctneibtab$gene),]

wctneibtabs=merge(wctneibtab,survneibtable,by="gene",all.x=T)
wctneibtabs[is.na(wctneibtabs)]=0

bctneibtabs=merge(bctneibtab,survneibtable,by="gene",all.x=T)
bctneibtabs[is.na(bctneibtabs)]=0


wctneibtabs$SA=wctneibtabs$allstrict
wctneibtabs$SA[wctneibtabs$allstrict==6]=5
wctneibtabs$SA[wctneibtabs$allstrict==7]=5
wctneibtabs$SA[wctneibtabs$allstrict==8]=5

wctneibtabs$SApos=wctneibtabs$posstrict
wctneibtabs$SApos[wctneibtabs$posstrict==4]=3
wctneibtabs$SApos[wctneibtabs$posstrict==5]=3
wctneibtabs$SApos[wctneibtabs$posstrict==6]=3
wctneibtabs$SApos[wctneibtabs$posstrict==7]=3
wctneibtabs$SApos[wctneibtabs$posstrict==8]=3

wctneibtabs$SAneg=wctneibtabs$negstrict
wctneibtabs$SAneg[wctneibtabs$negstrict==3]=2
wctneibtabs$SAneg[wctneibtabs$negstrict==4]=2
wctneibtabs$SAneg[wctneibtabs$negstrict==5]=2
wctneibtabs$SAneg[wctneibtabs$negstrict==6]=2
wctneibtabs$SAneg[wctneibtabs$negstrict==7]=2
wctneibtabs$SAneg[wctneibtabs$negstrict==8]=2

bctneibtabs$SA=bctneibtabs$allstrict
bctneibtabs$SA[bctneibtabs$allstrict==6]=5
bctneibtabs$SA[bctneibtabs$allstrict==7]=5
bctneibtabs$SA[bctneibtabs$allstrict==8]=5

bctneibtabs$SApos=bctneibtabs$posstrict
bctneibtabs$SApos[bctneibtabs$posstrict==4]=3
bctneibtabs$SApos[bctneibtabs$posstrict==5]=3
bctneibtabs$SApos[bctneibtabs$posstrict==6]=3
bctneibtabs$SApos[bctneibtabs$posstrict==7]=3
bctneibtabs$SApos[bctneibtabs$posstrict==8]=3

bctneibtabs$SAneg=bctneibtabs$negstrict
bctneibtabs$SAneg[bctneibtabs$negstrict==3]=2
bctneibtabs$SAneg[bctneibtabs$negstrict==4]=2
bctneibtabs$SAneg[bctneibtabs$negstrict==5]=2
bctneibtabs$SAneg[bctneibtabs$negstrict==6]=2
bctneibtabs$SAneg[bctneibtabs$negstrict==7]=2
bctneibtabs$SAneg[bctneibtabs$negstrict==8]=2

wctneibtabs$type="WCT"
bctneibtabs$type="BCT"

ORwctneib=logicOR((wctneibtabs$enriched==1),(wctneibtabs$SA>1))

ORbctneib=logicOR((bctneibtabs$enriched==1),(bctneibtabs$SA>1))

# associations between DA sign patterns and SA sign patterns

wctneibtabs$DAsign="mixed"
wctneibtabs$DAsign[wctneibtabs$posfrac>=0.8]="positive"
wctneibtabs$DAsign[wctneibtabs$posfrac<=0.2]="negative"

bctneibtabs$DAsign="mixed"
bctneibtabs$DAsign[bctneibtabs$posfrac>=0.8]="positive"
bctneibtabs$DAsign[bctneibtabs$posfrac<=0.2]="negative"

wctneibtabs$SAsign="none"
wctneibtabs$SAsign[wctneibtabs$posstrict>0 & wctneibtabs$negstrict>0]="both"
wctneibtabs$SAsign[wctneibtabs$posstrict>0 & wctneibtabs$negstrict==0]="positive"
wctneibtabs$SAsign[wctneibtabs$posstrict==0 & wctneibtabs$negstrict>0]="negative"

bctneibtabs$SAsign="none"
bctneibtabs$SAsign[bctneibtabs$posstrict>0 & bctneibtabs$negstrict>0]="both"
bctneibtabs$SAsign[bctneibtabs$posstrict>0 & bctneibtabs$negstrict==0]="positive"
bctneibtabs$SAsign[bctneibtabs$posstrict==0 & bctneibtabs$negstrict>0]="negative"


#driver analysis

cancerdrivers <- read.delim("./data/raw/NCG_cancerdrivers_annotation_supporting_evidence.tsv",stringsAsFactors = F)

canonical=unique(cancerdrivers$symbol[cancerdrivers$type=="Canonical Cancer Driver"])
oncogenes=unique(cancerdrivers$symbol[cancerdrivers$NCG_oncogene==1])
oncogenes=oncogenes[!is.na(oncogenes)]
tsgenes=unique(cancerdrivers$symbol[cancerdrivers$NCG_tsg==1])
tsgenes=tsgenes[!is.na(tsgenes)]

bctdrivertab=bctbydriv$stat
wctdrivertab=wctbydriv$stat

bctdrivertab$canonical=0
bctdrivertab$canonical[is.element(bctdrivertab$gene,canonical)]=1
bctdrivertab$oncogenes=0
bctdrivertab$oncogenes[is.element(bctdrivertab$gene,oncogenes)]=1
bctdrivertab$tsgenes=0
bctdrivertab$tsgenes[is.element(bctdrivertab$gene,tsgenes)]=1


wctdrivertab$canonical=0
wctdrivertab$canonical[is.element(wctdrivertab$gene,canonical)]=1
wctdrivertab$oncogenes=0
wctdrivertab$oncogenes[is.element(wctdrivertab$gene,oncogenes)]=1
wctdrivertab$tsgenes=0
wctdrivertab$tsgenes[is.element(wctdrivertab$gene,tsgenes)]=1



bctdrivertab$role="unknown"
bctdrivertab$role[bctdrivertab$oncogenes==1]="OG"
bctdrivertab$role[bctdrivertab$tsgenes==1]="TSG"

wctdrivertab$role="unknown"
wctdrivertab$role[wctdrivertab$oncogenes==1]="OG"
wctdrivertab$role[wctdrivertab$tsgenes==1]="TSG"

wctdrivertab$type="WCT"
bctdrivertab$type="BCT"


save(wctdrivertab,file="./data/processed/wctdrivertab.RData")
save(bctdrivertab,file="./data/processed/bctdrivertab.RData")



#add driver info to neib dataframes

unidrivers=unique(cancerdrivers$symbol)
bctneibtabs$driver=0
bctneibtabs$driver[is.element(bctneibtabs$gene,unidrivers)]=1

bctneibtabs$canonical=0
bctneibtabs$canonical[is.element(bctneibtabs$gene,canonical)]=1

bctneibtabs$oncogene=0
bctneibtabs$oncogene[is.element(bctneibtabs$gene,oncogenes)]=1

bctneibtabs$tsgene=0
bctneibtabs$tsgene[is.element(bctneibtabs$gene,tsgenes)]=1

wctneibtabs$driver=0
wctneibtabs$driver[is.element(wctneibtabs$gene,unidrivers)]=1

wctneibtabs$canonical=0
wctneibtabs$canonical[is.element(wctneibtabs$gene,canonical)]=1

wctneibtabs$oncogene=0
wctneibtabs$oncogene[is.element(wctneibtabs$gene,oncogenes)]=1

wctneibtabs$tsgene=0
wctneibtabs$tsgene[is.element(wctneibtabs$gene,tsgenes)]=1

bctPcand=bctneibtabs[bctneibtabs$enriched==1 & bctneibtabs$DAsign=="positive" & bctneibtabs$SAsign=="positive" & bctneibtabs$posstrict>1,]
bctNcand=bctneibtabs[bctneibtabs$enriched==1 & bctneibtabs$DAsign=="negative" & bctneibtabs$SAsign=="negative" & bctneibtabs$negstrict>1,]

wctMcand=wctneibtabs[wctneibtabs$enriched==1 & wctneibtabs$DAsign=="mixed" & wctneibtabs$SAsign=="both" & wctneibtabs$SA>2,]
wctPcand=wctneibtabs[wctneibtabs$enriched==1 & wctneibtabs$DAsign=="positive" & wctneibtabs$SAsign=="positive",]
wctNcand=wctneibtabs[wctneibtabs$enriched==1 & wctneibtabs$DAsign=="negative" & wctneibtabs$SAsign=="negative",]


write.csv(wctneibtabs,file="./data/processed/wctneibtabs.csv",row.names=F)
write.csv(bctneibtabs,file="./data/processed/bctneibtabs.csv",row.names=F)

write.csv(wctdrivertab,file="./data/processed/wctdrivertab.csv",row.names=F)
write.csv(bctdrivertab,file="./data/processed/bctdrivertab.csv",row.names=F)
