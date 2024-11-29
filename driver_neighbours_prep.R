#Preparation of necessary variables for analysis

driver_neib_pairs <- read.csv("driver_neib_pairs.csv")

cancerdrivers <- read.delim("NCG_cancerdrivers_annotation_supporting_evidence.tsv")

mutationtab=read.delim("mc3.v0.2.8.PUBLIC.nonsilentGene.xena",header=T,stringsAsFactors = F)
mutationtab$sample[7224]="XXX" #identify duplicate gene to eliminate

exptab=read.delim("RNASeqV2.geneExp.xena",header=T,stringsAsFactors = F)

which(duplicated(exptab$sample))
exptab$sample[16303]
exptab=exptab[!duplicated(exptab$sample),] #remove duplicated gene

for (i in 1:20530){ #substitute NaN values by minimum observed for the respective gene
  minexp=min(exptab[i,2:11070],na.rm=T)
  exptab[i,is.na(exptab[i,])]=minexp
  print(i)
}
save(exptab,file="exptab.RData")
load("exptab.RData")

clintab=read.delim("Survival_table",header=T,stringsAsFactors = F)
clintabsamples=gsub("-",".",clintab$sample,fixed=T) #to make patient identifiers uniform across data frames
clintab$sample=clintabsamples

phenotab=read.delim("TCGA_phenotype.tsv",header=T,stringsAsFactors = F)
phenosamples=gsub("-",".",phenotab$sample,fixed=T) #to make patient identifiers uniform across data frames
phenotab$sample=phenosamples


#variables for paired Tumour-Normal analysis

# identify pairs of tumour and healthy tissue samples from the same patient
phenotabT=phenotab[phenotab$sample_type!="Solid Tissue Normal",]

normalsamples=phenotab[phenotab$sample_type_id==11,]
expnames=names(exptab)

length(intersect(expnames,normalsamples$sample))

length(unique(expnames))

normal_exp=intersect(expnames,normalsamples$sample)
normal_exp_ind=substr(normal_exp,1,12)

exp_ind=substr(expnames,1,12)

normal_match=list()
normal_tab=data.frame(normid=rep('',737),tumourid=rep('',737))
for (i in 1:737){
  normal_match[[i]]=expnames[exp_ind==normal_exp_ind[i]]
  stype=substr(normal_match[[i]],14,15)
  a=which(stype=="11")
  b=which(stype=="01")
  if (length(a)==1 & length(b)==1){
    
    normal_tab[i,1]=normal_match[[i]][a]
    normal_tab[i,2]=normal_match[[i]][b]
  }
}

normal_tab=normal_tab[normal_tab$normid!='',]
normal_tab$individualid=substr(normal_tab$normid,1,12)

mutnames=names(mutationtab)

normal_tab=normal_tab[is.element(normal_tab$tumourid,mutnames),]

exptab0=exptab[,normal_tab$normid] #rna_seq expression table for healthy samples with paired tumour sample
rownames(exptab0)=exptab$sample

exptab1=exptab[,normal_tab$tumourid] #rna_seq expression table for corresponding tumor samples
rownames(exptab1)=exptab$sample

exptab1_0=exptab1-exptab0 #difference tumour - normal tissue 665 individuals

save(exptab1_0,file="exptab1_0.RData")


mutationtab1=mutationtab[,normal_tab$tumourid]
rownames(mutationtab1)=mutationtab$sample
unidrivers=unique(driver_neib_pairs$driver)
mutationtab1=mutationtab1[unidrivers,]

save(mutationtab1,file="mutationtab1.RData")

cancertype=clintab[,c(1,3,4,5,7)]
names(cancertype)=c("patient","cancer_type","age","gender","stage")

rownames(cancertype)=cancertype$patient
rownames(clintab)=cancertype$patient

cancertype1=cancertype[normal_tab$tumourid,] # cancer type info for paired comparison

save(cancertype1,file="cancertype1.RData")

#variables for BCT and WCT analyses

#get patient identifiers that are present in the mutation, expression, clinical data and are tumour samples
expnames=names(exptab)
mutnames=names(mutationtab)
samenames=intersect(mutnames,expnames)
samenames=intersect(samenames,phenotabT$sample)
samenames=union("sample",intersect(samenames,clintab$sample)) # tumour samples with exp, mutation, phenotype and clinical information

mutationtab2=mutationtab[,samenames]
rownames(mutationtab2)=mutationtab2$sample

exptab2=exptab[,samenames]
rownames(exptab2)=exptab2$sample
exptab2=exptab2[,2:8766]

cancertype2=cancertype[samenames[2:8766],]
clintab2=clintab[samenames[2:8766],]

unidrivers=unique(driver_neib_pairs$driver)
mutationtab2=mutationtab2[is.element(mutationtab2$sample,unidrivers),] #mutation table with driver genes only
rownames(mutationtab2)=mutationtab2$sample
mutationtab2=mutationtab2[,2:8766]

tmbd=colSums(mutationtab2) #sum of mutations in driver genes per patient
Q1=quantile(tmbd,0.25)
Q3=quantile(tmbd,0.75)
TF=Q3+1.5*(Q3-Q1)

mutationtab3=mutationtab2[,tmbd<TF] # filter out patientes with an outlier number of mutations in driver genes q3+1.5*(q3-q1)=85.5
exptab3=exptab2[,tmbd<TF] # filter for number of mutations
cancertype3=cancertype2[tmbd<TF,] # filter for number of mutations
clintab3=clintab2[tmbd<TF,]

save(mutationtab2,file="mutationtab2.RData")
save(mutationtab3,file="mutationtab3.RData")
save(exptab2,file="exptab2.RData")
save(exptab3,file="exptab3.RData")
save(cancertype2,file="cancertype2.RData")
save(cancertype3,file="cancertype3.RData")
save(clintab2,file="clintab2.RData")
save(clintab3,file="clintab3.RData")
