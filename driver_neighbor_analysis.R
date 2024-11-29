

load("dfrespair2.RData")
load("dfreslmd3.RData")

load("ctdtab3.RData")
load("ctetab3.RData")
load("cancertype3.RData")
load("mutationtab3.RData")
load("exptab3.RData")
source("driver_neighbour_functions.R")


dfreslmd3p=dfreslmd3[dfrespair2$paircoef_p>0.05,] #dataframe with all interactions but excluding interactions where the paired data associations were significative 

plotdf=mkdf2plotSpearman("ATG3","RAN",ctdtab3,ctetab3,1)

fig2A=ggplot(plotdf,aes(x=meanexp,y=mutrfreq,colour=cancer_type))+
  geom_point(size=3)+
  xlab("RAN expression")+
  ylab("ATG3 mutation frequency")+
  labs(colour="Cancer type")

plotdf=mkdf2plotLM("WWTR1","CDK1",cancertype3,mutationtab3,exptab3,ctdtab3,1)

fig2B=ggplot(plotdf,aes(x=cancer_type,y=exp,colour=as.factor(mutated)))+
  geom_boxplot()+
  xlab("Cancer Type")+
  ylab("CDK1 expression")+
  labs(colour="mutated\nWWTR1")+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))


#concatenate p values in 101 random permutations of expdata
rhopmat=matrix(data=1,nrow=nrow(dfreslmd3),ncol=101)
coefpmat=matrix(data=1,nrow=nrow(dfreslmd3),ncol=101)
for (i in 1:101){
  load(paste0("dfreslmd3_r",i,".RData"))
  dfreslmd3_r$rho[is.na(dfreslmd3_r$rho)]=0
  dfreslmd3_r$coef[is.na(dfreslmd3_r$coef)]=0
  rhopmat[,i]=dfreslmd3_r$rho_p
  coefpmat[,i]=dfreslmd3_r$coef_p
  
}


rhopvec=as.vector(rhopmat)
coefpvec=as.vector(coefpmat)

rhopdf1=data.frame(data="observed", p_value=dfreslmd3$rho_p)
rhopdf2=data.frame(data="permuted", p_value=rhopvec)
rhopdf3=rbind(rhopdf1,rhopdf2)

fig2C=ggplot(rhopdf3, aes(x=p_value, color=data, fill=data)) + 
  geom_density(alpha=0.3)+
  xlab("BCT p-value")+
  ylab("Probability density")+
  theme(legend.title = element_blank())+
  theme(legend.position="inside")+
  theme(legend.position.inside=c(0.7,0.2))

coefpdf1=data.frame(data="observed", p_value=dfreslmd3$coef_p)
coefpdf2=data.frame(data="permuted", p_value=coefpvec)
coefpdf3=rbind(coefpdf1,coefpdf2)

fig2D=ggplot(coefpdf3, aes(x=p_value, color=data, fill=data)) + 
  geom_density(alpha=0.3) +
  xlab("WCT p-value")+
  ylab("Probability density")+
  theme(legend.title = element_blank())+
  theme(legend.position="inside")+
  theme(legend.position.inside=c(0.7,0.2))

library(ggpubr)

ggarrange(fig2A, fig2B, fig2C, fig2D, 
          labels = c("A", "B", "C","D"),
          ncol = 2, nrow = 2)


#significant (p<=0.05) interactions in the 101 random permutations 
rhopdf=cbind(dfreslmd3[,5:6],as.data.frame(1*(rhopmat<=0.05)))
rhopdf=rhopdf[dfrespair2$paircoef_p>0.05,] # eliminate interactions with significant differences in the paired comparison

coefpdf=cbind(dfreslmd3[,5:6],as.data.frame(1*(coefpmat<=0.05)))
coefpdf=coefpdf[dfrespair2$paircoef_p>0.05,] # eliminate interactions with significant differences in the paired comparison


#significant interactions (p<=0.05) in the real data
drhopdf=cbind(dfreslmd3[,5:6],1*(dfreslmd3$rho_p<=0.05))
drhopdf=drhopdf[dfrespair2$paircoef_p>0.05,] # eliminate interactions with significant differences in the paired comparison
dcoefpdf=cbind(dfreslmd3[,5:6],1*(dfreslmd3$coef_p<=0.05))
dcoefpdf=dcoefpdf[dfrespair2$paircoef_p>0.05,] # eliminate interactions with significant differences in the paired comparison


#count significant interactions per neighbour
neibrho=aggregate(drhopdf[,3],by=list(drhopdf$neighbour),FUN="sum")
neibrhor=aggregate(rhopdf[,3:103],by=list(rhopdf$neighbour),FUN="sum")

neibcoef=aggregate(dcoefpdf[,3],by=list(dcoefpdf$neighbour),FUN="sum")
neibcoefr=aggregate(coefpdf[,3:103],by=list(coefpdf$neighbour),FUN="sum")


#compute z-score for the number of significant interactions

neibrho$n=0
neibcoef$n=0
neibrho$bin=0
neibcoef$bin=0
neibrho$sd=0
neibcoef$sd=0
for (i in 1:15071){
  neibrho$n[i]=sum(drhopdf$neighbour==neibrho$Group.1[i])
  neibcoef$n[i]=sum(dcoefpdf$neighbour==neibcoef$Group.1[i])
  neibrho$bin[i]=sum(neibrhor[i,2:102])/101
  neibcoef$bin[i]=sum(neibcoefr[i,2:102])/101
  neibrho$sd[i]=sd(neibrhor[i,2:102])
  neibcoef$sd[i]=sd(neibcoefr[i,2:102])
}
neibrho$prob=neibrho$bin/neibrho$n
neibcoef$prob=neibcoef$bin/neibcoef$n
neibrho$z=(neibrho$x-neibrho$bin)/neibrho$sd
neibcoef$z=(neibcoef$x-neibcoef$bin)/neibcoef$sd


#computing z-scores for the random permutations
neibrhorz=neibrhor
neibcoefrz=neibcoefr
for (i in 1:101){
  neibrhorz[,i+1]=(neibrhor[,i+1]-neibrho$bin)/neibrho$sd
  neibcoefrz[,i+1]=(neibcoefr[,i+1]-neibcoef$bin)/neibcoef$sd
}

#exclude drivers and neighbors that in the real data do not have at least 2 significant interactions

neibrho1=neibrho[neibrho$x>1,]
neibrhor1=neibrhor[neibrho$x>1,]
neibrhorz1=neibrhorz[neibrho$x>1,]

neibcoef1=neibcoef[neibcoef$x>1,]
neibcoefr1=neibcoefr[neibcoef$x>1,]
neibcoefrz1=neibcoefrz[neibcoef$x>1,]

# finding the pvalue thresholds for adequate fdr


rhothres=2
fdrneibrhoz=1
while (fdrneibrhoz>0.05){
  rhothres=rhothres+0.05
  fdrneibrhoz=mean(colSums(neibrhorz1[,2:102]>=rhothres))/sum(neibrho1$z>=rhothres)  
}
print(rhothres)
sum(neibrho1$z>=rhothres)

coefthres=2
fdrneibcoefz=1
while (fdrneibcoefz>0.05){
  coefthres=coefthres+0.05
  fdrneibcoefz=mean(colSums(neibcoefrz1[,2:102]>=coefthres))/sum(neibcoef1$z>=coefthres)
}
print(coefthres)
sum(neibcoef1$z>=coefthres)

#prepare dataframe with neighbour statistics


for (i in 1:nrow(neibrho)){
  lines=dfreslmd3p[dfreslmd3p$neighbour==neibrho$Group.1[i],]
  neibrho$pos[i]=sum(lines$rho>0 & lines$rho_p<=0.05)
  neibrho$neg[i]=sum(lines$rho<=0 & lines$rho_p<=0.05)
}  

for (i in 1:nrow(neibcoef)){
  lines=dfreslmd3p[dfreslmd3p$neighbour==neibcoef$Group.1[i],]
  neibcoef$pos[i]=sum(lines$coef>0 & lines$coef_p<=0.05)
  neibcoef$neg[i]=sum(lines$coef<=0 & lines$coef_p<=0.05)
}

neibrho$relx=neibrho$x/neibrho$n
neibrho$posfrac=neibrho$pos/neibrho$x
neibrho$sign="mixed"
neibrho$sign[neibrho$posfrac<=0.1 & neibrho$x>1]="negative"
neibrho$sign[neibrho$posfrac>=0.9 & neibrho$x>1]="positive"
neibrho$drivsig="non significant"
neibrho$drivsig[neibrho$z>=rhothres & neibrho$x>1]="significant"

neibcoef$relx=neibcoef$x/neibcoef$n
neibcoef$posfrac=neibcoef$pos/neibcoef$x
neibcoef$sign="mixed"
neibcoef$sign[neibcoef$posfrac<=0.1 & neibcoef$x>1]="negative"
neibcoef$sign[neibcoef$posfrac>=0.9 & neibcoef$x>1]="positive"
neibcoef$drivsig="non significant"
neibcoef$drivsig[neibcoef$z>=coefthres & neibcoef$x>1]="significant"

save(neibrho,file="neibrho.RData")

#Figure 3

fig3A=ggplot(neibrho[neibrho$x>=0,], aes(x=n, y=x,colour=drivsig)) +
  geom_point(alpha=0.5) +
  #scale_size(range = c(.1, 24), name="Population (M)") +
  #scale_fill_viridis(discrete=TRUE, guide=FALSE, option="A") +
  #scale_colour_gradient2(low="red",mid="white",high="blue",midpoint=0.5) + 
  #theme_ipsum() +
  #theme(legend.position="bottom") +
  ylab("BCT Driver associations") +
  xlab("Driver interactions") +
  labs(colour="Enrichment")+
  theme(legend.position="inside")+
  theme(legend.position.inside=c(0.8,0.2))

fig3B=ggplot(neibcoef[neibcoef$x>=0,], aes(x=n, y=x,colour=drivsig)) +
  geom_point(alpha=0.5) +
  #scale_size(range = c(.1, 24), name="Population (M)") +
  #scale_fill_viridis(discrete=TRUE, guide=FALSE, option="A") +
  #scale_colour_gradient2(low="red",mid="white",high="blue",midpoint=0.5) + 
  #theme_ipsum() +
  #theme(legend.position="bottom") +
  ylab("WCT Driver associations") +
  xlab("Driver interactions") +
  labs(colour="Enrichment")+
  theme(legend.position="inside")+
  theme(legend.position.inside=c(0.8,0.2))


posfracdf1=data.frame(type="BCT",posfrac=neibrho$posfrac[neibrho$drivsig=="significant"])
posfracdf2=data.frame(type="WCT",posfrac=neibcoef$posfrac[neibcoef$drivsig=="significant"])
posfracdf=rbind(posfracdf1,posfracdf2)

fig3C=ggplot(data=posfracdf,mapping=aes(x=posfrac,fill=type)) + 
  geom_density(alpha=0.5,position = 'identity') + 
  xlab("Fraction of positive associations") +
  ylab("Probability density")+
  theme(legend.title = element_blank())+
  theme(legend.position="inside")+
  theme(legend.position.inside=c(0.2,0.8))


sigposfrac1=data.frame(type="BCT",sign=neibrho$sign[neibrho$drivsig=="significant"],z=neibrho$z[neibrho$drivsig=="significant"])
sigposfrac2=data.frame(type="WCT",sign=neibcoef$sign[neibcoef$drivsig=="significant"],z=neibcoef$z[neibcoef$drivsig=="significant"])
sigposfrac=rbind(sigposfrac1,sigposfrac2)

fig3D=ggplot(sigposfrac,aes(x=type,y=z,fill=sign))+
  geom_boxplot()+
  xlab("Driver association type")+
  ylab("Enrichment z-score")+
  labs(fill="Association signs")+
  ylim(c(2.5,10))


ggarrange(fig3A, fig3B, fig3C, fig3D, 
          labels = c("A", "B", "C","D"),
          ncol = 2, nrow = 2)


#Add survival analysis results to neibrho and neibcoef dataframes

#run driver_neighbour_survival.R if neibsurv.RData is not available

load("neibsurv.RData")

for (i in 1:nrow(neibrho)){
    sig=which(neibsurv$p[i,]<=0.05)
    possurv=which(neibsurv$coef[i,]>0)
    negsurv=which(neibsurv$coef[i,]<=0)
    neibrho$psurv[i]=length(intersect(sig,possurv))
    neibrho$nsurv[i]=length(intersect(sig,negsurv))
}

neibcoef$psurv=neibrho$psurv
neibcoef$nsurv=neibrho$nsurv

neibrho$psurvp=pbinom(neibrho$psurv-1,32,0.025,lower.tail=F)
neibrho$nsurvp=pbinom(neibrho$nsurv-1,32,0.025,lower.tail=F)
neibrho$pnsurvp=pbinom(neibrho$psurv+neibrho$nsurv-1,32,0.05,lower.tail=F)

neibrho$psurvpadj=p.adjust(neibrho$psurvp,method="BY")
neibrho$nsurvpadj=p.adjust(neibrho$nsurvp,method="BY")
neibrho$pnsurvpadj=p.adjust(neibrho$pnsurvp,method="BY")

neibcoef$psurvp=neibrho$psurvp
neibcoef$nsurvp=neibrho$nsurvp
neibcoef$pnsurvp=neibrho$pnsurvp
neibcoef$psurvpadj=neibrho$psurvpadj
neibcoef$nsurvpadj=neibrho$nsurvpadj
neibcoef$pnsurvpadj=neibrho$pnsurvpadj

rownames(neibrho)=neibrho$Group.1
rownames(neibcoef)=neibcoef$Group.1


neibrho$survsig="non significant"
neibrho$survsig[neibrho$pnsurvpadj<=0.05]="significant"
neibrho$psurvsig="non significant"
neibrho$psurvsig[neibrho$psurvpadj<=0.05]="significant"
neibrho$nsurvsig="non significant"
neibrho$nsurvsig[neibrho$nsurvpadj<=0.05]="significant"
neibrho$survsign="none"
neibrho$survsign[neibrho$psurvsig=="significant" & neibrho$nsurvsig=="non significant"]="positive"
neibrho$survsign[neibrho$psurvsig=="non significant" & neibrho$nsurvsig=="significant"]="negative"
neibrho$survsign[neibrho$psurvsig=="significant" & neibrho$nsurvsig=="significant"]="mixed"
neibrho$survposfrac=neibrho$psurv/(neibrho$psurv+neibrho$nsurv)

neibcoef$survsig="non significant"
neibcoef$survsig[neibcoef$pnsurvpadj<=0.05]="significant"
neibcoef$psurvsig="non significant"
neibcoef$psurvsig[neibcoef$psurvpadj<=0.05]="significant"
neibcoef$nsurvsig="non significant"
neibcoef$nsurvsig[neibcoef$nsurvpadj<=0.05]="significant"
neibcoef$survsign="none"
neibcoef$survsign[neibcoef$psurvsig=="significant" & neibcoef$nsurvsig=="non significant"]="positive"
neibcoef$survsign[neibcoef$psurvsig=="non significant" & neibcoef$nsurvsig=="significant"]="negative"
neibcoef$survsign[neibcoef$psurvsig=="significant" & neibcoef$nsurvsig=="significant"]="mixed"
neibcoef$survposfrac=neibcoef$psurv/(neibcoef$psurv+neibcoef$nsurv)


#Figure 4

library(boot)
# Function to compute Spearman's correlation
spearman_cor <- function(data, indices) {
  resampled_data <- data[indices, ]
  cor(resampled_data[,1], resampled_data[,2], method = "spearman")
}

spdf=data.frame(type=c(rep(c("BCT"),4),rep(c("WCT"),4)),enriched=rep(c("No","Yes"),4),sign=rep(c("+","+","-","-"),2),rho=0,ciI=0,ciS=0)

for (i in 1:nrow(spdf)){
  if (spdf$type[i]=="BCT"){
    if (spdf$enriched[i]=="No"){
      data=neibrho[neibrho$drivsig=="non significant",]
    } else {
      data=neibrho[neibrho$drivsig=="significant",]
    }
  } else {
    if (spdf$enriched[i]=="No"){
      data=neibcoef[neibcoef$drivsig=="non significant",]
    } else {
      data=neibcoef[neibcoef$drivsig=="significant",]
    }
  }
  if (spdf$sign[i]=="+"){
    data=data[,c("pos","psurv")]
  } else {
    data=data[,c("neg","nsurv")]
  }   
  spdf$rho[i]=cor.test(data[,1],data[,2],method="spearman",exact=F)$estimate
  # Perform bootstrap
  set.seed(123) # For reproducibility
  bootstrap_results <- boot(data, statistic = spearman_cor, R = 1000)
  # Confidence interval (default percentile method)
  conf_interval <- boot.ci(bootstrap_results, type = "perc")
  spdf$ciI[i]=conf_interval$percent[4]
  spdf$ciS[i]=conf_interval$percent[5]
}

fig4A=ggplot(spdf,aes(x=interaction(type,sign),y=rho,fill=enriched))+
  geom_bar(stat = "identity", position = position_dodge())+
  geom_errorbar(aes(ymax = ciS, ymin = ciI),position = position_dodge(width = 0.9), width = 0.2)+
  ylab(expression(rho~"(#DA,#SA)"))+
  xlab("Type and sign of associations")+
  labs(fill="DA enriched")+
  scale_x_discrete(labels=c("BCT.-" = "BCT-", "WCT.-" = "WCT-","BCT.+" = "BCT+","WCT.+" = "WCT+"))

typevec=c("BCT","BCT","WCT","WCT")
drivvec=c("No","Yes","No","Yes") #,"NotPos","Pos","NotPos","Pos","NotNeg","Neg","NotNeg","Neg")
propvec=vector()
ciIvec=vector()
ciSvec=vector()


h=table(neibrho[,c("survsig","drivsig")])
h
p1=prop.test(h[1,2],sum(h[1,]))
p2=prop.test(h[2,2],sum(h[2,]))
propvec[1]=p1$estimate
propvec[2]=p2$estimate
ciIvec[1]=p1$conf.int[1]
ciSvec[1]=p1$conf.int[2]
ciIvec[2]=p2$conf.int[1]
ciSvec[2]=p2$conf.int[2]

fisher.test(h)

h=table(neibcoef[,c("survsig","drivsig")])
h
p1=prop.test(h[1,2],sum(h[1,]))
p2=prop.test(h[2,2],sum(h[2,]))
propvec[3]=p1$estimate
propvec[4]=p2$estimate
ciIvec[3]=p1$conf.int[1]
ciSvec[3]=p1$conf.int[2]
ciIvec[4]=p2$conf.int[1]
ciSvec[4]=p2$conf.int[2]
fisher.test(h)

propsigdf=data.frame(type=typevec,drivsig=drivvec,prop=propvec,ciI=ciIvec,ciS=ciSvec)


fig4B=ggplot(propsigdf,aes(x=type,y=prop,fill=drivsig))+
  geom_bar(stat = "identity", position = position_dodge())+
  geom_errorbar(aes(ymax = ciS, ymin = ciI),position = position_dodge(width = 0.9), width = 0.2)+
  xlab("Driver association type")+
  ylab("Proportion enriched in SA")+
  labs(fill="DA enriched")


typevec=c("BCT","BCT","WCT","WCT")
drivvec=c("No","Yes","No","Yes") #,"NotPos","Pos","NotPos","Pos","NotNeg","Neg","NotNeg","Neg")
propvec=vector()
ciIvec=vector()
ciSvec=vector()

h=table(data.frame(dpos=neibrho$sign=="positive",spos=neibrho$survsign=="positive"))  
h
p1=prop.test(h[1,2],sum(h[1,]))
p2=prop.test(h[2,2],sum(h[2,]))
propvec[1]=p1$estimate
propvec[2]=p2$estimate
ciIvec[1]=p1$conf.int[1]
ciSvec[1]=p1$conf.int[2]
ciIvec[2]=p2$conf.int[1]
ciSvec[2]=p2$conf.int[2]
fisher.test(h)

h=table(data.frame(dpos=neibcoef$sign=="positive",spos=neibcoef$survsign=="positive"))  
h
p1=prop.test(h[1,2],sum(h[1,]))
p2=prop.test(h[2,2],sum(h[2,]))
propvec[3]=p1$estimate
propvec[4]=p2$estimate
ciIvec[3]=p1$conf.int[1]
ciSvec[3]=p1$conf.int[2]
ciIvec[4]=p2$conf.int[1]
ciSvec[4]=p2$conf.int[2]

fisher.test(h)

propsignposdf=data.frame(type=typevec,drivsig=drivvec,prop=propvec,ciI=ciIvec,ciS=ciSvec)


fig4C=ggplot(propsignposdf,aes(x=type,y=prop,fill=drivsig))+
  geom_bar(stat = "identity", position = position_dodge())+
  geom_errorbar(aes(ymax = ciS, ymin = ciI),position = position_dodge(width = 0.9), width = 0.2)+
  xlab("Driver association type")+
  ylab("Proportion enriched in SA+")+
  labs(fill="DA+")


typevec=c("BCT","BCT","WCT","WCT")
drivvec=c("No","Yes","No","Yes") #,"NotPos","Pos","NotPos","Pos","NotNeg","Neg","NotNeg","Neg")
propvec=vector()
ciIvec=vector()
ciSvec=vector()

h=table(data.frame(dneg=neibrho$sign=="negative",spos=neibrho$survsign=="negative"))  
h
p1=prop.test(h[1,2],sum(h[1,]))
p2=prop.test(h[2,2],sum(h[2,]))
propvec[1]=p1$estimate
propvec[2]=p2$estimate
ciIvec[1]=p1$conf.int[1]
ciSvec[1]=p1$conf.int[2]
ciIvec[2]=p2$conf.int[1]
ciSvec[2]=p2$conf.int[2]
fisher.test(h)


h=table(data.frame(dneg=neibcoef$sign=="negative",spos=neibcoef$survsign=="negative"))  
h
p1=prop.test(h[1,2],sum(h[1,]))
p2=prop.test(h[2,2],sum(h[2,]))
propvec[3]=p1$estimate
propvec[4]=p2$estimate
ciIvec[3]=p1$conf.int[1]
ciSvec[3]=p1$conf.int[2]
ciIvec[4]=p2$conf.int[1]
ciSvec[4]=p2$conf.int[2]

fisher.test(h)

propsignnegdf=data.frame(type=typevec,drivsig=drivvec,prop=propvec,ciI=ciIvec,ciS=ciSvec)
propsignnegdf$drivsig <- factor(propsignnegdf$drivsig, levels = c("No", "Yes"))

fig4D=ggplot(propsignnegdf,aes(x=type,y=prop,fill=drivsig))+
  geom_bar(stat = "identity", position = position_dodge())+
  geom_errorbar(aes(ymax = ciS, ymin = ciI),position = position_dodge(width = 0.9), width = 0.2)+
  xlab("Driver association type")+
  ylab("Proportion enriched in SA-")+
  labs(fill="DA-")

fig4=ggarrange(fig4A, fig4B, fig4C, fig4D, 
          labels = c("A", "B", "C","D"),
          ncol = 2, nrow = 2)
annotate_figure(fig4,
                bottom="DA - Driver Association, SA - Survival Association")

# produce final tables with candidate neighbours

neibrhoPsel=neibrho[neibrho$drivsig=="significant" & neibrho$sign=="positive" & neibrho$survsign=="positive",]
neibrhoNsel=neibrho[neibrho$drivsig=="significant" & neibrho$sign=="negative" & neibrho$survsign=="negative",]

neibcoefPsel=neibcoef[neibcoef$drivsig=="significant" & neibcoef$sign=="positive" & neibcoef$survsign=="positive",]
neibcoefNsel=neibcoef[neibcoef$drivsig=="significant" & neibcoef$sign=="negative" & neibcoef$survsign=="negative",]

intersect(neibrhoPsel$Group.1,neibcoefPsel$Group.1)
intersect(neibrhoNsel$Group.1,neibcoefNsel$Group.1)

intersect(neibrhoPsel$Group.1,neibcoefNsel$Group.1)
intersect(neibrhoNsel$Group.1,neibcoefPsel$Group.1)

canonical=unique(cancerdrivers$symbol[cancerdrivers$type=="Canonical Cancer Driver"])
oncogenes=unique(cancerdrivers$symbol[cancerdrivers$NCG_oncogene==1])
tsg=unique(cancerdrivers$symbol[cancerdrivers$NCG_tsg==1])
knowndriver=unique(cancerdrivers$symbol)

neibrhoPsel$kdriver=1*is.element(neibrhoPsel$Group.1,knowndriver)
neibrhoNsel$kdriver=1*is.element(neibrhoNsel$Group.1,knowndriver)
neibcoefPsel$kdriver=1*is.element(neibcoefPsel$Group.1,knowndriver)
neibcoefNsel$kdriver=1*is.element(neibcoefNsel$Group.1,knowndriver)

neibrhoPsel$canonical=1*is.element(neibrhoPsel$Group.1,canonical)
neibrhoNsel$canonical=1*is.element(neibrhoNsel$Group.1,canonical)
neibcoefPsel$canonical=1*is.element(neibcoefPsel$Group.1,canonical)
neibcoefNsel$canonical=1*is.element(neibcoefNsel$Group.1,canonical)

neibrhoPsel$oncogene=1*is.element(neibrhoPsel$Group.1,oncogenes)
neibrhoNsel$oncogene=1*is.element(neibrhoNsel$Group.1,oncogenes)
neibcoefPsel$oncogene=1*is.element(neibcoefPsel$Group.1,oncogenes)
neibcoefNsel$oncogene=1*is.element(neibcoefNsel$Group.1,oncogenes)

neibrhoPsel$tsg=1*is.element(neibrhoPsel$Group.1,tsg)
neibrhoNsel$tsg=1*is.element(neibrhoNsel$Group.1,tsg)
neibcoefPsel$tsg=1*is.element(neibcoefPsel$Group.1,tsg)
neibcoefNsel$tsg=1*is.element(neibcoefNsel$Group.1,tsg)

neibrhoPsel$type="BCT"
neibrhoNsel$type="BCT"
neibcoefPsel$type="WCT"
neibcoefNsel$type="WCT"

suptab1=rbind(neibrhoPsel,neibcoefPsel)
suptab2=rbind(neibrhoNsel,neibcoefNsel)

write.csv(suptab1,file="supptable1.csv",row.names = F)
write.csv(suptab2,file="supptable2.csv",row.names = F)



