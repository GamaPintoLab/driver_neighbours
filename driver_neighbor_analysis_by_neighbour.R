

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
save(neibcoef,file="neibcoef.RData")

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




