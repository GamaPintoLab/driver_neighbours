
#Add survival analysis results to neibrho and neibcoef dataframes

#run driver_neighbour_survival.R if neibsurv.RData is not available

load("neibsurv.RData")
load("neibrho.RData")
load("neibcoef.RData")
source("driver_neighbour_functions.R")

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

