
library(arrow)
library(ggplot2)
library(ggpubr)
source("bct_wct_functions.R")


#figure 2
mutdf=read_feather(file="./data/mutation.feather")
expdf=read_feather(file="./data/expression.feather")
maingraph=read.csv("./data/main_graph.csv",header=T,stringsAsFactors = F)
mtb=read.csv("./data/mutation_burden.csv",header=T,stringsAsFactors = F)

ctdtab=mkctdtab(mtb,mutdf[,1:3081])
ctmtbtab=mkctmtbtab(mtb)
ctetab=mkctetab(mtb,expdf[,1:15464])

plotdf=mkdf2plotSpearmanMB("DDX3X","PPIA",ctdtab,ctetab,ctmtbtab,1)

fig2A=ggplot(plotdf,aes(x=meanexp,y=mutfreqmb,colour=cancer_type))+
  geom_point(size=2.5)+
  theme_minimal() +
  xlab("PPIA expression")+
  ylab("DDX3X mutation frequency\n(TMB normalized)")+
  theme(legend.text=element_text(size=8),legend.key.size = unit(0.3, "lines"))+
  theme(legend.title = element_text(size=8))+
  labs(colour="Cancer type")

plotdf=mkdf2plotLM("WWTR1","CDK1",mtb,mutdf,expdf,ctdtab,1)

fig2B=ggplot(plotdf,aes(x=cancer_type,y=exp,colour=as.factor(mutated)))+
  geom_boxplot()+
  theme_minimal() +
  xlab("Cancer type")+
  ylab("CDK1 expression")+
  labs(colour="mutated\nWWTR1")+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))

#pvalue distribution plots
bct_obs=read_feather("./data/bct_obs.feather")
wct_obs=read_feather("./data/wct_obs.feather")

bctrpmat=read_feather("./data/bctrpmat.feather")
wctrpmat=read_feather("./data/wctrpmat.feather")

rhopvec=as.vector(as.matrix(bctrpmat))
coefpvec=as.vector(as.matrix(wctrpmat))

rhopdf1=data.frame(data="observed", p_value=bct_obs$rho_pval)
rhopdf2=data.frame(data="permuted", p_value=rhopvec)
rhopdf3=rbind(rhopdf1,rhopdf2)

fig2C=ggplot(rhopdf3, aes(x=p_value, color=data, fill=data)) + 
  geom_density(alpha=0.3)+
  xlab("BCT p-value")+
  ylab("Probability density")+
  theme_minimal() +
  theme(legend.title = element_blank())+
  theme(legend.position="inside")+
  theme(legend.position.inside=c(0.7,0.25))

coefpdf1=data.frame(data="observed", p_value=wct_obs$coef_p)
coefpdf2=data.frame(data="permuted", p_value=coefpvec)
coefpdf3=rbind(coefpdf1,coefpdf2)

fig2D=ggplot(coefpdf3, aes(x=p_value, color=data, fill=data)) + 
  geom_density(alpha=0.3) +
  xlab("WCT p-value")+
  ylab("Probability density")+
  theme_minimal() +
  theme(legend.title = element_blank())+
  theme(legend.position="inside")+
  theme(legend.position.inside=c(0.7,0.25))


ggarrange(fig2A, fig2B, fig2C, fig2D, 
          labels = c("A", "B", "C","D"),
          ncol = 2, nrow = 2)


#Figure 3
# figures with distribution of fraction of positives

bctneibtabs=read.csv("./data/bctneibtabs.csv",header=T,stringsAsFactors = F)
wctneibtabs=read.csv("./data/wctneibtabs.csv",header=T,stringsAsFactors = F)

bctdrivertab=read.csv("./data/bctdrivertab.csv",header=T,stringsAsFactors = F)
wctdrivertab=read.csv("./data/wctdrivertab.csv",header=T,stringsAsFactors = F)

load("./data/bctbyneib.RData")
load("./data/wctbyneib.RData")
load("./data/bctbydriv.RData")
load("./data/wctbydriv.RData")

posfracdf1=data.frame(type="observed",posfrac=bctneibtabs$posfrac[bctneibtabs$enriched==1])
posfracdf2=data.frame(type="permuted",posfrac=bctbyneib$sample)

posfracdf=rbind(posfracdf1,posfracdf2)

#figure 3A
fig3A=ggplot(data=posfracdf,mapping=aes(x=posfrac,fill=type)) + 
  geom_density(alpha=0.3,position = 'identity') + 
  xlab("Fraction of positive BCT DA") +
  ylab("Probability density")+
  theme_minimal() +
  theme(legend.title = element_blank())+
  theme(legend.position="inside")+
  theme(legend.position.inside=c(0.3,0.8))


posfracdf1=data.frame(type="observed",posfrac=bctdrivertab$posfrac[bctdrivertab$enriched==1])
posfracdf2=data.frame(type="permuted",posfrac=bctbydriv$sample)

posfracdf=rbind(posfracdf1,posfracdf2)

#figure 3C
fig3C=ggplot(data=posfracdf,mapping=aes(x=posfrac,fill=type)) + 
  geom_density(alpha=0.3,position = 'identity') + 
  xlab("Fraction of positive BCT NA") +
  ylab("Probability density")+
  theme_minimal() +
  theme(legend.title = element_blank())+
  theme(legend.position="inside")+
  theme(legend.position.inside=c(0.15,0.8))


posfracdf1=data.frame(type="observed",posfrac=wctneibtabs$posfrac[wctneibtabs$enriched==1])
posfracdf2=data.frame(type="permuted",posfrac=wctbyneib$sample)

posfracdf=rbind(posfracdf1,posfracdf2)

#figure 3b
fig3B=ggplot(data=posfracdf,mapping=aes(x=posfrac,fill=type)) + 
  geom_density(alpha=0.3,position = 'identity') + 
  xlab("Fraction of positive WCT DA") +
  ylab("Probability density")+
  theme_minimal() +
  theme(legend.title = element_blank())+
  theme(legend.position="inside")+
  theme(legend.position.inside=c(0.2,0.8))


posfracdf1=data.frame(type="observed",posfrac=wctdrivertab$posfrac[wctdrivertab$enriched==1])
posfracdf2=data.frame(type="permuted",posfrac=wctbydriv$sample)

posfracdf=rbind(posfracdf1,posfracdf2)

#figure 3D
fig3D=ggplot(data=posfracdf,mapping=aes(x=posfrac,fill=type)) + 
  geom_density(alpha=0.3,position = 'identity') + 
  xlab("Fraction of positive WCT NA") +
  ylab("Probability density")+
  theme_minimal() + 
  theme(legend.title = element_blank())+
  theme(legend.position="inside")+
  theme(legend.position.inside=c(0.2,0.8))

ggarrange(fig3A, fig3B, fig3C, fig3D, 
          labels = c("A", "B", "C","D"),
          ncol = 2, nrow = 2)

#figure 4

ctneibtabs=rbind(wctneibtabs,bctneibtabs)

#figure 4A
fig4A=ggplot(ctneibtabs[ctneibtabs$enriched==1,],mapping=aes(y=npos,x=as.factor(SApos),fill=type))+
  geom_boxplot()+
  theme_minimal() +
  scale_y_continuous(trans='log10') + 
  xlab("Unfavourable SA") + 
  ylab("Positive DA") + 
  labs(fill="") +
  scale_x_discrete(labels=c("0" = "0", "1" = "1","2" = "2", "3" = ">2"))

#Figure 4B
fig4B=ggplot(ctneibtabs[ctneibtabs$enriched==1,],mapping=aes(y=nneg,x=as.factor(SAneg),fill=type))+
  geom_boxplot()+
  theme_minimal() + 
  scale_y_continuous(trans='log10') +
  xlab("Favourable SA") + 
  ylab("Negative DA") + 
  labs(fill="") +
  scale_x_discrete(labels=c("0" = "0", "1" = "1","2" = ">1"))

DAsignvec=c("mixed","mixed","mixed","mixed","positive","positive","positive","positive","negative","negative","negative","negative")
SAsignvec=c("positive","negative","both","none","positive","negative","both","none","positive","negative","both","none")
wORdf=data.frame(DAsign=DAsignvec,SAsign=SAsignvec,OR=1,ORinf=1,ORsup=1,pval=1,stringsAsFactors = F)
bORdf=data.frame(DAsign=DAsignvec,SAsign=SAsignvec,OR=1,ORinf=1,ORsup=1,pval=1,stringsAsFactors = F)


wORdf[1,3:6]=logicOR((wctneibtabs$enriched==1 & wctneibtabs$DAsign=="mixed"),(wctneibtabs$SAsign=="positive"))
wORdf[2,3:6]=logicOR((wctneibtabs$enriched==1 & wctneibtabs$DAsign=="mixed"),(wctneibtabs$SAsign=="negative"))
wORdf[3,3:6]=logicOR((wctneibtabs$enriched==1 & wctneibtabs$DAsign=="mixed"),(wctneibtabs$SAsign=="both"))
wORdf[4,3:6]=logicOR((wctneibtabs$enriched==1 & wctneibtabs$DAsign=="mixed"),(wctneibtabs$SAsign=="none"))


wORdf[5,3:6]=logicOR((wctneibtabs$enriched==1 & wctneibtabs$DAsign=="positive"),(wctneibtabs$SAsign=="positive"))
wORdf[6,3:6]=logicOR((wctneibtabs$enriched==1 & wctneibtabs$DAsign=="positive"),(wctneibtabs$SAsign=="negative"))
wORdf[7,3:6]=logicOR((wctneibtabs$enriched==1 & wctneibtabs$DAsign=="positive"),(wctneibtabs$SAsign=="both"))
wORdf[8,3:6]=logicOR((wctneibtabs$enriched==1 & wctneibtabs$DAsign=="positive"),(wctneibtabs$SAsign=="none"))


wORdf[9,3:6]=logicOR((wctneibtabs$enriched==1 & wctneibtabs$DAsign=="negative"),(wctneibtabs$SAsign=="positive"))
wORdf[10,3:6]=logicOR((wctneibtabs$enriched==1 & wctneibtabs$DAsign=="negative"),(wctneibtabs$SAsign=="negative"))
wORdf[11,3:6]=logicOR((wctneibtabs$enriched==1 & wctneibtabs$DAsign=="negative"),(wctneibtabs$SAsign=="both"))
wORdf[12,3:6]=logicOR((wctneibtabs$enriched==1 & wctneibtabs$DAsign=="negative"),(wctneibtabs$SAsign=="none"))


bORdf[1,3:6]=logicOR((bctneibtabs$enriched==1 & bctneibtabs$DAsign=="mixed"),(bctneibtabs$SAsign=="positive"))
bORdf[2,3:6]=logicOR((bctneibtabs$enriched==1 & bctneibtabs$DAsign=="mixed"),(bctneibtabs$SAsign=="negative"))
bORdf[3,3:6]=logicOR((bctneibtabs$enriched==1 & bctneibtabs$DAsign=="mixed"),(bctneibtabs$SAsign=="both"))
bORdf[4,3:6]=logicOR((bctneibtabs$enriched==1 & bctneibtabs$DAsign=="mixed"),(bctneibtabs$SAsign=="none"))

bORdf[5,3:6]=logicOR((bctneibtabs$enriched==1 & bctneibtabs$DAsign=="positive"),(bctneibtabs$SAsign=="positive"))
bORdf[6,3:6]=logicOR((bctneibtabs$enriched==1 & bctneibtabs$DAsign=="positive"),(bctneibtabs$SAsign=="negative"))
bORdf[7,3:6]=logicOR((bctneibtabs$enriched==1 & bctneibtabs$DAsign=="positive"),(bctneibtabs$SAsign=="both"))
bORdf[8,3:6]=logicOR((bctneibtabs$enriched==1 & bctneibtabs$DAsign=="positive"),(bctneibtabs$SAsign=="none"))

bORdf[9,3:6]=logicOR((bctneibtabs$enriched==1 & bctneibtabs$DAsign=="negative"),(bctneibtabs$SAsign=="positive"))
bORdf[10,3:6]=logicOR((bctneibtabs$enriched==1 & bctneibtabs$DAsign=="negative"),(bctneibtabs$SAsign=="negative"))
bORdf[11,3:6]=logicOR((bctneibtabs$enriched==1 & bctneibtabs$DAsign=="negative"),(bctneibtabs$SAsign=="both"))
bORdf[12,3:6]=logicOR((bctneibtabs$enriched==1 & bctneibtabs$DAsign=="negative"),(bctneibtabs$SAsign=="none"))


wORdf$SAsign=factor(wORdf$SAsign,levels=c("positive", "negative", "both", "none"))
bORdf$SAsign=factor(bORdf$SAsign,levels=c("positive", "negative", "both", "none"))

#figure 4D
fig4D=ggplot(wORdf, aes(x = OR, y = DAsign, color = SAsign)) + 
  geom_vline(aes(xintercept = 1), size = .5, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = ORsup, xmin = ORinf), size = .75, height = .3, position = position_dodge(0.7)) +
  geom_point(size = 2, position = position_dodge(0.7)) +
  #coord_trans(x ="log10") +
  #scale_x_continuous(breaks = log10(seq(0.1, 2.5, 0.1)), labels = seq(0.1, 2.5, 0.1), limits = log10(c(0.05,3.5))) +
  theme_minimal()+
  theme(panel.grid.minor = element_blank()) +
  ylab("WCT DA pattern") +
  xlab("Odds ratio") +
  labs(color="SA") +
  annotate(geom = "text", y =1.04, x = 2.1, label = "*", size = 6, hjust = 0) +
  annotate(geom = "text", y =1.2, x = 0.15, label = "*", size = 6, hjust = 0) +
  scale_color_discrete(labels = c("Unfavourable", "Favourable", "Both", "None"))

# figure 4C
fig4C=ggplot(bORdf, aes(x = OR, y = DAsign, color = SAsign)) + 
  geom_vline(aes(xintercept = 1), size = .5, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = ORsup, xmin = ORinf), size = .75, height = .3, position = position_dodge(0.7)) +
  geom_point(size = 2, position = position_dodge(0.7)) +
  #coord_trans(x ="log10") +
  #scale_x_continuous(breaks = log10(seq(0.1, 2.5, 0.1)), labels = seq(0.1, 2.5, 0.1), limits = log10(c(0.05,3.5))) +
  theme_minimal()+
  theme(panel.grid.minor = element_blank()) +
  ylab("BCT DA pattern") +
  xlab("Odds ratio") +
  labs(color="SA") +
  annotate(geom = "text", y =1.67, x = 1.25, label = "***", size = 6, hjust = 0) +
  annotate(geom = "text", y =1.85 , x = 5.25, label = "***", size = 6, hjust = 0) +
  annotate(geom = "text", y =2.2, x = 1.25, label = "***", size = 6, hjust = 0) +
  annotate(geom = "text", y =2.66, x = 2.8, label = "***", size = 6, hjust = 0) +
  annotate(geom = "text", y =2.84, x = 1.25, label = "***", size = 6, hjust = 0) +
  annotate(geom = "text", y =3.2, x = 1.25, label = "***", size = 6, hjust = 0) +
  scale_color_discrete(labels = c("Unfavourable", "Favourable", "Both", "None"))

ggarrange(fig4A, fig4B, fig4C, fig4D, 
          labels = c("A", "B", "C","D"),
          ncol = 2, nrow = 2)


#####
#supplementary figures

#tau of bct enriched neighbours

tauhpa=read_feather("./data/consensus_grouped.feather")
tauhpa=tauhpa[,1:2]
tauhpa=tauhpa[!duplicated(tauhpa),]
names(tauhpa)=c("gene","tau")

wctneibtabstau=merge(wctneibtabs,tauhpa)
bctneibtabstau=merge(bctneibtabs,tauhpa)

ggplot(data=bctneibtabstau,mapping=aes(x=tau,fill=as.factor(enriched)))+
  geom_density(alpha=0.3)+
  labs(fill="BCT DA\nenrichment")+
  xlab("Tissue specificity (tau index)")+
  ylab("Probability density")+
  theme_minimal()+
  scale_fill_discrete(labels = c("No", "Yes"))+
  theme(legend.position="inside")+
  theme(legend.position.inside=c(0.8,0.8))
  

#z-score of oncogenes and tumour suppressor genes

ctdrivertab=rbind(wctdrivertab,bctdrivertab)

ggplot(ctdrivertab[ctdrivertab$enriched==1,],aes(x=role,y=z,fill=type))+
  geom_boxplot()+
  theme_minimal()+
  xlab("Driver role")+
  ylab("NA z-score")+
  labs(fill="")





