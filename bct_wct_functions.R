



###################

###################
mkctdtab=function(ctt,mutt){ #auxiliar function to compute fraction of mutated individuals per cancer type
  a=data.frame(table(ctt$cancer_type))
  names(a)=c("cancer_type","Freq")
  b=matrix(nrow=nrow(a),ncol=ncol(mutt))
  for (i in 1:nrow(a)){
    lines=which(ctt$cancer_type==a$cancer_type[i])
    if (length(lines)>1){
      b[i,]=colSums(as.matrix(mutt[lines,]))  
    } else {
      if (length(lines)==1){
        b[i,]=as.matrix(mutt[lines,])
      } else {
        b[i,]=rep(0,ncol(mutt))
      }
    }
    
  }
  b=as.data.frame(b)
  names(b)=names(mutt)
  a=cbind(a,b)
  a
}


###################

###################
mkctmtbtab=function(ctt){
  a=data.frame(table(ctt$cancer_type))
  names(a)=c("cancer_type","Freq")
  a$mtb=0
  for (i in 1:nrow(a)){
    a$mtb[i]=sum(ctt$mutation_burden[ctt$cancer_type==a$cancer_type[i]])/a$Freq[i]
  }
  a
}

###################

###################
mkctetab=function(ctt,expt){ #auxiliary function to compute average expression per cancer type
  a=data.frame(table(ctt$cancer_type))
  names(a)=c("cancer_type","Freq")
  b=matrix(nrow=nrow(a),ncol=ncol(expt))
  for (i in 1:nrow(a)){
    lines=which(ctt$cancer_type==a$cancer_type[i])
    if (length(lines)>1){
      b[i,]=colMeans(as.matrix(expt[lines,]))  
    } else {
      if (length(lines)==1){
        b[i,]=as.matrix(expt[lines,])
      } else {
        b[i,]=rep(0,ncol(mutt))
      }
    }
    
  }
  b=as.data.frame(b)
  names(b)=names(expt)
  a=cbind(a,b)
  a
}

###################

###################
mkpairedtest=function(dname,nname,mutt,expt,ctt,ctdtab,minmut){ #function to test if neighbour expression is influenced by driver mutation through the comparison of matched tumour and normal samples
  outt=data.frame(ctt,stringsAsFactors = F)
  outt$mutated=as.numeric(mutt[,dname])
  outt$exp=as.numeric(expt[,nname])
  a=ctdtab[,c(1,2)]
  a$mutfreq=ctdtab[,dname]
  a$mutrfreq=a$mutfreq/a$Freq
  b=a[a$mutfreq>=minmut,]
  if (nrow(b)>1){ #at least 2 cancer types with mutated individuals
    d=data.frame(outt[is.element(outt$cancer_type,b$cancer_type),])
    fitglm=lm(as.numeric(exp)~mutated+cancer_type,data=d)
    sumtab=summary(fitglm)
    output=c(paircoef=sumtab$coefficients["mutated",1],paircoef_p=sumtab$coefficients["mutated",4])  
  } else {
    output=c(paircoef=0,paircoef_p=1)
  }
  
}



###################

###################
mkwcttests=function(dname,nnames,mutt,expt,ctt,ctdtab,minmut){ #function to test if neighbour expression is associated with driver mutation status both within and between cancer types
  outt=data.frame(ctt,stringsAsFactors = F)
  outt$mutated=as.numeric(mutt[,dname])
  neibexp=(expt[,nnames])
  outt=cbind(outt,(neibexp))
  a=ctdtab[,c(1,2)]
  a$mutfreq=ctdtab[,dname]
  a$mutrfreq=a$mutfreq/a$Freq
  b=a[a$mutfreq>=minmut,]
  output=data.frame(coef=rep(0,length(nnames)),coef_p=1)
  if (nrow(b)>2){ #at least 3 cancer types with mutated individuals
    d=data.frame(outt[is.element(outt$cancer_type,b$cancer_type),])
    for (i in 1:length(nnames)){
      d$exp=d[,i+4]
      fitglm=lm(as.numeric(exp)~mutated+cancer_type+log10(mutation_burden+1),data=d)
      sumtab=summary(fitglm)
      output[i,]=c(sumtab$coefficients["mutated",1],sumtab$coefficients["mutated",4])
    }
  } 
  output$driver=dname
  output$neighbour=nnames
  output
}

###################

###################
mkrandexpwct=function(exptab){ #auxiliary function to generate random permutation of expression table
  rexptab=exptab[sample(nrow(exptab),nrow(exptab),replace=F),]
  names(rexptab)=names(exptab)
  rexptab
}


###################

###################
xct_by_dn=function(xct_o,xctrsigmat,xctrsignmat,type,tt,maxfdr){
  if (type=="d"){
    unitab=as.data.frame(table(xct_o$driver))
  } else {
    unitab=as.data.frame(table(xct_o$neighbour))
  }
  unitab=unitab[unitab$Freq>1,]
  xctneibstat=data.frame(gene=unitab[,1],nint=0,nsig=0,z=0,npos=0,nneg=0,posfrac=0,stringsAsFactors = F)
  xctrneibcount=matrix(data=0,nrow=nrow(unitab),ncol=ncol(xctrsigmat))
  xctrposcount=matrix(data=0,nrow=nrow(unitab),ncol=ncol(xctrsignmat))
  xctrsignmat=xctrsignmat*xctrsigmat
  
  for (i in 1:nrow(unitab)){
    if (type=="d"){
      lines=which(xct_o$driver==unitab[i,1])  
    } else {
      lines=which(xct_o$neighbour==unitab[i,1])
    }
    xctneibstat$nint[i]=length(lines)
    siglines=xct_o[lines,]
    if (tt=="coef"){
      xctneibstat$nsig[i]=sum(siglines$coef_p<=0.05,na.rm=T)
      xctneibstat$npos[i]=sum(siglines$coef>0 & siglines$coef_p<=0.05,na.rm=T)
      xctneibstat$nneg[i]=sum(siglines$coef<0 & siglines$coef_p<=0.05,na.rm=T)
      xctneibstat$posfrac[i]=xctneibstat$npos[i]/xctneibstat$nsig[i]  
    } else {
      xctneibstat$nsig[i]=sum(siglines$rho_pval<=0.05,na.rm=T)
      xctneibstat$npos[i]=sum(siglines$rho>0 & siglines$rho_pval<=0.05,na.rm=T)
      xctneibstat$nneg[i]=sum(siglines$rho<0 & siglines$rho_pval<=0.05,na.rm=T)
      xctneibstat$posfrac[i]=xctneibstat$npos[i]/xctneibstat$nsig[i]
    }
    
    if (length(lines)>1){
      sigcount=colSums(xctrsigmat[lines,],na.rm=T)
      poscount=colSums(xctrsignmat[lines,],na.rm=T)
      xctrneibcount[i,]=sigcount
      xctrposcount[i,]=poscount
      xctneibstat$z[i]=(xctneibstat$nsig[i]-mean(sigcount))/sd(sigcount)
    } else {
      xctrneibcount[i,]=xctrsigmat[lines,]
      xctrposcount[i,]=xctrsignmat[lines,]
      xctneibstat$z[i]=0
    }
  }
  xctneibstat1=xctneibstat[xctneibstat$nsig>1,]
  xctrneibcount[is.na(xctrneibcount)]=0
  xctrposcount[is.na(xctrposcount)]=0
  xctrneibz=matrix(data=0,nrow=nrow(unitab),ncol=ncol(xctrsigmat))
  
  for (i in 1:nrow(xctrneibz)){
    m=mean(xctrneibcount[i,])
    s=sd(xctrneibcount[i,])
    xctrneibz[i,]=(xctrneibcount[i,]-m)/s
  }
  xctrneibz[is.na(xctrneibz)]=0
  
  xctthres=1
  fdrxct=1
  while (fdrxct>maxfdr){
    xctthres=xctthres+0.05
    fdrxct=mean(colSums((xctrneibz>=xctthres)*(xctrneibcount>1),na.rm=T))/sum(xctneibstat1$z>=xctthres,na.rm=T)
  }
  xctneibstat$enriched=1*(xctneibstat$z>=xctthres & xctneibstat$nsig>1)
  
  rposfrac=xctrposcount/xctrneibcount
  posfracdist=vector()
  nsigdist=vector()
  for (i in 1:100){
    lines=which(xctrneibcount[,i]>1) #xctrneibz[,i]>=xctthres &
    posfracdist=c(posfracdist,rposfrac[lines,i])
    nsigdist=c(nsigdist,xctrneibcount[lines,i])
  } 
  posfracsample=vector()
  enrsigdist=xctneibstat$nsig[xctneibstat$enriched==1]
  for (i in 2:20){
    ni=sum(enrsigdist==i)
    rvalues=posfracdist[nsigdist==i]
    
    if (ni>0){
      samplevalues=sample(rvalues,ni*5,replace=(length(rvalues)<ni*5))
      posfracsample=c(posfracsample,samplevalues)
    } 
    
  }
  
  for (i in 3:10){
    ni=sum(enrsigdist>(i-1)*10+1 & enrsigdist<=10*i)
    rvalues=posfracdist[nsigdist>(i-1)*10+1 & nsigdist<=i*10]
    if (ni>0 & length(rvalues)>0){
      samplevalues=sample(rvalues,ni*5,replace=(length(rvalues)<ni*5))
      posfracsample=c(posfracsample,samplevalues)
    }
    
  }
  ni=sum(enrsigdist>100)
  rvalues=posfracdist[nsigdist>100]
  if (ni>0 & length(rvalues)>0){
    samplevalues=sample(rvalues,ni*5,replace=(length(rvalues)<ni*5))
    posfracsample=c(posfracsample,samplevalues)
  }
  
  
  list(stat=xctneibstat,fdrz=xctthres,rposfrac=posfracdist,rnsig=nsigdist,sample=posfracsample)
}



###################

###################
sigdist=function(xct_obs,rep,pcol){
  sigdf=xct_obs[which(xct_obs[,pcol]<=0.05),]
  rmeann=vector()
  rn=vector()
  rmeandn=vector()
  rdn=vector()
  neibdist=as.data.frame(table(sigdf$neighbour))
  drivdist=as.data.frame(table(sigdf$driver))
  
  for (i in 1:rep){
    rsig=xct_obs[sample.int(nrow(xct_obs),nrow(sigdf),replace=F),]
    rneibdist=as.data.frame(table(rsig$neighbour))
    rmeann[i]=mean(rneibdist$Freq)
    rn[i]=nrow(rneibdist)
    
    rdrivdist=as.data.frame(table(rsig$driver))
    rmeandn[i]=mean(rdrivdist$Freq)
    rdn[i]=nrow(rdrivdist)
  }
  meann=mean(neibdist$Freq)
  nr=nrow(neibdist)
  meandn=mean(drivdist$Freq)
  ndr=nrow(drivdist)
  zmeann=(meann-mean(rmeann))/sd(rmeann)
  zmeandn=(meandn-mean(rmeandn))/sd(rmeandn)
  znr=(nr-mean(rn))/sd(rn)
  zndr=(ndr-mean(rdn))/sd(rdn)
  outlist=list(zmeann=zmeann,meann=meann,rmeann=rmeann, zmeandn=zmeandn,meandn=meandn,rmeandn=rmeandn, znr=znr,nr=nr,rnr=rn, zndr=zndr,ndr=ndr,rndr=rdn )
}



###################

###################
logicOR=function(logicA,logicB){
  a=sum(logicA*logicB)
  b=sum(logicA*(!logicB))
  c=sum((!logicA)*logicB)
  d=sum((!logicA)*(!logicB))
  abcd=matrix(data=c(a,b,c,d),nrow=2,ncol=2,byrow=T)
  h=fisher.test(abcd)
  output=c(h$estimate, h$conf.int[1], h$conf.int[2], h$p.value)
}


###################

###################
mkdf2plotSpearman=function(dname,nname,ctdtab,ctetab,minmut){
  a=ctdtab[,c(1,2)]
  a$mutfreq=ctdtab[,dname]
  a$mutrfreq=a$mutfreq/a$Freq
  a$meanexp=ctetab[,nname]
  b=a[a$mutfreq>=minmut,]
  b
}


###################

###################
mkdf2plotSpearmanMB=function(dname,nname,ctdtab,ctetab,ctmbtab,minmut){
  a=ctdtab[,c(1,2)]
  a$mtb=log10(ctmbtab$mtb+1)
  a$mutfreq=ctdtab[,dname]
  a$mutrfreq=a$mutfreq/a$Freq
  a$meanexp=ctetab[,nname]
  a$mutfreqmb=a$mutrfreq/a$mtb
  b=a[a$mutfreq>=minmut,]
  b
}


###################

###################
mkdf2plotLM=function(dname,nname,ctt,mutt,expt,ctdtab,minmut){
  outt=data.frame(ctt,stringsAsFactors = F)
  outt$mutated=unlist(mutt[,dname])
  outt$exp=unlist(expt[,nname])
  a=ctdtab[,c(1,2)]
  a$mutfreq=ctdtab[,dname]
  b=a[a$mutfreq>=minmut,]
  d=data.frame(outt[is.element(outt$cancer_type,b$cancer_type),])
  d
}



