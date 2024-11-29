
mkctdtab=function(ctt,mutt){ #auxiliar function to compute fraction of mutated individuals per cancer type
  a=data.frame(table(ctt$cancer_type))
  names(a)=c("cancer_type","Freq")
  b=matrix(nrow=nrow(a),ncol=nrow(mutt))
  for (i in 1:nrow(a)){
    lines=which(ctt$cancer_type==a$cancer_type[i])
    if (length(lines)>1){
      b[i,]=rowSums(mutt[,ctt$patient[lines]])  
    } else {
      if (length(lines)==1){
        b[i,]=mutt[,ctt$patient[lines]]
      } else {
        b[i,]=rep(0,nrow(mutt))
      }
    }
    
  }
  b=as.data.frame(b)
  names(b)=row.names(mutt)
  a=cbind(a,b)
  a
}

mkctetab=function(ctt,expt){ #auxiliary function to compute average expression per cancer type
  a=data.frame(table(ctt$cancer_type))
  names(a)=c("cancer_type","Freq")
  b=matrix(nrow=nrow(a),ncol=nrow(expt))
  for (i in 1:nrow(a)){
    lines=which(ctt$cancer_type==a$cancer_type[i])
    if (length(lines)>1){
      b[i,]=rowMeans(expt[,ctt$patient[lines]])  
    } else {
      if (length(lines)==1){
        b[i,]=expt[,ctt$patient[lines]]
      } else {
        b[i,]=rep(0,nrow(mutt))
      }
    }
    
  }
  b=as.data.frame(b)
  names(b)=row.names(expt)
  a=cbind(a,b)
  a
}

mkrandexp=function(exptab,cct){ #auxiliary function to generate random permutation of expression table
  #rexptab=exptab[sample(nrow(exptab),nrow(exptab),replace=F),sample(ncol(exptab),ncol(exptab),replace=F)]
  rexptab=exptab[,sample(ncol(exptab),ncol(exptab),replace=F)]
  names(rexptab)=names(exptab)
  #row.names(rexptab)=row.names(exptab)
  rctetab=mkctetab(cct,rexptab)
  list(rexptab=rexptab,rctetab=rctetab)
}


mkpairedtest=function(dname,nname,mutt,expt,ctt,ctdtab,minmut){ #function to test if neighbour expression is influenced by driver mutation through the comparison of matched tumour and normal samples
  outt=data.frame(ctt,stringsAsFactors = F)
  outt$mutated=as.numeric(mutt[dname,])
  outt$exp=as.numeric(expt[nname,])
  a=ctdtab[,c(1,2)]
  a$mutfreq=ctdtab[,dname]
  a$mutrfreq=a$mutfreq/a$Freq
  b=a[a$mutfreq>=minmut,]
  if (nrow(b)>1){ #at least 2 cancer types with mutated individuals
    d=data.frame(outt[is.element(outt$cancer_type,b$cancer_type),])
    fitglm=glm(as.numeric(exp)~mutated+cancer_type,data=d)
    sumtab=summary(fitglm)
    output=c(paircoef=sumtab$coefficients["mutated",1],paircoef_p=sumtab$coefficients["mutated",4])  
  } else {
    output=c(paircoef=0,paircoef_p=1)
  }
  
}

mktestslmd=function(dname,nnames,mutt,expt,ctt,ctdtab,ctetab,minmut){ #function to test if neighbour expression is associated with driver mutation status both within and between cancer types
  outt=data.frame(ctt,stringsAsFactors = F)
  outt$mutated=as.numeric(mutt[dname,])
  neibexp=(expt[nnames,])
  outt=cbind(outt,t(neibexp))
  a=ctdtab[,c(1,2)]
  a$mutfreq=ctdtab[,dname]
  a=cbind(a,ctetab[,nnames])
  a$mutrfreq=a$mutfreq/a$Freq
  b=a[a$mutfreq>=minmut,]
  output=data.frame(rho=rep(0,length(nnames)),rho_p=1,coef=0,coef_p=1)
  if (nrow(b)>2){ #at least 3 cancer types with mutated individuals
    d=data.frame(outt[is.element(outt$cancer_type,b$cancer_type),])
    for (i in 1:length(nnames)){
      c=cor.test(b[,i+3],b$mutrfreq,method="spearman",exact=F)
      d$exp=d[,i+6]
      fitglm=lm(as.numeric(exp)~mutated+cancer_type,data=d)
      sumtab=summary(fitglm)  
      output[i,]=c(c$estimate,c$p.value,sumtab$coefficients["mutated",1],sumtab$coefficients["mutated",4])
    }
  } 
  output$driver=dname
  output$neighbour=nnames
  output
}

mkdf2plotSpearman=function(dname,nname,ctdtab,ctetab,minmut){
  a=ctdtab[,c(1,2)]
  a$mutfreq=ctdtab[,dname]
  a$mutrfreq=a$mutfreq/a$Freq
  a$meanexp=ctetab[,nname]
  b=a[a$mutfreq>=minmut,]
  b
}

mkdf2plotLM=function(dname,nname,ctt,mutt,expt,ctdtab,minmut){
  outt=data.frame(ctt,stringsAsFactors = F)
  outt$mutated=as.numeric(mutt[dname,])
  outt$exp=as.numeric(expt[nname,])
  a=ctdtab[,c(1,2)]
  a$mutfreq=ctdtab[,dname]
  b=a[a$mutfreq>=minmut,]
  d=data.frame(outt[is.element(outt$cancer_type,b$cancer_type),])
  d
}


survtest=function(neib,ct,clindf,expdf){
  outdf=clindf
  outdf$neibexp=as.numeric(expdf[neib,])
  outdf=outdf[outdf$cancer.type.abbreviation==ct,]
  outdf$highexp=1*(outdf$neibexp>=median(outdf$neibexp))
  h=coxph(Surv(OS.time,OS)~highexp,data=outdf)
  hh=summary(h)
  list(survcoef=h$coefficients[[1]],survp=hh$coefficients[[5]],df=outdf)
}

