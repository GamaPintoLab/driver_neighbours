mutdata=read.delim("./data/raw/mc3.v0.2.8.PUBLIC.xena",header=T,stringsAsFactors = F)

library(arrow)
expdf=read_feather(file="./data/processed/expression.feather")
maingraph=read.csv("./data/processed/main_graph.csv",header=T,stringsAsFactors = F)
mtb=read.csv("./data/processed/mutation_burden.csv",header=T,stringsAsFactors = F)

library(coselens)

mutdata=mutdata[,c(1,2,3,5,6)]
names(mutdata)=c("sampleID", "chr", "pos", "ref", "alt")

table(mtb$cancer_type)

coselens_prep=function(dname,nname,ctype,ctt,mutall,expmat){
  ctpatients=ctt$patient[ctt$cancer_type==ctype]
  ctsamples=paste(ctpatients,"01",sep="-")
  mut_ct=mutall[is.element(mutall$sampleID,ctsamples),]
  exp_ct=expdf[ctt$cancer_type==ctype,nname]
  med_exp=median(exp_ct)
  upsamples=ctsamples[exp_ct>=med_exp]
  mut_ct_up=mut_ct[is.element(mut_ct$sampleID,upsamples),]
  mut_ct_down=mut_ct[!is.element(mut_ct$sampleID,upsamples),]
  coselens_res = coselens(mut_ct_up, mut_ct_down, subset.genes.by = dname)
  coselens_summary=coselens_res$substitutions
  coselens_summary=rbind(coselens_summary,coselens_res$indels)
  coselens_summary=rbind(coselens_summary,coselens_res$missense_sub)
  coselens_summary=rbind(coselens_summary,coselens_res$truncating_sub)
  coselens_summary
}

#1
cos_cic_hdac1_lgg=coselens_prep("CIC","HDAC1","LGG",mtb,mutdata,expdf)
#2
cos_cic_lsm14a_brca=coselens_prep("CIC","LSM14A","BRCA",mtb,mutdata,expdf)
#3
cos_fubp1_trim67_lgg=coselens_prep("FUBP1","TRIM67","LGG",mtb,mutdata,expdf)
#4
cos_gtf2i_mab21l2_thym=coselens_prep("GTF2I","MAB21L2","THYM",mtb,mutdata,expdf)
#5
cos_nfe2l2_ugt1a6_cesc=coselens_prep("NFE2L2","UGT1A6","CESC",mtb,mutdata,expdf)
#6
cos_idh1_serpinh1_lgg=coselens_prep("IDH1","SERPINH1","LGG",mtb,mutdata,expdf)
#7
cos_egfr_sec61g_gbm=coselens_prep("EGFR","SEC61G","GBM",mtb,mutdata,expdf)
#8
cos_stk11_sik1_luad=coselens_prep("STK11","SIK1","LUAD",mtb,mutdata,expdf)
#9
cos_ctnnb1_ezr_blca=coselens_prep("CTNNB1","EZR","BLCA",mtb,mutdata,expdf)
#10
cos_atrx_cct8l2_lgg=coselens_prep("ATRX","CCT8L2","LGG",mtb,mutdata,expdf)


cos_10=rbind(cos_cic_hdac1_lgg,cos_tp53_abcf1_brca)
cos_10=rbind(cos_10,cos_fubp1_trim67_lgg)
cos_10=rbind(cos_10,cos_gtf2i_mab21l2_thym)
cos_10=rbind(cos_10,cos_nfe2l2_ugt1a6_cesc)
cos_10=rbind(cos_10,cos_idh1_serpinh1_lgg)
cos_10=rbind(cos_10,cos_egfr_sec61g_gbm)
cos_10=rbind(cos_10,cos_stk11_sik1_luad)
cos_10=rbind(cos_10,cos_ctnnb1_ezr_blca)
cos_10=rbind(cos_10,cos_atrx_cct8l2_lgg)
write.csv(cos_10,file="./data/processed/cos10.csv",row.names = F)
