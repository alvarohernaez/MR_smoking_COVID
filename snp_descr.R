rm(list=ls())

library(chron)
library(colorspace)
library(mime)
library(dichromat)
library(munsell)
library(labeling)
library(rlang)
library(stringi)
library(evaluate)
library(highr)
library(markdown)
library(yaml)
library(backports)
library(jsonlite)
library(digest)
library(plyr)
library(reshape2)
library(scales)
library(tibble)
library(lazyeval)
library(RColorBrewer)
library(stringr)
library(knitr)
library(magrittr)
library(checkmate)
library(htmlwidgets)
library(viridisLite)
library(Rcpp)
library(Formula)
library(ggplot2)
library(latticeExtra)
library(acepack)
library(gtable)
library(gridEXtra)
library(data.table)
library(htmlTable)
library(viridis)
library(htmltools)
library(base64enc)
library(minqa)
library(RcppEigen)
library(lme4)
library(SparseM)
library(MatrixModels)
library(pbkrtest)
library(quantreg)
library(car)
library(htmlTable)
library(Hmisc)
library(survival)
library(foreign)
library(bitops)
library(caTools)
library(gplots)
library(ROCR)
library(RODBC)
library(compareGroups)
library(nlme)
library(vcd)
library(psy)
library(irr)
library(boot)
library(tibble)
library(haven)
library(icenReg)
library(arm)
library(standardize)
library(MASS)
library(sandwich)   
library(lmtest)
library(gam)
library(smoothHR)
library(meta)
library(metafor)
library(mgcv)
library(gratia)
library(MuMIn)
library(plotrix)
library(tidyr)
library(nephro)
library(miceadds)


### GUAPAS ###
##############

guapa<-function(x)
{
  redondeo<-ifelse(abs(x)<0.00001,signif(x,1),
                   ifelse(abs(x)<0.0001,signif(x,1),
                          ifelse(abs(x)<0.001,signif(x,1),
                                 ifelse(abs(x)<0.1,sprintf("%.3f",round(x,3)),
                                        ifelse(abs(x)<1,sprintf("%.2f",round(x,2)),
                                               ifelse(abs(x)<10,sprintf("%.2f",round(x,2)),
                                                      ifelse(abs(x)<100,sprintf("%.1f",round(x,1)),
                                                             ifelse(abs(x)>=100,round(x,0),round(x,0)))))))))
  return(redondeo)
}


ic_guapa<-function(x,y,z)
{
  ic<-paste(x," [",y,"; ",z,"]",sep="")
  return(ic)
}

ic_guapa2<-function(x,y,z)
{
  ic<-paste(x," (",y," to ",z,")",sep="")
  return(ic)
}

pval_guapa<-function(x)
{
  pval<-ifelse(x<0.00001,"<0.00001",
               ifelse(x<0.001,"<0.001",guapa(x)))
  return(pval)
}

pval_guapa2<-function(x)
{
  pval<-ifelse(x<0.00001," < 0.00001",
               ifelse(x<0.001," < 0.001",paste(" = ",guapa(x),sep="")))
  return(pval)
}

mean_ic_guapa <- function(x, na.rm=FALSE) 
{
  if (na.rm) x <- na.omit(x)
  se<-sqrt(var(x)/length(x))
  z<-qnorm(1-0.05/2)
  media<-mean(x)
  ic95a<-guapa(media-(z*se))
  ic95b<-guapa(media+(z*se))
  media<-guapa(media)
  ic_ok<-ic_guapa(media,ic95a,ic95b)
  return(ic_ok)
}

mean_sd_guapa <- function(x) 
{
  media<-guapa(mean(x, na.rm=TRUE))
  sd<-guapa(sd(x, na.rm=TRUE))
  end<-paste(media," (",sd,")",sep="")
  return(end)
}

beta_se_ic_guapa <- function(x, y) 
{
  z<-qnorm(1-0.05/2)
  ic95a<-guapa(x-(z*y))
  ic95b<-guapa(x+(z*y))
  media<-guapa(x)
  ic_ok<-ic_guapa(media,ic95a,ic95b)
  return(ic_ok)
}

risk_se_ic_guapa <- function(x,y) 
{
  z<-qnorm(1-0.05/2)
  hr<-guapa(exp(x))
  ic95a<-guapa(exp(x-(z*y)))
  ic95b<-guapa(exp(x+(z*y)))
  ic_ok<-ic_guapa(hr,ic95a,ic95b)
  return(ic_ok)
}

risk_se_ic_guapa2 <- function(x,y) 
{
  z<-qnorm(1-0.05/2)
  hr<-round(exp(x),5)
  ic95a<-round(exp(x-(z*y)),5)
  ic95b<-round(exp(x+(z*y)),5)
  ic_ok<-ic_guapa(hr,ic95a,ic95b)
  return(ic_ok)
}

header.true <- function(df)
{
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}

z<-qnorm(1-0.05/2)

setwd("C:/Users/ALHE/OneDrive - Folkehelseinstituttet/Articles/MR_smoking_COVID/Results")


### smkinit ###

smkinit<-read.csv2("./one_sample_MR/gwas_smkinit_liu2019.csv",header=TRUE,sep=";",dec=".")
smkinit_gwas<-as.data.frame(smkinit$rsid)
names(smkinit_gwas)<-c("rsid")
smkinit_1smr<-read.csv2("./one_sample_MR/snps_smkinit_1smr.csv",header=FALSE,sep=";",dec=".")
names(smkinit_1smr)<-c("rsid")
smkinit_1smr$onesmr<-1
smkinit_2smr<-read.csv2("./two_sample_MR/snps_smkinit_2smr.csv",header=FALSE,sep=";",dec=".")
names(smkinit_2smr)<-c("rsid")
smkinit_2smr$twosmr<-1
smkinit_gwas<-merge2(smkinit_gwas,smkinit_1smr,by.id=c("rsid"),all.x=TRUE,sort=FALSE)
smkinit_gwas<-merge2(smkinit_gwas,smkinit_2smr,by.id=c("rsid"),all.x=TRUE,sort=FALSE)
smkinit_gwas$onesmr<-with(smkinit_gwas,ifelse(is.na(onesmr),"No","Yes"))
smkinit_gwas$twosmr<-with(smkinit_gwas,ifelse(is.na(twosmr),"No","Yes"))
smkinit<-smkinit[,c("rsid","chr","pos","other_allele","effect_allele","eaf","beta","se","pval")]
smkinit<-merge2(smkinit,smkinit_gwas,by.id=c("rsid"),all.x=TRUE,sort=FALSE)
write.table(smkinit,file="./snps_smkinit.csv",sep=";",col.names=TRUE,row.names=FALSE)


### cigday ###

cigday<-read.csv2("./one_sample_MR/gwas_cigday_liu2019.csv",header=TRUE,sep=";",dec=".")
cigday_gwas<-as.data.frame(cigday$rsid)
names(cigday_gwas)<-c("rsid")
cigday_1smr<-read.csv2("./one_sample_MR/snps_cigday_1smr.csv",header=FALSE,sep=";",dec=".")
names(cigday_1smr)<-c("rsid")
cigday_1smr$onesmr<-1
cigday_2smr<-read.csv2("./two_sample_MR/snps_cigday_2smr.csv",header=FALSE,sep=";",dec=".")
names(cigday_2smr)<-c("rsid")
cigday_2smr$twosmr<-1
cigday_gwas<-merge2(cigday_gwas,cigday_1smr,by.id=c("rsid"),all.x=TRUE,sort=FALSE)
cigday_gwas<-merge2(cigday_gwas,cigday_2smr,by.id=c("rsid"),all.x=TRUE,sort=FALSE)
cigday_gwas$onesmr<-with(cigday_gwas,ifelse(is.na(onesmr),"No","Yes"))
cigday_gwas$twosmr<-with(cigday_gwas,ifelse(is.na(twosmr),"No","Yes"))
cigday<-cigday[,c("rsid","chr","pos","other_allele","effect_allele","eaf","beta","se","pval")]
cigday<-merge2(cigday,cigday_gwas,by.id=c("rsid"),all.x=TRUE,sort=FALSE)
write.table(cigday,file="./snps_cigday.csv",sep=";",col.names=TRUE,row.names=FALSE)


### agesmk ###

agesmk<-read.csv2("./one_sample_MR/gwas_agesmk_liu2019.csv",header=TRUE,sep=";",dec=".")
agesmk_gwas<-as.data.frame(agesmk$rsid)
names(agesmk_gwas)<-c("rsid")
agesmk_1smr<-read.csv2("./one_sample_MR/snps_agesmk_1smr.csv",header=FALSE,sep=";",dec=".")
names(agesmk_1smr)<-c("rsid")
agesmk_1smr$onesmr<-1
agesmk_2smr<-read.csv2("./two_sample_MR/snps_agesmk_2smr.csv",header=FALSE,sep=";",dec=".")
names(agesmk_2smr)<-c("rsid")
agesmk_2smr$twosmr<-1
agesmk_gwas<-merge2(agesmk_gwas,agesmk_1smr,by.id=c("rsid"),all.x=TRUE,sort=FALSE)
agesmk_gwas<-merge2(agesmk_gwas,agesmk_2smr,by.id=c("rsid"),all.x=TRUE,sort=FALSE)
agesmk_gwas$onesmr<-with(agesmk_gwas,ifelse(is.na(onesmr),"No","Yes"))
agesmk_gwas$twosmr<-with(agesmk_gwas,ifelse(is.na(twosmr),"No","Yes"))
agesmk<-agesmk[,c("rsid","chr","pos","other_allele","effect_allele","eaf","beta","se","pval")]
agesmk<-merge2(agesmk,agesmk_gwas,by.id=c("rsid"),all.x=TRUE,sort=FALSE)
write.table(agesmk,file="./snps_agesmk.csv",sep=";",col.names=TRUE,row.names=FALSE)


### smkces ###

smkces<-read.csv2("./one_sample_MR/gwas_smkces_liu2019.csv",header=TRUE,sep=";",dec=".")
smkces_gwas<-as.data.frame(smkces$rsid)
names(smkces_gwas)<-c("rsid")
smkces_1smr<-read.csv2("./one_sample_MR/snps_smkces_1smr.csv",header=FALSE,sep=";",dec=".")
names(smkces_1smr)<-c("rsid")
smkces_1smr$onesmr<-1
smkces_2smr<-read.csv2("./two_sample_MR/snps_smkces_2smr.csv",header=FALSE,sep=";",dec=".")
names(smkces_2smr)<-c("rsid")
smkces_2smr$twosmr<-1
smkces_gwas<-merge2(smkces_gwas,smkces_1smr,by.id=c("rsid"),all.x=TRUE,sort=FALSE)
smkces_gwas<-merge2(smkces_gwas,smkces_2smr,by.id=c("rsid"),all.x=TRUE,sort=FALSE)
smkces_gwas$onesmr<-with(smkces_gwas,ifelse(is.na(onesmr),"No","Yes"))
smkces_gwas$twosmr<-with(smkces_gwas,ifelse(is.na(twosmr),"No","Yes"))
smkces<-smkces[,c("rsid","chr","pos","other_allele","effect_allele","eaf","beta","se","pval")]
smkces<-merge2(smkces,smkces_gwas,by.id=c("rsid"),all.x=TRUE,sort=FALSE)
write.table(smkces,file="./snps_smkces.csv",sep=";",col.names=TRUE,row.names=FALSE)



### FOREST PLOTS ###
####################

# Suppress scientific notation for coefficients (no "4e-04", instead "0.0004")
options(scipen=999)

tab1<-NULL
tab2<-NULL
tab3<-NULL
dat01<-read.csv2("./one_sample_MR/main.csv",header=TRUE,sep=";",dec=".")[1:3,]
dat02<-read.csv2("./one_sample_MR/mvmr.csv",header=TRUE,sep=";",dec=".")[1:3,]
dat03<-read.csv2("./two_sample_MR/results_2smr.csv",header=TRUE,sep=";",dec=".")[1:3,]

tab1<-rbind(tab1,cbind(dat01$Raw.OR[1],dat01$Raw.lo[1],dat01$Raw.hi[1]),
            cbind(dat01$MV.OR[1],dat01$MV.lo[1],dat01$MV.hi[1]),
            cbind(dat01$MR.OR[1],dat01$MR.lo[1],dat01$MR.hi[1]),
            cbind(dat02$MVMR.OR[1],dat02$MVMR.lo[1],dat02$MVMR.hi[1]),
            cbind(dat03$ivw_or[1],dat03$ivw_or_lo[1],dat03$ivw_or_hi[1]))
tab1<-cbind(tab1,c("05. Raw logistic regression","04. Multivariable regression","03. One-sample MR","02. One-sample multivariable MR","01. Two-sample MR"))
tab1<-as.data.frame(tab1)
tab1$group<-c(" All participants")
tab2<-rbind(tab2,cbind(dat01$Raw.OR[2],dat01$Raw.lo[2],dat01$Raw.hi[2]),
            cbind(dat01$MV.OR[2],dat01$MV.lo[2],dat01$MV.hi[2]),
            cbind(dat01$MR.OR[2],dat01$MR.lo[2],dat01$MR.hi[2]),
            cbind(dat02$MVMR.OR[2],dat02$MVMR.lo[2],dat02$MVMR.hi[2]),
            cbind(dat03$ivw_or[1],dat03$ivw_or_lo[1],dat03$ivw_or_hi[1]))
tab2<-cbind(tab2,c("05. Raw logistic regression","04. Multivariable regression","03. One-sample MR","02. One-sample multivariable MR","01. Two-sample MR"))
tab2<-as.data.frame(tab2)
tab2$group<-c(" Women")
tab3<-rbind(tab3,cbind(dat01$Raw.OR[3],dat01$Raw.lo[3],dat01$Raw.hi[3]),
            cbind(dat01$MV.OR[3],dat01$MV.lo[3],dat01$MV.hi[3]),
            cbind(dat01$MR.OR[3],dat01$MR.lo[3],dat01$MR.hi[3]),
            cbind(dat02$MVMR.OR[3],dat02$MVMR.lo[3],dat02$MVMR.hi[3]),
            cbind(dat03$ivw_or[1],dat03$ivw_or_lo[1],dat03$ivw_or_hi[1]))
tab3<-cbind(tab3,c("05. Raw logistic regression","04. Multivariable regression","03. One-sample MR","02. One-sample multivariable MR","01. Two-sample MR"))
tab3<-as.data.frame(tab3)
tab3$group<-c("Men")

tab<-rbind(tab1,tab2,tab3)
names(tab)<-c("or","lci","uci","class","group")
tab$coef<-ic_guapa2(guapa(as.numeric(tab$or)),guapa(as.numeric(tab$lci)),guapa(as.numeric(tab$uci)))
tab<-tab[c(1,2,3,4,5,6,7,8,9,11,12,13,14),]
write.table(tab,file="./f_smkinit.csv",sep=";",col.names=TRUE,row.names=FALSE)

tab1<-NULL
tab2<-NULL
tab3<-NULL
dat01<-read.csv2("./one_sample_MR/main.csv",header=TRUE,sep=";",dec=".")[4:6,]
dat02<-read.csv2("./one_sample_MR/mvmr.csv",header=TRUE,sep=";",dec=".")[4:6,]
dat03<-read.csv2("./two_sample_MR/results_2smr.csv",header=TRUE,sep=";",dec=".")[4:6,]

tab1<-rbind(tab1,cbind(dat01$Raw.OR[1],dat01$Raw.lo[1],dat01$Raw.hi[1]),
            cbind(dat01$MV.OR[1],dat01$MV.lo[1],dat01$MV.hi[1]),
            cbind(dat01$MR.OR[1],dat01$MR.lo[1],dat01$MR.hi[1]),
            cbind(dat02$MVMR.OR[1],dat02$MVMR.lo[1],dat02$MVMR.hi[1]),
            cbind(dat03$ivw_or[1],dat03$ivw_or_lo[1],dat03$ivw_or_hi[1]))
tab1<-cbind(tab1,c("05. Raw logistic regression","04. Multivariable regression","03. One-sample MR","02. One-sample multivariable MR","01. Two-sample MR"))
tab1<-as.data.frame(tab1)
tab1$group<-c(" All participants")
tab2<-rbind(tab2,cbind(dat01$Raw.OR[2],dat01$Raw.lo[2],dat01$Raw.hi[2]),
            cbind(dat01$MV.OR[2],dat01$MV.lo[2],dat01$MV.hi[2]),
            cbind(dat01$MR.OR[2],dat01$MR.lo[2],dat01$MR.hi[2]),
            cbind(dat02$MVMR.OR[2],dat02$MVMR.lo[2],dat02$MVMR.hi[2]),
            cbind(dat03$ivw_or[1],dat03$ivw_or_lo[1],dat03$ivw_or_hi[1]))
tab2<-cbind(tab2,c("05. Raw logistic regression","04. Multivariable regression","03. One-sample MR","02. One-sample multivariable MR","01. Two-sample MR"))
tab2<-as.data.frame(tab2)
tab2$group<-c(" Women")
tab3<-rbind(tab3,cbind(dat01$Raw.OR[3],dat01$Raw.lo[3],dat01$Raw.hi[3]),
            cbind(dat01$MV.OR[3],dat01$MV.lo[3],dat01$MV.hi[3]),
            cbind(dat01$MR.OR[3],dat01$MR.lo[3],dat01$MR.hi[3]),
            cbind(dat02$MVMR.OR[3],dat02$MVMR.lo[3],dat02$MVMR.hi[3]),
            cbind(dat03$ivw_or[1],dat03$ivw_or_lo[1],dat03$ivw_or_hi[1]))
tab3<-cbind(tab3,c("05. Raw logistic regression","04. Multivariable regression","03. One-sample MR","02. One-sample multivariable MR","01. Two-sample MR"))
tab3<-as.data.frame(tab3)
tab3$group<-c("Men")

tab<-rbind(tab1,tab2,tab3)
names(tab)<-c("or","lci","uci","class","group")
tab$coef<-ic_guapa2(guapa(as.numeric(tab$or)),guapa(as.numeric(tab$lci)),guapa(as.numeric(tab$uci)))
tab<-tab[c(1,2,3,4,5,6,7,8,9,11,12,13,14),]
write.table(tab,file="./f_cigday.csv",sep=";",col.names=TRUE,row.names=FALSE)

tab1<-NULL
tab2<-NULL
tab3<-NULL
dat01<-read.csv2("./one_sample_MR/main.csv",header=TRUE,sep=";",dec=".")[7:9,]
dat02<-read.csv2("./one_sample_MR/mvmr.csv",header=TRUE,sep=";",dec=".")[7:9,]
dat03<-read.csv2("./two_sample_MR/results_2smr.csv",header=TRUE,sep=";",dec=".")[10:12,]

tab1<-rbind(tab1,cbind(dat01$Raw.OR[1],dat01$Raw.lo[1],dat01$Raw.hi[1]),
            cbind(dat01$MV.OR[1],dat01$MV.lo[1],dat01$MV.hi[1]),
            cbind(dat01$MR.OR[1],dat01$MR.lo[1],dat01$MR.hi[1]),
            cbind(dat02$MVMR.OR[1],dat02$MVMR.lo[1],dat02$MVMR.hi[1]),
            cbind(dat03$ivw_or[1],dat03$ivw_or_lo[1],dat03$ivw_or_hi[1]))
tab1<-cbind(tab1,c("05. Raw logistic regression","04. Multivariable regression","03. One-sample MR","02. One-sample multivariable MR","01. Two-sample MR"))
tab1<-as.data.frame(tab1)
tab1$group<-c(" All participants")
tab2<-rbind(tab2,cbind(dat01$Raw.OR[2],dat01$Raw.lo[2],dat01$Raw.hi[2]),
            cbind(dat01$MV.OR[2],dat01$MV.lo[2],dat01$MV.hi[2]),
            cbind(dat01$MR.OR[2],dat01$MR.lo[2],dat01$MR.hi[2]),
            cbind(dat02$MVMR.OR[2],dat02$MVMR.lo[2],dat02$MVMR.hi[2]),
            cbind(dat03$ivw_or[1],dat03$ivw_or_lo[1],dat03$ivw_or_hi[1]))
tab2<-cbind(tab2,c("05. Raw logistic regression","04. Multivariable regression","03. One-sample MR","02. One-sample multivariable MR","01. Two-sample MR"))
tab2<-as.data.frame(tab2)
tab2$group<-c(" Women")
tab3<-rbind(tab3,cbind(dat01$Raw.OR[3],dat01$Raw.lo[3],dat01$Raw.hi[3]),
            cbind(dat01$MV.OR[3],dat01$MV.lo[3],dat01$MV.hi[3]),
            cbind(dat01$MR.OR[3],dat01$MR.lo[3],dat01$MR.hi[3]),
            cbind(dat02$MVMR.OR[3],dat02$MVMR.lo[3],dat02$MVMR.hi[3]),
            cbind(dat03$ivw_or[1],dat03$ivw_or_lo[1],dat03$ivw_or_hi[1]))
tab3<-cbind(tab3,c("05. Raw logistic regression","04. Multivariable regression","03. One-sample MR","02. One-sample multivariable MR","01. Two-sample MR"))
tab3<-as.data.frame(tab3)
tab3$group<-c("Men")

tab<-rbind(tab1,tab2,tab3)
names(tab)<-c("or","lci","uci","class","group")
tab$coef<-ic_guapa2(guapa(as.numeric(tab$or)),guapa(as.numeric(tab$lci)),guapa(as.numeric(tab$uci)))
tab<-tab[c(1,2,3,4,5,6,7,8,9,11,12,13,14),]
write.table(tab,file="./f_smkces.csv",sep=";",col.names=TRUE,row.names=FALSE)

tab1<-NULL
tab2<-NULL
tab3<-NULL
dat01<-read.csv2("./one_sample_MR/main.csv",header=TRUE,sep=";",dec=".")[10,]
dat02<-read.csv2("./one_sample_MR/mvmr.csv",header=TRUE,sep=";",dec=".")[10,]
dat03<-read.csv2("./two_sample_MR/results_2smr.csv",header=TRUE,sep=";",dec=".")[7:9,]

tab1<-rbind(tab1,cbind(dat01$Raw.OR[1],dat01$Raw.lo[1],dat01$Raw.hi[1]),
            cbind(dat01$MV.OR[1],dat01$MV.lo[1],dat01$MV.hi[1]),
            cbind(dat01$MR.OR[1],dat01$MR.lo[1],dat01$MR.hi[1]),
            cbind(dat02$MVMR.OR[1],dat02$MVMR.lo[1],dat02$MVMR.hi[1]),
            cbind(dat03$ivw_or[1],dat03$ivw_or_lo[1],dat03$ivw_or_hi[1]))
tab1<-cbind(tab1,c("05. Raw logistic regression","04. Multivariable regression","03. One-sample MR","02. One-sample multivariable MR","01. Two-sample MR* (women + men)"))
tab1<-as.data.frame(tab1)
tab1$group<-c("aaa")
tab2<-tab1
tab2$group<-c("bbb")
tab3<-tab1
tab3$group<-c("Women")

tab<-rbind(tab1,tab2,tab3)
names(tab)<-c("or","lci","uci","class","group")
tab$coef<-ic_guapa2(guapa(as.numeric(tab$or)),guapa(as.numeric(tab$lci)),guapa(as.numeric(tab$uci)))
tab<-tab[c(1,2,3,4,6,7,8,9,11,12,13,14,15),]
write.table(tab,file="./f_agesmk.csv",sep=";",col.names=TRUE,row.names=FALSE)


vars01<-c("smkinit","cigday","smkces","agesmk")
vars02<-c(1.35,2.3,2,2.5)

for(i in 1:length(vars01))
  
{
  filename<-paste("./f_",vars01[i],".csv",sep="")
  dat<-read.csv2(filename,header=TRUE,sep=";",dec=".")

  figure<-ggplot(data=dat,
                 aes(x=class, y=or, ymin=lci, ymax=uci)) +
    geom_hline(aes(fill=class), yintercept=1, linetype=2) +
    geom_pointrange(aes(col=class), size=0.7, shape=15) +
    geom_text(data=dat, size=4, aes(y=max(uci)*vars02[i], x=class, label=coef, hjust='inward')) +
    xlab(" ") + ylab("Risk of COVID-19 infection (odds ratio, 95% confidence interval)") +
    scale_y_continuous(trans = log2_trans()) +
    geom_errorbar(aes(ymin=lci, ymax=uci, col=class), width=0.5, cex=1) + 
    facet_wrap(~group, strip.position="left", nrow=9,scales="free_y") +
    coord_flip() +
    theme_minimal() +
    theme(legend.position="none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.background = element_rect(fill = "white", color = NA),
          axis.line = element_line(colour = "black"),
          axis.text.y=element_text(size=11),
          axis.text.x=element_text(size=11),
          axis.ticks.y=element_line(),
          axis.ticks.x=element_line(),
          axis.title.x=element_text(size=11),
          axis.title.y=element_text(size=11))
  
  namefile<-paste("./",vars01[i],".tiff",sep="")
  ggsave(filename=namefile, units="px", width=9000, height=6000, dpi=1200, bg="transparent")
  figure
}




