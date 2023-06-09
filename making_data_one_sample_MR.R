rm(list=ls())

library(chron)
library(colorspace)
library(mime)
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
library(mice)
library(officer)
library(uuid)
library(writexl)
library(HardyWeinberg)
library(compareGroups)
library(nlme)
library(vcd)
library(boot)
library(tibble)
library(haven)
library(MASS)
library(sandwich)   
library(lmtest)
library(gam)
library(smoothHR)
library(metafor)
library(DBI)
library(mitools)
library(RcppArmadillo)
library(miceadds)
library(dplyr)
library(estimatr)
library(lubridate)
library(snakecase)
library(janitor)
library(fmsb)


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
               ifelse(x<0.001,"<0.001",
                      ifelse(abs(x)<0.01,sprintf("%.3f",round(x,3)),
                             ifelse(abs(x)<0.1,sprintf("%.3f",round(x,3)),
                                    ifelse(abs(x)<1,sprintf("%.3f",round(x,3)),guapa(x))))))
  return(pval)
}

pval_guapa2<-function(x)
{
  pval<-ifelse(x<0.00001," < 0.00001",
               ifelse(x<0.001," < 0.001",
                      ifelse(abs(x)<0.01,sprintf("%.3f",round(x,3)),
                             ifelse(abs(x)<0.1,sprintf("%.3f",round(x,3)),
                                    ifelse(abs(x)<1,sprintf("%.3f",round(x,3)),guapa(x))))))
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

beta_se_ic_guapa2 <- function(x, y) 
{
  z<-qnorm(1-0.05/2)
  ic95a<-guapa(x-(z*y))
  ic95b<-guapa(x+(z*y))
  media<-guapa(x)
  ic_ok<-ic_guapa2(media,ic95a,ic95b)
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
  hr<-guapa(exp(x))
  ic95a<-guapa(exp(x-(z*y)))
  ic95b<-guapa(exp(x+(z*y)))
  ic_ok<-ic_guapa2(hr,ic95a,ic95b)
  return(ic_ok)
}

risk_se_ic_guapa3 <- function(x,y) 
{
  z<-qnorm(1-0.05/2)
  hr<-round(exp(x),3)
  ic95a<-round(exp(x-(z*y)),3)
  ic95b<-round(exp(x+(z*y)),3)
  ic_ok<-ic_guapa2(hr,ic95a,ic95b)
  return(ic_ok)
}

header.true <- function(df)
{
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}

z<-qnorm(1-0.05/2)

dir.create("N:/durable/Projects/ALHE_smoking_covid")
dir.create("N:/durable/Projects/ALHE_smoking_covid/Data")
dir.create("N:/durable/Projects/ALHE_smoking_covid/Outputs")
dir.create("N:/durable/Projects/ALHE_smoking_covid/Outputs/Descriptive")
dir.create("N:/durable/Projects/ALHE_smoking_covid/Outputs/Results")

setwd("N:/durable/Projects/ALHE_smoking_covid")


###########################
### CALCULATION OF GRSs ###
###########################

### SmkInit ###
###############

dat<-as.data.frame(fread("./Data/smkinit_liu2019.raw",header=TRUE,sep="\t",sep2="\t"))
dat<-dat[,c(2,7:dim(dat)[2])]

gwas<-read.csv2("./Data/gwas_smkinit_liu2019.csv",header=TRUE,sep=";",dec=".")
names(gwas)<-tolower(names(gwas))
gwas<-rename.vars(gwas,
                  from=c("effect_allele_freq"),
                  to=c("eaf"))
gwas$effect_allele<-tolower(gwas$effect_allele)


# 351 SNPs IN MoBa AND GWAS #

bbb<-as.character(sort(gwas$rsid,decreasing=FALSE))
dat_inv<-header.true(data.frame(t(dat)))
ccc<-as.character(sort(sub("\\_.*", "", rownames(dat_inv)),decreasing=FALSE))
write.table(ccc,file="./Outputs/Descriptive/snps_smkinit_1smr.csv",sep=";", col.names=FALSE ,row.names=FALSE)
length(which(ccc%in%bbb)) ### VARIABLES IN dat_inv AND gwas MATCH


# GRS CALCULATION #

rsid<-as.character(sub("\\_.*", "", rownames(dat_inv)))
tested_allele_moba<-sub('.*_', '', rownames(dat_inv))
moba_info<-as.data.frame(cbind(rsid,tested_allele_moba))
moba_info$tested_allele_moba<-tolower(moba_info$tested_allele_moba)
gwas<-merge2(gwas,moba_info,by.id=c("rsid"),all.x=TRUE,sort=FALSE)
gwas<-na.omit(gwas)

gwas$same_coding<-with(gwas,ifelse(tested_allele_moba==effect_allele,1,
                                   ifelse(tested_allele_moba!=effect_allele,0,NA)))

gwas$flip<-with(gwas,ifelse(beta>=0 & same_coding==1,0,
                            ifelse(beta<0 & same_coding==0,0,
                                   ifelse(beta>=0 & same_coding==0,1,
                                          ifelse(beta<0 & same_coding==1,1,NA)))))

gwas$maf<-with(gwas,ifelse(same_coding==1,eaf,
                           ifelse(same_coding==0,1-eaf,NA)))
gwas$coef<-with(gwas,ifelse(beta>0,beta,
                            ifelse(beta<0,-beta,NA)))

n_snps<-dim(dat)[2]-1
colnames(dat)<-sub("\\_.*", "", colnames(dat))
dat<-rename.vars(dat,from=c("IID"),to=c("id"))

vars01<-colnames(dat)[2:dim(dat)[2]]
vars02<-paste(vars01,"_flip",sep="")
vars03<-paste(vars01,"_coef",sep="")
vars04<-paste(vars01,"_weight",sep="")

for(i in 1:length(vars01))
  
{
  dat[,vars02[i]]<-gwas[gwas$rsid==vars01[i],"flip"]
  dat[,vars03[i]]<-gwas[gwas$rsid==vars01[i],"coef"]
  dat[,vars04[i]]<-with(dat,ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==0,0*dat[,vars03[i]],
                                   ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==1,1*dat[,vars03[i]],
                                          ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==2,2*dat[,vars03[i]],
                                                 ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==0,2*dat[,vars03[i]],
                                                        ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==1,1*dat[,vars03[i]],
                                                               ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==2,0*dat[,vars03[i]],NA)))))))
}

dat_coef<-dat[,vars03]
dat$coef_sum<-(dat_coef %>%
                 replace(is.na(.), 0) %>%
                 mutate(coef_sum = rowSums(.)))$coef_sum[1]

vars05<-append("id",vars04,after=1)
dat_weight<-dat[,vars05]
dat_weight$weight_sum<-(dat_weight[,2:dim(dat_weight)[2]] %>%
                          replace(is.na(.), 0) %>%
                          mutate(weight_sum = rowSums(.)))$weight_sum
dat_weight<-dat_weight[,c("id","weight_sum")]
dat<-merge(dat,dat_weight,by="id")
dat$smkinit_grs<-(dat$weight_sum/dat$coef_sum)*n_snps
dat<-dat[,c("id","smkinit_grs")]
save(dat,file="./Data/smkinit_grs.RData")


### CigDay ###
##############

dat<-as.data.frame(fread("./Data/cigday_liu2019.raw",header=TRUE,sep="\t",sep2="\t"))
dat<-dat[,c(2,7:dim(dat)[2])]

gwas<-read.csv2("./Data/gwas_cigday_liu2019.csv",header=TRUE,sep=";",dec=".")
names(gwas)<-tolower(names(gwas))
gwas<-rename.vars(gwas,
                  from=c("effect_allele_freq"),
                  to=c("eaf"))
gwas$effect_allele<-tolower(gwas$effect_allele)


# 53 SNPs IN MoBa AND GWAS #

bbb<-as.character(sort(gwas$rsid,decreasing=FALSE))
dat_inv<-header.true(data.frame(t(dat)))
ccc<-as.character(sort(sub("\\_.*", "", rownames(dat_inv)),decreasing=FALSE))
write.table(ccc,file="./Outputs/Descriptive/snps_cigday_1smr.csv",sep=";", col.names=FALSE ,row.names=FALSE)
length(which(ccc%in%bbb)) ### VARIABLES IN dat_inv AND gwas MATCH


# GRS CALCULATION #

rsid<-as.character(sub("\\_.*", "", rownames(dat_inv)))
tested_allele_moba<-sub('.*_', '', rownames(dat_inv))
moba_info<-as.data.frame(cbind(rsid,tested_allele_moba))
moba_info$tested_allele_moba<-tolower(moba_info$tested_allele_moba)
gwas<-merge2(gwas,moba_info,by.id=c("rsid"),all.x=TRUE,sort=FALSE)
gwas<-na.omit(gwas)

gwas$same_coding<-with(gwas,ifelse(tested_allele_moba==effect_allele,1,
                                   ifelse(tested_allele_moba!=effect_allele,0,NA)))

gwas$flip<-with(gwas,ifelse(beta>=0 & same_coding==1,0,
                            ifelse(beta<0 & same_coding==0,0,
                                   ifelse(beta>=0 & same_coding==0,1,
                                          ifelse(beta<0 & same_coding==1,1,NA)))))

gwas$maf<-with(gwas,ifelse(same_coding==1,eaf,
                           ifelse(same_coding==0,1-eaf,NA)))
gwas$coef<-with(gwas,ifelse(beta>0,beta,
                            ifelse(beta<0,-beta,NA)))

n_snps<-dim(dat)[2]-1
colnames(dat)<-sub("\\_.*", "", colnames(dat))
dat<-rename.vars(dat,from=c("IID"),to=c("id"))

vars01<-colnames(dat)[2:dim(dat)[2]]
vars02<-paste(vars01,"_flip",sep="")
vars03<-paste(vars01,"_coef",sep="")
vars04<-paste(vars01,"_weight",sep="")

for(i in 1:length(vars01))
  
{
  dat[,vars02[i]]<-gwas[gwas$rsid==vars01[i],"flip"]
  dat[,vars03[i]]<-gwas[gwas$rsid==vars01[i],"coef"]
  dat[,vars04[i]]<-with(dat,ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==0,0*dat[,vars03[i]],
                                   ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==1,1*dat[,vars03[i]],
                                          ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==2,2*dat[,vars03[i]],
                                                 ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==0,2*dat[,vars03[i]],
                                                        ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==1,1*dat[,vars03[i]],
                                                               ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==2,0*dat[,vars03[i]],NA)))))))
}

dat_coef<-dat[,vars03]
dat$coef_sum<-(dat_coef %>%
                 replace(is.na(.), 0) %>%
                 mutate(coef_sum = rowSums(.)))$coef_sum[1]

vars05<-append("id",vars04,after=1)
dat_weight<-dat[,vars05]
dat_weight$weight_sum<-(dat_weight[,2:dim(dat_weight)[2]] %>%
                          replace(is.na(.), 0) %>%
                          mutate(weight_sum = rowSums(.)))$weight_sum
dat_weight<-dat_weight[,c("id","weight_sum")]
dat<-merge(dat,dat_weight,by="id")
dat$cigday_grs<-(dat$weight_sum/dat$coef_sum)*n_snps
dat<-dat[,c("id","cigday_grs")]
save(dat,file="./Data/cigday_grs.RData")


### AgeSmk ###
##############

dat<-as.data.frame(fread("./Data/agesmk_liu2019.raw",header=TRUE,sep="\t",sep2="\t"))
dat<-dat[,c(2,7:dim(dat)[2])]

gwas<-read.csv2("./Data/gwas_agesmk_liu2019.csv",header=TRUE,sep=";",dec=".")
names(gwas)<-tolower(names(gwas))
gwas<-rename.vars(gwas,
                  from=c("effect_allele_freq"),
                  to=c("eaf"))
gwas$effect_allele<-tolower(gwas$effect_allele)


# 9 SNPs IN MoBa AND GWAS #

bbb<-as.character(sort(gwas$rsid,decreasing=FALSE))
dat_inv<-header.true(data.frame(t(dat)))
ccc<-as.character(sort(sub("\\_.*", "", rownames(dat_inv)),decreasing=FALSE))
write.table(ccc,file="./Outputs/Descriptive/snps_agesmk_1smr.csv",sep=";", col.names=FALSE ,row.names=FALSE)
length(which(ccc%in%bbb)) ### VARIABLES IN dat_inv AND gwas MATCH


# GRS CALCULATION #

rsid<-as.character(sub("\\_.*", "", rownames(dat_inv)))
tested_allele_moba<-sub('.*_', '', rownames(dat_inv))
moba_info<-as.data.frame(cbind(rsid,tested_allele_moba))
moba_info$tested_allele_moba<-tolower(moba_info$tested_allele_moba)
gwas<-merge2(gwas,moba_info,by.id=c("rsid"),all.x=TRUE,sort=FALSE)
gwas<-na.omit(gwas)

gwas$same_coding<-with(gwas,ifelse(tested_allele_moba==effect_allele,1,
                                   ifelse(tested_allele_moba!=effect_allele,0,NA)))

gwas$flip<-with(gwas,ifelse(beta>=0 & same_coding==1,0,
                            ifelse(beta<0 & same_coding==0,0,
                                   ifelse(beta>=0 & same_coding==0,1,
                                          ifelse(beta<0 & same_coding==1,1,NA)))))

gwas$maf<-with(gwas,ifelse(same_coding==1,eaf,
                           ifelse(same_coding==0,1-eaf,NA)))
gwas$coef<-with(gwas,ifelse(beta>0,beta,
                            ifelse(beta<0,-beta,NA)))

n_snps<-dim(dat)[2]-1
colnames(dat)<-sub("\\_.*", "", colnames(dat))
dat<-rename.vars(dat,from=c("IID"),to=c("id"))

vars01<-colnames(dat)[2:dim(dat)[2]]
vars02<-paste(vars01,"_flip",sep="")
vars03<-paste(vars01,"_coef",sep="")
vars04<-paste(vars01,"_weight",sep="")

for(i in 1:length(vars01))
  
{
  dat[,vars02[i]]<-gwas[gwas$rsid==vars01[i],"flip"]
  dat[,vars03[i]]<-gwas[gwas$rsid==vars01[i],"coef"]
  dat[,vars04[i]]<-with(dat,ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==0,0*dat[,vars03[i]],
                                   ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==1,1*dat[,vars03[i]],
                                          ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==2,2*dat[,vars03[i]],
                                                 ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==0,2*dat[,vars03[i]],
                                                        ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==1,1*dat[,vars03[i]],
                                                               ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==2,0*dat[,vars03[i]],NA)))))))
}

dat_coef<-dat[,vars03]
dat$coef_sum<-(dat_coef %>%
                 replace(is.na(.), 0) %>%
                 mutate(coef_sum = rowSums(.)))$coef_sum[1]

vars05<-append("id",vars04,after=1)
dat_weight<-dat[,vars05]
dat_weight$weight_sum<-(dat_weight[,2:dim(dat_weight)[2]] %>%
                          replace(is.na(.), 0) %>%
                          mutate(weight_sum = rowSums(.)))$weight_sum
dat_weight<-dat_weight[,c("id","weight_sum")]
dat<-merge(dat,dat_weight,by="id")
dat$agesmk_grs<-(dat$weight_sum/dat$coef_sum)*n_snps
dat<-dat[,c("id","agesmk_grs")]
save(dat,file="./Data/agesmk_grs.RData")


### SmkCes ###
##############

dat<-as.data.frame(fread("./Data/smkces_liu2019.raw",header=TRUE,sep="\t",sep2="\t"))
dat<-dat[,c(2,7:dim(dat)[2])]

gwas<-read.csv2("./Data/gwas_smkces_liu2019.csv",header=TRUE,sep=";",dec=".")
names(gwas)<-tolower(names(gwas))
gwas<-rename.vars(gwas,
                  from=c("effect_allele_freq"),
                  to=c("eaf"))
gwas$effect_allele<-tolower(gwas$effect_allele)


# 22 SNPs IN MoBa AND GWAS #

bbb<-as.character(sort(gwas$rsid,decreasing=FALSE))
dat_inv<-header.true(data.frame(t(dat)))
ccc<-as.character(sort(sub("\\_.*", "", rownames(dat_inv)),decreasing=FALSE))
write.table(ccc,file="./Outputs/Descriptive/snps_smkces_1smr.csv",sep=";", col.names=FALSE ,row.names=FALSE)
length(which(ccc%in%bbb)) ### VARIABLES IN dat_inv AND gwas MATCH


# GRS CALCULATION #

rsid<-as.character(sub("\\_.*", "", rownames(dat_inv)))
tested_allele_moba<-sub('.*_', '', rownames(dat_inv))
moba_info<-as.data.frame(cbind(rsid,tested_allele_moba))
moba_info$tested_allele_moba<-tolower(moba_info$tested_allele_moba)
gwas<-merge2(gwas,moba_info,by.id=c("rsid"),all.x=TRUE,sort=FALSE)
gwas<-na.omit(gwas)

gwas$same_coding<-with(gwas,ifelse(tested_allele_moba==effect_allele,1,
                                   ifelse(tested_allele_moba!=effect_allele,0,NA)))

gwas$flip<-with(gwas,ifelse(beta>=0 & same_coding==1,0,
                            ifelse(beta<0 & same_coding==0,0,
                                   ifelse(beta>=0 & same_coding==0,1,
                                          ifelse(beta<0 & same_coding==1,1,NA)))))

gwas$maf<-with(gwas,ifelse(same_coding==1,eaf,
                           ifelse(same_coding==0,1-eaf,NA)))
gwas$coef<-with(gwas,ifelse(beta>0,beta,
                            ifelse(beta<0,-beta,NA)))

n_snps<-dim(dat)[2]-1
colnames(dat)<-sub("\\_.*", "", colnames(dat))
dat<-rename.vars(dat,from=c("IID"),to=c("id"))

vars01<-colnames(dat)[2:dim(dat)[2]]
vars02<-paste(vars01,"_flip",sep="")
vars03<-paste(vars01,"_coef",sep="")
vars04<-paste(vars01,"_weight",sep="")

for(i in 1:length(vars01))
  
{
  dat[,vars02[i]]<-gwas[gwas$rsid==vars01[i],"flip"]
  dat[,vars03[i]]<-gwas[gwas$rsid==vars01[i],"coef"]
  dat[,vars04[i]]<-with(dat,ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==0,0*dat[,vars03[i]],
                                   ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==1,1*dat[,vars03[i]],
                                          ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==2,2*dat[,vars03[i]],
                                                 ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==0,2*dat[,vars03[i]],
                                                        ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==1,1*dat[,vars03[i]],
                                                               ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==2,0*dat[,vars03[i]],NA)))))))
}

dat_coef<-dat[,vars03]
dat$coef_sum<-(dat_coef %>%
                 replace(is.na(.), 0) %>%
                 mutate(coef_sum = rowSums(.)))$coef_sum[1]

vars05<-append("id",vars04,after=1)
dat_weight<-dat[,vars05]
dat_weight$weight_sum<-(dat_weight[,2:dim(dat_weight)[2]] %>%
                          replace(is.na(.), 0) %>%
                          mutate(weight_sum = rowSums(.)))$weight_sum
dat_weight<-dat_weight[,c("id","weight_sum")]
dat<-merge(dat,dat_weight,by="id")
dat$smkces_grs<-(dat$weight_sum/dat$coef_sum)*n_snps
dat<-dat[,c("id","smkces_grs")]
save(dat,file="./Data/smkces_grs.RData")


### Risktk ###
##############

dat<-as.data.frame(fread("./Data/risktk_strawbridge2018.raw",header=TRUE,sep="\t",sep2="\t"))
dat<-dat[,c(2,7:dim(dat)[2])]

gwas<-read.csv2("./Data/gwas_strawbridge_2018.csv",header=TRUE,sep=";",dec=".")
names(gwas)<-tolower(names(gwas))
gwas<-rename.vars(gwas,
                  from=c("effect_allele_freq"),
                  to=c("eaf"))
gwas$effect_allele<-tolower(gwas$effect_allele)


# 2 SNPs IN MoBa AND GWAS #

bbb<-as.character(sort(gwas$rsid,decreasing=FALSE))
dat_inv<-header.true(data.frame(t(dat)))
ccc<-as.character(sort(sub("\\_.*", "", rownames(dat_inv)),decreasing=FALSE))
write.table(ccc,file="./Outputs/Descriptive/snps_risktk_1smr.csv",sep=";", col.names=FALSE ,row.names=FALSE)
length(which(ccc%in%bbb)) ### VARIABLES IN dat_inv AND gwas MATCH


# GRS CALCULATION #

rsid<-as.character(sub("\\_.*", "", rownames(dat_inv)))
tested_allele_moba<-sub('.*_', '', rownames(dat_inv))
moba_info<-as.data.frame(cbind(rsid,tested_allele_moba))
moba_info$tested_allele_moba<-tolower(moba_info$tested_allele_moba)
gwas<-merge2(gwas,moba_info,by.id=c("rsid"),all.x=TRUE,sort=FALSE)
gwas<-na.omit(gwas)

gwas$same_coding<-with(gwas,ifelse(tested_allele_moba==effect_allele,1,
                                   ifelse(tested_allele_moba!=effect_allele,0,NA)))

gwas$flip<-with(gwas,ifelse(beta>=0 & same_coding==1,0,
                            ifelse(beta<0 & same_coding==0,0,
                                   ifelse(beta>=0 & same_coding==0,1,
                                          ifelse(beta<0 & same_coding==1,1,NA)))))

gwas$maf<-with(gwas,ifelse(same_coding==1,eaf,
                           ifelse(same_coding==0,1-eaf,NA)))
gwas$coef<-with(gwas,ifelse(beta>0,beta,
                            ifelse(beta<0,-beta,NA)))

n_snps<-dim(dat)[2]-1
colnames(dat)<-sub("\\_.*", "", colnames(dat))
dat<-rename.vars(dat,from=c("IID"),to=c("id"))

vars01<-colnames(dat)[2:dim(dat)[2]]
vars02<-paste(vars01,"_flip",sep="")
vars03<-paste(vars01,"_coef",sep="")
vars04<-paste(vars01,"_weight",sep="")

for(i in 1:length(vars01))
  
{
  dat[,vars02[i]]<-gwas[gwas$rsid==vars01[i],"flip"]
  dat[,vars03[i]]<-gwas[gwas$rsid==vars01[i],"coef"]
  dat[,vars04[i]]<-with(dat,ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==0,0*dat[,vars03[i]],
                                   ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==1,1*dat[,vars03[i]],
                                          ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==2,2*dat[,vars03[i]],
                                                 ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==0,2*dat[,vars03[i]],
                                                        ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==1,1*dat[,vars03[i]],
                                                               ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==2,0*dat[,vars03[i]],NA)))))))
}

dat_coef<-dat[,vars03]
dat$coef_sum<-(dat_coef %>%
                 replace(is.na(.), 0) %>%
                 mutate(coef_sum = rowSums(.)))$coef_sum[1]

vars05<-append("id",vars04,after=1)
dat_weight<-dat[,vars05]
dat_weight$weight_sum<-(dat_weight[,2:dim(dat_weight)[2]] %>%
                          replace(is.na(.), 0) %>%
                          mutate(weight_sum = rowSums(.)))$weight_sum
dat_weight<-dat_weight[,c("id","weight_sum")]
dat<-merge(dat,dat_weight,by="id")
dat$risktk_grs<-(dat$weight_sum/dat$coef_sum)*n_snps
dat<-dat[,c("id","risktk_grs")]
save(dat,file="./Data/risktk_grs.RData")


### BMI ###
###########

dat<-as.data.frame(fread("./Data/bmi_yengo2018.raw",header=TRUE,sep="\t",sep2="\t"))
dat<-dat[,c(2,7:dim(dat)[2])]

gwas<-read.csv2("./Data/gwas_bmi_yengo2018.csv",header=TRUE,sep=",",dec=".")
names(gwas)<-tolower(names(gwas))
gwas$beta<-NULL
gwas<-rename.vars(gwas,
                  from=c("snp","tested_allele","freq_tested_allele_in_hrs","beta_cojo"),
                  to=c("rsid","effect_allele","eaf","beta"))
gwas$effect_allele<-tolower(gwas$effect_allele)


# 908 SNPs IN MoBa AND GWAS #

bbb<-as.character(sort(gwas$rsid,decreasing=FALSE))
dat_inv<-header.true(data.frame(t(dat)))
ccc<-as.character(sort(sub("\\_.*", "", rownames(dat_inv)),decreasing=FALSE))
write.table(ccc,file="./Outputs/Descriptive/snps_bmi_1smr.csv",sep=";", col.names=FALSE ,row.names=FALSE)
length(which(ccc%in%bbb)) ### VARIABLES IN dat_inv AND gwas MATCH


# GRS CALCULATION #

rsid<-as.character(sub("\\_.*", "", rownames(dat_inv)))
tested_allele_moba<-sub('.*_', '', rownames(dat_inv))
moba_info<-as.data.frame(cbind(rsid,tested_allele_moba))
moba_info$tested_allele_moba<-tolower(moba_info$tested_allele_moba)
gwas<-merge2(gwas,moba_info,by.id=c("rsid"),all.x=TRUE,sort=FALSE)
gwas<-na.omit(gwas)

gwas$same_coding<-with(gwas,ifelse(tested_allele_moba==effect_allele,1,
                                   ifelse(tested_allele_moba!=effect_allele,0,NA)))

gwas$flip<-with(gwas,ifelse(beta>=0 & same_coding==1,0,
                            ifelse(beta<0 & same_coding==0,0,
                                   ifelse(beta>=0 & same_coding==0,1,
                                          ifelse(beta<0 & same_coding==1,1,NA)))))

gwas$maf<-with(gwas,ifelse(same_coding==1,eaf,
                           ifelse(same_coding==0,1-eaf,NA)))
gwas$coef<-with(gwas,ifelse(beta>0,beta,
                            ifelse(beta<0,-beta,NA)))

n_snps<-dim(dat)[2]-1
colnames(dat)<-sub("\\_.*", "", colnames(dat))
dat<-rename.vars(dat,from=c("IID"),to=c("id"))

vars01<-colnames(dat)[2:dim(dat)[2]]
vars02<-paste(vars01,"_flip",sep="")
vars03<-paste(vars01,"_coef",sep="")
vars04<-paste(vars01,"_weight",sep="")

for(i in 1:length(vars01))
  
{
  dat[,vars02[i]]<-gwas[gwas$rsid==vars01[i],"flip"]
  dat[,vars03[i]]<-gwas[gwas$rsid==vars01[i],"coef"]
  dat[,vars04[i]]<-with(dat,ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==0,0*dat[,vars03[i]],
                                   ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==1,1*dat[,vars03[i]],
                                          ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==2,2*dat[,vars03[i]],
                                                 ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==0,2*dat[,vars03[i]],
                                                        ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==1,1*dat[,vars03[i]],
                                                               ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==2,0*dat[,vars03[i]],NA)))))))
}

dat_coef<-dat[,vars03]
dat$coef_sum<-(dat_coef %>%
                 replace(is.na(.), 0) %>%
                 mutate(coef_sum = rowSums(.)))$coef_sum[1]

vars05<-append("id",vars04,after=1)
dat_weight<-dat[,vars05]
dat_weight$weight_sum<-(dat_weight[,2:dim(dat_weight)[2]] %>%
                          replace(is.na(.), 0) %>%
                          mutate(weight_sum = rowSums(.)))$weight_sum
dat_weight<-dat_weight[,c("id","weight_sum")]
dat<-merge(dat,dat_weight,by="id")
dat$bmi_grs<-(dat$weight_sum/dat$coef_sum)*n_snps
dat<-dat[,c("id","bmi_grs")]
save(dat,file="./Data/bmi_grs.RData")
dat<-NULL
dat_coef<-NULL
dat_inv<-NULL
dat_weight<-NULL


### EduYears ###
################

dat<-as.data.frame(fread("./Data/eduyears_liu2019.raw",header=TRUE,sep="\t",sep2="\t"))
dat<-dat[,c(2,7:dim(dat)[2])]

gwas<-read.csv2("./Data/gwas_eduyears_liu2019.csv",header=TRUE,sep=";",dec=".")
names(gwas)<-tolower(names(gwas))
gwas<-rename.vars(gwas,
                  from=c("effect_allele_freq"),
                  to=c("eaf"))
gwas$effect_allele<-tolower(gwas$effect_allele)


# 1176 SNPs IN MoBa AND GWAS #

bbb<-as.character(sort(gwas$rsid,decreasing=FALSE))
dat_inv<-header.true(data.frame(t(dat)))
ccc<-as.character(sort(sub("\\_.*", "", rownames(dat_inv)),decreasing=FALSE))
write.table(ccc,file="./Outputs/Descriptive/snps_eduyears_1smr.csv",sep=";", col.names=FALSE ,row.names=FALSE)
length(which(ccc%in%bbb)) ### VARIABLES IN dat_inv AND gwas MATCH


# GRS CALCULATION #

rsid<-as.character(sub("\\_.*", "", rownames(dat_inv)))
tested_allele_moba<-sub('.*_', '', rownames(dat_inv))
moba_info<-as.data.frame(cbind(rsid,tested_allele_moba))
moba_info$tested_allele_moba<-tolower(moba_info$tested_allele_moba)
gwas<-merge2(gwas,moba_info,by.id=c("rsid"),all.x=TRUE,sort=FALSE)
gwas<-na.omit(gwas)

gwas$same_coding<-with(gwas,ifelse(tested_allele_moba==effect_allele,1,
                                   ifelse(tested_allele_moba!=effect_allele,0,NA)))

gwas$flip<-with(gwas,ifelse(beta>=0 & same_coding==1,0,
                            ifelse(beta<0 & same_coding==0,0,
                                   ifelse(beta>=0 & same_coding==0,1,
                                          ifelse(beta<0 & same_coding==1,1,NA)))))

gwas$maf<-with(gwas,ifelse(same_coding==1,eaf,
                           ifelse(same_coding==0,1-eaf,NA)))
gwas$coef<-with(gwas,ifelse(beta>0,beta,
                            ifelse(beta<0,-beta,NA)))

n_snps<-dim(dat)[2]-1
colnames(dat)<-sub("\\_.*", "", colnames(dat))
dat<-rename.vars(dat,from=c("IID"),to=c("id"))

vars01<-colnames(dat)[2:dim(dat)[2]]
vars02<-paste(vars01,"_flip",sep="")
vars03<-paste(vars01,"_coef",sep="")
vars04<-paste(vars01,"_weight",sep="")

for(i in 1:length(vars01))
  
{
  dat[,vars02[i]]<-gwas[gwas$rsid==vars01[i],"flip"]
  dat[,vars03[i]]<-gwas[gwas$rsid==vars01[i],"coef"]
  dat[,vars04[i]]<-with(dat,ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==0,0*dat[,vars03[i]],
                                   ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==1,1*dat[,vars03[i]],
                                          ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==2,2*dat[,vars03[i]],
                                                 ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==0,2*dat[,vars03[i]],
                                                        ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==1,1*dat[,vars03[i]],
                                                               ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==2,0*dat[,vars03[i]],NA)))))))
}

dat_coef<-dat[,vars03]
dat$coef_sum<-(dat_coef %>%
                 replace(is.na(.), 0) %>%
                 mutate(coef_sum = rowSums(.)))$coef_sum[1]

vars05<-append("id",vars04,after=1)
dat_weight<-dat[,vars05]
dat_weight$weight_sum<-(dat_weight[,2:dim(dat_weight)[2]] %>%
                          replace(is.na(.), 0) %>%
                          mutate(weight_sum = rowSums(.)))$weight_sum
dat_weight<-dat_weight[,c("id","weight_sum")]
dat<-merge(dat,dat_weight,by="id")
dat$eduyears_grs<-(dat$weight_sum/dat$coef_sum)*n_snps
dat<-dat[,c("id","eduyears_grs")]
save(dat,file="./Data/eduyears_grs.RData")
dat<-NULL
dat_coef<-NULL
dat_inv<-NULL
dat_weight<-NULL
gwas<-NULL
moba_info<-NULL


######################################
### GENERATION OF WORKING DATABASE ###
######################################

### BASIC FORMAT OF VARIABLES ###

load("./Data/d_moba_mr.RData")
load("./Data/d_moba_bmi_2021.RData")
d_moba_bmi<-rename.vars(d_moba_bmi,
                        from=c("bmi"),
                        to=c("bmi_ok"))
dat<-merge(d_moba_mr,d_moba_bmi,by=c("P_ID_2824"),all.x=TRUE,sort=FALSE)

names(dat)<-tolower(names(dat))
dat$role<-with(dat,ifelse(role=="Mother",0,1))
attr(dat$role,"value.labels")<-c("Mother"=0, "Father"=1)
dat<-rename.vars(dat,
                 from=c("p_id_2824","smoking_number_daily_cig","covid_diagnosis","number_cohabitants","work_situation"),
                 to=c("id_2824","smoking_number","covid","cohabs","work"))
dat$edu<-with(dat,ifelse(edu=="< High school",0,
                         ifelse(edu=="High school",1,
                                ifelse(edu=="College ???4 years",2,
                                       ifelse(edu==">4 years college",3,NA)))))
dat$edu<-with(dat,ifelse(is.na(edu),2,edu))
attr(dat$edu,"value.labels")<-c("< High school"=0, "High school"=1, "College <=4 years"=2, ">4 years college"=3)
dat$eduyears<-with(dat,ifelse(edu==0,10,
                              ifelse(edu==1,13,
                                     ifelse(edu==2,19,
                                            ifelse(edu==3,20,NA)))))
dat$region<-with(dat,ifelse(region=="OsloViken",1,
                            ifelse(region=="South",2,
                                   ifelse(region=="West",3,
                                          ifelse(region=="Middle",4,
                                                 ifelse(region=="North",5,NA))))))
dat$region<-with(dat,ifelse(is.na(region),6,region))
attr(dat$region,"value.labels")<-c("Oslo/Viken"=1, "South"=2, "West"=3, "Middle"=4, "North"=5, "Unknown"=6)
dat$work<-with(dat,ifelse(work=="no_or_other",1,
                          ifelse(work=="home_office",2,
                                 ifelse(work=="lostjob_or_sickleave",3,NA))))
dat$work<-with(dat,ifelse(is.na(work),4,work))
attr(dat$work,"value.labels")<-c("no_or_other"=1, "home_office"=2, "lostjob_or_sickleave"=3, "unknown"=4)
dat$covid<-with(dat,ifelse(is.na(covid),0,covid))

moms<-dat[dat$role==0,]
moms<-rename.vars(moms, from=c("id_2824"), to=c("m_id_2824"))
moms<-moms[which(!duplicated(moms$m_id_2824)),]
dads<-dat[dat$role==1,]
dads<-rename.vars(dads, from=c("id_2824"), to=c("f_id_2824"))
dads<-dads[which(!duplicated(dads$f_id_2824)),]

moba<-spss.get("N:/durable/Data/PDB2824_transferedfiles/MoBaVersjon12_d30092021/PDB2824_SV_INFO_v12.sav",
               use.value.labels=FALSE,to.data.frame=TRUE,allow="_")
names(moba)<-tolower(names(moba))
moba$m_id_2824<-gsub(" ", "", moba$m_id_2824)
moba$f_id_2824<-gsub(" ", "", moba$f_id_2824)
moba_mom<-moba[,c("m_id_2824","preg_id_2824","faar")]
moba_mom<-moba_mom[order(moba_mom$m_id_2824,-abs(moba_mom$faar)),]
moba_mom<-moba_mom[!duplicated(moba_mom$m_id_2824),]
moba_mom$faar<-NULL
moba_dad<-moba[,c("f_id_2824","preg_id_2824","faar")]
moba_dad<-moba_dad[order(moba_dad$f_id_2824,-abs(moba_dad$faar)),]
moba_dad<-moba_dad[!duplicated(moba_dad$f_id_2824),]
moba_dad$faar<-NULL
moms<-merge(moms,moba_mom,by=c("m_id_2824"),all.x=TRUE,sort=FALSE)
dads<-merge(dads,moba_dad,by=c("f_id_2824"),all.x=TRUE,sort=FALSE)


### BMI AND AGESMK FOR MOTHERS FROM MOBA DATABASE (BMI IN THE LATEST PREGNANCY) ###

moba<-spss.get("N:/durable/Data/PDB2824_transferedfiles/MoBaVersjon12_d30092021/PDB2824_Skjema1_v12.sav",
               use.value.labels=FALSE,to.data.frame=TRUE,allow="_")
names(moba)<-tolower(names(moba))
moba$agesmk<-moba$aa1362
moba$bmi_mom<-moba$aa85/((moba$aa87/100)^2)
moba$bmi_dad<-moba$aa89/((moba$aa88/100)^2)
moba_mom<-moba[,c("preg_id_2824","bmi_mom","agesmk")]
moms<-merge(moms,moba_mom,by=c("preg_id_2824"),all.x=TRUE,sort=FALSE)
moba_dad<-moba[,c("preg_id_2824","bmi_dad","agesmk")]
moba_dad$agesmk<-NA
dads<-merge2(dads,moba_dad,by.id=c("preg_id_2824"),all.x=TRUE,sort=FALSE)

moms$bmi_mom<-with(moms,ifelse(!is.na(bmi_ok),bmi_ok,bmi_mom))
dads$bmi_dad<-with(dads,ifelse(!is.na(bmi_ok),bmi_ok,bmi_dad))
moms$bmi_mom<-with(moms,ifelse(bmi_mom<13 | bmi_mom>60,NA,bmi_mom))
dads$bmi_dad<-with(dads,ifelse(bmi_dad<13 | bmi_dad>60,NA,bmi_dad))
moms$bmi_ok<-NULL
dads$bmi_ok<-NULL

moba<-NULL
moba_mom<-NULL
moba_dad<-NULL
d_moba_mr<-NULL
d_moba_bmi<-NULL


### MERGING GRSs ###

load("./Data/smkinit_grs.RData")
smkinit_grs<-dat
smkinit_grs$genotype<-1
genotype<-smkinit_grs[,c("id","genotype")]

mom_gen<-spss.get("N:/durable/Data/PDB2824_updated_linkagefiles/pdb2824_gwas_20220516/20220516_MoBaGeneticsTot_Mother_PDB2824.sav",
                  use.value.labels=FALSE,to.data.frame=TRUE,allow="_")
names(mom_gen)<-tolower(names(mom_gen))
mom_gen$m_id_2824<-gsub(" ", "", mom_gen$m_id_2824)
mom_gen$sentrixid_mom<-gsub(" ", "", mom_gen$sentrix_id)
mom_gen$batch_mom<-gsub(" ", "", mom_gen$batch)
mom_gen<-mom_gen[,c("m_id_2824","sentrixid_mom","batch_mom")]
names(genotype)<-c("sentrixid_mom","genotype")
mom_gen<-merge2(mom_gen,genotype,by.id=c("sentrixid_mom"),all.x=TRUE,sort=FALSE)
mom_gen$genotype<-with(mom_gen,ifelse(is.na(genotype),0,genotype))
mom_gen<-mom_gen[order(mom_gen$m_id_2824,-abs(mom_gen$genotype)),]
mom_gen<-mom_gen[!duplicated(mom_gen$m_id_2824),]
mom_gen$genotype<-NULL
moms<-merge2(moms,mom_gen,by.id=c("m_id_2824"),all.x=TRUE,sort=FALSE)
mom_gen<-NULL

dad_gen<-spss.get("N:/durable/Data/PDB2824_updated_linkagefiles/pdb2824_gwas_20220516/20220516_MoBaGeneticsTot_Father_PDB2824.sav",
                  use.value.labels=FALSE,to.data.frame=TRUE,allow="_")
names(dad_gen)<-tolower(names(dad_gen))
dad_gen$f_id_2824<-gsub(" ", "", dad_gen$f_id_2824)
dad_gen$sentrixid_dad<-gsub(" ", "", dad_gen$sentrix_id)
dad_gen$batch_dad<-gsub(" ", "", dad_gen$batch)
dad_gen<-dad_gen[,c("f_id_2824","sentrixid_dad","batch_dad")]
names(genotype)<-c("sentrixid_dad","genotype")
dad_gen<-merge2(dad_gen,genotype,by.id=c("sentrixid_dad"),all.x=TRUE,sort=FALSE)
dad_gen$genotype<-with(dad_gen,ifelse(is.na(genotype),0,genotype))
dad_gen<-dad_gen[order(dad_gen$f_id_2824,-abs(dad_gen$genotype)),]
dad_gen<-dad_gen[!duplicated(dad_gen$f_id_2824),]
dad_gen$genotype<-NULL
dads<-merge2(dads,dad_gen,by.id=c("f_id_2824"),all.x=TRUE,sort=FALSE)
dad_gen<-NULL
genotype<-NULL

load("./Data/smkinit_grs.RData")
smkinit_grs<-dat
load("./Data/cigday_grs.RData")
cigday_grs<-dat
load("./Data/agesmk_grs.RData")
agesmk_grs<-dat
load("./Data/smkces_grs.RData")
smkces_grs<-dat
load("./Data/lifesmk_grs.RData")
lifesmk_grs<-dat
load("./Data/risktk_grs.RData")
risktk_grs<-dat
load("./Data/bmi_grs.RData")
bmi_grs<-dat
load("./Data/eduyears_grs.RData")
eduyears_grs<-dat
dat<-NULL

smkinit_grs<-rename.vars(smkinit_grs,from=c("id","smkinit_grs"),to=c("sentrixid_mom","smkinit_grs_mom"))
moms<-merge2(moms,smkinit_grs,by.id=c("sentrixid_mom"),all.x=TRUE,sort=FALSE)
smkinit_grs<-rename.vars(smkinit_grs,from=c("sentrixid_mom","smkinit_grs_mom"),to=c("sentrixid_dad","smkinit_grs_dad"))
dads<-merge2(dads,smkinit_grs,by.id=c("sentrixid_dad"),all.x=TRUE,sort=FALSE)
smkinit_grs<-NULL

cigday_grs<-rename.vars(cigday_grs,from=c("id","cigday_grs"),to=c("sentrixid_mom","cigday_grs_mom"))
moms<-merge2(moms,cigday_grs,by.id=c("sentrixid_mom"),all.x=TRUE,sort=FALSE)
cigday_grs<-rename.vars(cigday_grs,from=c("sentrixid_mom","cigday_grs_mom"),to=c("sentrixid_dad","cigday_grs_dad"))
dads<-merge2(dads,cigday_grs,by.id=c("sentrixid_dad"),all.x=TRUE,sort=FALSE)
cigday_grs<-NULL

agesmk_grs<-rename.vars(agesmk_grs,from=c("id","agesmk_grs"),to=c("sentrixid_mom","agesmk_grs_mom"))
moms<-merge2(moms,agesmk_grs,by.id=c("sentrixid_mom"),all.x=TRUE,sort=FALSE)
agesmk_grs<-rename.vars(agesmk_grs,from=c("sentrixid_mom","agesmk_grs_mom"),to=c("sentrixid_dad","agesmk_grs_dad"))
dads<-merge2(dads,agesmk_grs,by.id=c("sentrixid_dad"),all.x=TRUE,sort=FALSE)
agesmk_grs<-NULL

smkces_grs<-rename.vars(smkces_grs,from=c("id","smkces_grs"),to=c("sentrixid_mom","smkces_grs_mom"))
moms<-merge2(moms,smkces_grs,by.id=c("sentrixid_mom"),all.x=TRUE,sort=FALSE)
smkces_grs<-rename.vars(smkces_grs,from=c("sentrixid_mom","smkces_grs_mom"),to=c("sentrixid_dad","smkces_grs_dad"))
dads<-merge2(dads,smkces_grs,by.id=c("sentrixid_dad"),all.x=TRUE,sort=FALSE)
smkces_grs<-NULL

lifesmk_grs<-rename.vars(lifesmk_grs,from=c("id","lifesmk_grs"),to=c("sentrixid_mom","lifesmk_grs_mom"))
moms<-merge2(moms,lifesmk_grs,by.id=c("sentrixid_mom"),all.x=TRUE,sort=FALSE)
lifesmk_grs<-rename.vars(lifesmk_grs,from=c("sentrixid_mom","lifesmk_grs_mom"),to=c("sentrixid_dad","lifesmk_grs_dad"))
dads<-merge2(dads,lifesmk_grs,by.id=c("sentrixid_dad"),all.x=TRUE,sort=FALSE)
lifesmk_grs<-NULL

risktk_grs<-rename.vars(risktk_grs,from=c("id","risktk_grs"),to=c("sentrixid_mom","risktk_grs_mom"))
moms<-merge2(moms,risktk_grs,by.id=c("sentrixid_mom"),all.x=TRUE,sort=FALSE)
risktk_grs<-rename.vars(risktk_grs,from=c("sentrixid_mom","risktk_grs_mom"),to=c("sentrixid_dad","risktk_grs_dad"))
dads<-merge2(dads,risktk_grs,by.id=c("sentrixid_dad"),all.x=TRUE,sort=FALSE)
risktk_grs<-NULL

bmi_grs<-rename.vars(bmi_grs,from=c("id","bmi_grs"),to=c("sentrixid_mom","bmi_grs_mom"))
moms<-merge2(moms,bmi_grs,by.id=c("sentrixid_mom"),all.x=TRUE,sort=FALSE)
bmi_grs<-rename.vars(bmi_grs,from=c("sentrixid_mom","bmi_grs_mom"),to=c("sentrixid_dad","bmi_grs_dad"))
dads<-merge2(dads,bmi_grs,by.id=c("sentrixid_dad"),all.x=TRUE,sort=FALSE)
bmi_grs<-NULL

eduyears_grs<-rename.vars(eduyears_grs,from=c("id","eduyears_grs"),to=c("sentrixid_mom","eduyears_grs_mom"))
moms<-merge2(moms,eduyears_grs,by.id=c("sentrixid_mom"),all.x=TRUE,sort=FALSE)
eduyears_grs<-rename.vars(eduyears_grs,from=c("sentrixid_mom","eduyears_grs_mom"),to=c("sentrixid_dad","eduyears_grs_dad"))
dads<-merge2(dads,eduyears_grs,by.id=c("sentrixid_dad"),all.x=TRUE,sort=FALSE)
eduyears_grs<-NULL

length(which(!is.na(moms$smkinit_grs_mom)))
length(which(!is.na(moms$cigday_grs_mom)))
length(which(!is.na(moms$agesmk_grs_mom)))
length(which(!is.na(moms$smkces_grs_mom)))
length(which(!is.na(moms$lifesmk_grs_mom)))
length(which(!is.na(moms$risktk_grs_mom)))
length(which(!is.na(moms$bmi_grs_mom)))
length(which(!is.na(moms$eduyears_grs_mom)))
length(which(!is.na(dads$smkinit_grs_dad)))
length(which(!is.na(dads$cigday_grs_dad)))
length(which(!is.na(dads$agesmk_grs_dad)))
length(which(!is.na(dads$smkces_grs_dad)))
length(which(!is.na(dads$lifesmk_grs_dad)))
length(which(!is.na(dads$risktk_grs_dad)))
length(which(!is.na(dads$bmi_grs_dad)))
length(which(!is.na(dads$eduyears_grs_dad)))


### MERGING PRINCIPAL COMPONENTS ###

pcs<-as.data.frame(read.delim("N:/durable/Data/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc-cov-noMoBaIDs.txt",
                              header=TRUE,sep="\t"))[,c("SENTRIXID","genotyping_batch_num","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10",
                                                        "PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20")]
names(pcs)<-c("sentrixid_mom","genotype_batch_mom","pc01_mom","pc02_mom","pc03_mom","pc04_mom","pc05_mom","pc06_mom","pc07_mom","pc08_mom","pc09_mom","pc10_mom",
              "pc11_mom","pc12_mom","pc13_mom","pc14_mom","pc15_mom","pc16_mom","pc17_mom","pc18_mom","pc19_mom","pc20_mom")
pcs$sentrixid_mom<-as.character(pcs$sentrixid_mom)
moms<-merge2(moms,pcs,by.id=c("sentrixid_mom"),all.x=TRUE,sort=FALSE)

pcs<-rename.vars(pcs,
                 from=c("sentrixid_mom","genotype_batch_mom","pc01_mom","pc02_mom","pc03_mom","pc04_mom","pc05_mom","pc06_mom","pc07_mom","pc08_mom","pc09_mom","pc10_mom",
                        "pc11_mom","pc12_mom","pc13_mom","pc14_mom","pc15_mom","pc16_mom","pc17_mom","pc18_mom","pc19_mom","pc20_mom"),
                 to=c("sentrixid_dad","genotype_batch_dad","pc01_dad","pc02_dad","pc03_dad","pc04_dad","pc05_dad","pc06_dad","pc07_dad","pc08_dad","pc09_dad","pc10_dad",
                      "pc11_dad","pc12_dad","pc13_dad","pc14_dad","pc15_dad","pc16_dad","pc17_dad","pc18_dad","pc19_dad","pc20_dad"))
pcs$sentrixid_dad<-as.character(pcs$sentrixid_dad)
dads<-merge2(dads,pcs,by.id=c("sentrixid_dad"),all.x=TRUE,sort=FALSE)

length(which(!is.na(moms$pc01_mom)))
length(which(!is.na(moms$pc02_mom)))
length(which(!is.na(dads$pc01_dad)))
length(which(!is.na(dads$pc02_dad)))
pcs<-NULL


### PARTICIPANTS WITHOUT CONSENT ###

excl<-spss.get("N:/durable/Data/PDB2824_transferedfiles/MoBaVersjon12_d30092021/PDB2824_SV_INFO_v12.sav",
               use.value.labels=FALSE,to.data.frame=TRUE,allow="_")
names(excl)<-tolower(names(excl))
excl$consent<-1
excl<-excl[,c("preg_id_2824","consent")]
moms<-merge2(moms,excl,by.id=c("preg_id_2824"),all.x=TRUE,sort=FALSE)
moms$consent<-with(moms,ifelse(is.na(consent),0,consent))
dads<-merge2(dads,excl,by.id=c("preg_id_2824"),all.x=TRUE,sort=FALSE)
dads$consent<-with(dads,ifelse(is.na(consent),0,consent))
excl<-NULL


### FINAL FORMAT OF DATA ###

moms$sentrixid_mom<-NULL
moms$batch_mom<-NULL
names(moms)<-c("preg_id_2824","id_2824","role","smkinit","smkces","cigday","covid","age","edu","cohabs",
               "region","work","eduyears","bmi","agesmk",
               "smkinit_grs","cigday_grs","agesmk_grs","smkces_grs","lifesmk_grs","risktk_grs","bmi_grs","eduyears_grs","genotype_batch",
               "pc01","pc02","pc03","pc04","pc05","pc06","pc07","pc08","pc09","pc10","pc11","pc12",
               "pc13","pc14","pc15","pc16","pc17","pc18","pc19","pc20","consent")
dads$sentrixid_dad<-NULL
dads$batch_dad<-NULL
names(dads)<-c("preg_id_2824","id_2824","role","smkinit","smkces","cigday","covid","age","edu","cohabs",
               "region","work","eduyears","bmi","agesmk",
               "smkinit_grs","cigday_grs","agesmk_grs","smkces_grs","lifesmk_grs","risktk_grs","bmi_grs","eduyears_grs","genotype_batch",
               "pc01","pc02","pc03","pc04","pc05","pc06","pc07","pc08","pc09","pc10","pc11","pc12",
               "pc13","pc14","pc15","pc16","pc17","pc18","pc19","pc20","consent")

dat<-rbind(moms,dads)
attr(dat$role,"value.labels")<-c("Mother"=0, "Father"=1)
attr(dat$edu,"value.labels")<-c("< High school"=0, "High school"=1, "College <=4 years"=2, ">4 years college"=3)
attr(dat$region,"value.labels")<-c("Oslo/Viken"=1, "South"=2, "West"=3, "Middle"=4, "North"=5, "Unknown"=6)
attr(dat$work,"value.labels")<-c("no_or_other"=1, "home_office"=2, "lostjob_or_sickleave"=3, "unknown"=4)
attr(moms$role,"value.labels")<-c("Mother"=0, "Father"=1)
attr(moms$edu,"value.labels")<-c("< High school"=0, "High school"=1, "College <=4 years"=2, ">4 years college"=3)
attr(moms$region,"value.labels")<-c("Oslo/Viken"=1, "South"=2, "West"=3, "Middle"=4, "North"=5, "Unknown"=6)
attr(moms$work,"value.labels")<-c("no_or_other"=1, "home_office"=2, "lostjob_or_sickleave"=3, "unknown"=4)
attr(dads$role,"value.labels")<-c("Mother"=0, "Father"=1)
attr(dads$edu,"value.labels")<-c("< High school"=0, "High school"=1, "College <=4 years"=2, ">4 years college"=3)
attr(dads$region,"value.labels")<-c("Oslo/Viken"=1, "South"=2, "West"=3, "Middle"=4, "North"=5, "Unknown"=6)
attr(dads$work,"value.labels")<-c("no_or_other"=1, "home_office"=2, "lostjob_or_sickleave"=3, "unknown"=4)

dat<-dat[!is.na(dat$smkinit_grs) & !is.na(dat$pc01) & dat$consent==1,]
save(dat,file="./Data/mr_smoking_covid_all.RData")


### STUDY FLOW CHART ###

excl<-spss.get("N:/durable/Data/PDB2824_transferedfiles/MoBaVersjon12_d30092021/PDB2824_SV_INFO_v12.sav",
               use.value.labels=FALSE,to.data.frame=TRUE,allow="_")
names(excl)<-tolower(names(excl))
excl$m_id_2824<-gsub(" ", "", excl$m_id_2824)
excl$m_id_2824<-with(excl,ifelse(m_id_2824=="",NA,m_id_2824))
excl$f_id_2824<-gsub(" ", "", excl$f_id_2824)
excl$f_id_2824<-with(excl,ifelse(f_id_2824=="",NA,f_id_2824))

# Total unique mothers in the last MoBa wave, consent OK (n = 95135)
length(which(!is.na(excl$m_id_2824) & !duplicated(excl$m_id_2824)))

# Total unique mothers in smoking-COVID sub-study, consent OK (n = 57397)
moms<-moms[moms$consent==1,] 
dim(moms)[1]

# Total unique mothers in smoking-COVID sub-study with genotype data, consent OK (n = 47506)
moms<-moms[!is.na(moms$smkinit_grs) & !is.na(moms$pc01),]
dim(moms)[1]
save(moms,file="./Data/mr_smoking_covid_moms.RData")

# Total unique fathers in the last MoBa wave, consent OK (n = 75109)
length(which(!is.na(excl$f_id_2824) & !duplicated(excl$f_id_2824)))

# Total unique fathers in smoking-COVID sub-study, consent OK (n = 39801)
dads<-dads[dads$consent==1,] 
dim(dads)[1]

# Total unique fathers in smoking-COVID sub-study with genotype data, consent OK (n = 28229)
dads<-dads[!is.na(dads$smkinit_grs) & !is.na(dads$pc01),]
dim(dads)[1]
save(dads,file="./Data/mr_smoking_covid_dads.RData")

excl<-NULL


### DESCRIPTIVE OF GRSs ###

dir.create("N:/durable/Projects/ALHE_smoking_covid/Outputs/Descriptive/normality")
setwd("N:/durable/Projects/ALHE_smoking_covid/Outputs/Descriptive/normality")

pdf(file="./smkinit_grs_all.pdf")
par(las=1,cex=1,mar=c(5,5,2,2))
plot(compareGroups(~smkinit_grs,data=dat))
dev.off()

pdf(file="./smkinit_grs_moms.pdf")
par(las=1,cex=1,mar=c(5,5,2,2))
plot(compareGroups(~smkinit_grs,data=moms))
dev.off()

pdf(file="./smkinit_grs_dads.pdf")
par(las=1,cex=1,mar=c(5,5,2,2))
plot(compareGroups(~smkinit_grs,data=dads))
dev.off()

pdf(file="./cigday_grs_all.pdf")
par(las=1,cex=1,mar=c(5,5,2,2))
plot(compareGroups(~cigday_grs,data=dat))
dev.off()

pdf(file="./cigday_grs_moms.pdf")
par(las=1,cex=1,mar=c(5,5,2,2))
plot(compareGroups(~cigday_grs,data=moms))
dev.off()

pdf(file="./cigday_grs_dads.pdf")
par(las=1,cex=1,mar=c(5,5,2,2))
plot(compareGroups(~cigday_grs,data=dads))
dev.off()

pdf(file="./agesmk_grs_all.pdf")
par(las=1,cex=1,mar=c(5,5,2,2))
plot(compareGroups(~agesmk_grs,data=dat))
dev.off()

pdf(file="./agesmk_grs_moms.pdf")
par(las=1,cex=1,mar=c(5,5,2,2))
plot(compareGroups(~agesmk_grs,data=moms))
dev.off()

pdf(file="./agesmk_grs_dads.pdf")
par(las=1,cex=1,mar=c(5,5,2,2))
plot(compareGroups(~agesmk_grs,data=dads))
dev.off()

pdf(file="./smkces_grs_all.pdf")
par(las=1,cex=1,mar=c(5,5,2,2))
plot(compareGroups(~smkces_grs,data=dat))
dev.off()

pdf(file="./smkces_grs_moms.pdf")
par(las=1,cex=1,mar=c(5,5,2,2))
plot(compareGroups(~smkces_grs,data=moms))
dev.off()

pdf(file="./smkces_grs_dads.pdf")
par(las=1,cex=1,mar=c(5,5,2,2))
plot(compareGroups(~smkces_grs,data=dads))
dev.off()

pdf(file="./lifesmk_grs_all.pdf")
par(las=1,cex=1,mar=c(5,5,2,2))
plot(compareGroups(~lifesmk_grs,data=dat))
dev.off()

pdf(file="./lifesmk_grs_moms.pdf")
par(las=1,cex=1,mar=c(5,5,2,2))
plot(compareGroups(~lifesmk_grs,data=moms))
dev.off()

pdf(file="./lifesmk_grs_dads.pdf")
par(las=1,cex=1,mar=c(5,5,2,2))
plot(compareGroups(~lifesmk_grs,data=dads))
dev.off()

pdf(file="./risktk_grs_all.pdf")
par(las=1,cex=1,mar=c(5,5,2,2))
plot(compareGroups(~risktk_grs,data=dat))
dev.off()

pdf(file="./risktk_grs_moms.pdf")
par(las=1,cex=1,mar=c(5,5,2,2))
plot(compareGroups(~risktk_grs,data=moms))
dev.off()

pdf(file="./risktk_grs_dads.pdf")
par(las=1,cex=1,mar=c(5,5,2,2))
plot(compareGroups(~risktk_grs,data=dads))
dev.off()

pdf(file="./bmi_grs_all.pdf")
par(las=1,cex=1,mar=c(5,5,2,2))
plot(compareGroups(~bmi_grs,data=dat))
dev.off()

pdf(file="./bmi_grs_moms.pdf")
par(las=1,cex=1,mar=c(5,5,2,2))
plot(compareGroups(~bmi_grs,data=moms))
dev.off()

pdf(file="./bmi_grs_dads.pdf")
par(las=1,cex=1,mar=c(5,5,2,2))
plot(compareGroups(~bmi_grs,data=dads))
dev.off()

pdf(file="./eduyears_grs_all.pdf")
par(las=1,cex=1,mar=c(5,5,2,2))
plot(compareGroups(~eduyears_grs,data=dat))
dev.off()

pdf(file="./eduyears_grs_moms.pdf")
par(las=1,cex=1,mar=c(5,5,2,2))
plot(compareGroups(~eduyears_grs,data=moms))
dev.off()

pdf(file="./eduyears_grs_dads.pdf")
par(las=1,cex=1,mar=c(5,5,2,2))
plot(compareGroups(~eduyears_grs,data=dads))
dev.off()



