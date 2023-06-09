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

library(devtools)
library(googleAuthR)
library(MendelianRandomization)
library(mr.raps)
library(meta)
library(MRPRESSO)
library(MRInstruments)
library(MRMix)
library(RadialMR)
library(ieugwasr)
library(TwoSampleMR)


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

beta_se_ic_guapa3 <- function(x, y) 
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
  hr<-round(exp(x),5)
  ic95a<-round(exp(x-(z*y)),5)
  ic95b<-round(exp(x+(z*y)),5)
  ic_ok<-ic_guapa(hr,ic95a,ic95b)
  return(ic_ok)
}

risk_se_ic_guapa3 <- function(x,y) 
{
  z<-qnorm(1-0.05/2)
  hr<-guapa(exp(x))
  ic95a<-guapa(exp(x-(z*y)))
  ic95b<-guapa(exp(x+(z*y)))
  ic_ok<-ic_guapa2(hr,ic95a,ic95b)
  return(ic_ok)
}


header.true <- function(df)
{
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}

z<-qnorm(1-0.05/2)


dir.create("C:/Users/ALHE/OneDrive - Folkehelseinstituttet/Articles/MR_smoking_COVID/Results")
dir.create("C:/Users/ALHE/OneDrive - Folkehelseinstituttet/Articles/MR_smoking_COVID/Results/Data")
setwd("C:/Users/ALHE/OneDrive - Folkehelseinstituttet/Articles/MR_smoking_COVID/Results")


################################
### CLEAN SUMMARY STATISTICS ###
################################

### EXPOSURES ###
#################

smkinit<-read.csv2("C:/Users/ALHE/OneDrive - Folkehelseinstituttet/Data/MR_GWAS_sources/MR_smoking/gwas_smkinit_liu2019.csv",header=TRUE,sep=";",dec=".")
smkinit<-rename.vars(smkinit,
                     from=c("rsid","effect_allele_freq","n"),
                     to=c("SNP","eaf","samplesize"))
smkinit<-smkinit[,c("SNP","effect_allele","other_allele","eaf","beta","se","pval","samplesize")]
smkinit$Phenotype<-c("Ever being a regular smoker")
smkinit$units<-c("units")
fwrite(smkinit,"./Data/liu2019_smkinit.txt",append=FALSE, sep="\t", sep2=c("\t","|","\t"), row.names=FALSE, col.names=TRUE)
smkinit<-NULL

cigday<-read.csv2("C:/Users/ALHE/OneDrive - Folkehelseinstituttet/Data/MR_GWAS_sources/MR_smoking/gwas_cigday_liu2019.csv",header=TRUE,sep=";",dec=".")
cigday<-rename.vars(cigday,
                     from=c("rsid","effect_allele_freq","n"),
                     to=c("SNP","eaf","samplesize"))
cigday<-cigday[,c("SNP","effect_allele","other_allele","eaf","beta","se","pval","samplesize")]
cigday$Phenotype<-c("Cigarettes per day, ever smokers")
cigday$units<-c("units")
fwrite(cigday,"./Data/liu2019_cigday.txt",append=FALSE, sep="\t", sep2=c("\t","|","\t"), row.names=FALSE, col.names=TRUE)
cigday<-NULL

agesmk<-read.csv2("C:/Users/ALHE/OneDrive - Folkehelseinstituttet/Data/MR_GWAS_sources/MR_smoking/gwas_agesmk_liu2019.csv",header=TRUE,sep=";",dec=".")
agesmk<-rename.vars(agesmk,
                     from=c("rsid","effect_allele_freq","n"),
                     to=c("SNP","eaf","samplesize"))
agesmk<-agesmk[,c("SNP","effect_allele","other_allele","eaf","beta","se","pval","samplesize")]
agesmk$Phenotype<-c("Age of initiation of regular smoking")
agesmk$units<-c("units")
fwrite(agesmk,"./Data/liu2019_agesmk.txt",append=FALSE, sep="\t", sep2=c("\t","|","\t"), row.names=FALSE, col.names=TRUE)
agesmk<-NULL

smkces<-read.csv2("C:/Users/ALHE/OneDrive - Folkehelseinstituttet/Data/MR_GWAS_sources/MR_smoking/gwas_smkces_liu2019.csv",header=TRUE,sep=";",dec=".")
smkces<-rename.vars(smkces,
                     from=c("rsid","effect_allele_freq","n"),
                     to=c("SNP","eaf","samplesize"))
smkces<-smkces[,c("SNP","effect_allele","other_allele","eaf","beta","se","pval","samplesize")]
smkces$Phenotype<-c("Smoking cessation")
smkces$units<-c("units")
fwrite(smkces,"./Data/liu2019_smkces.txt",append=FALSE, sep="\t", sep2=c("\t","|","\t"), row.names=FALSE, col.names=TRUE)
smkces<-NULL


### OUTCOMES ###
################

### COVID INFECTION ###

covid<-fread("C:/Users/ALHE/OneDrive - Folkehelseinstituttet/Data/GWAS_COVID19/hgi_covid_vs_population_20220403_EUR.tsv",header=TRUE,sep="\t",sep2="\t")
covid$SNP<-NULL
covid<-rename.vars(covid,
                   from=c("rsid","#CHR","POS","REF","ALT","all_meta_AF","all_inv_var_meta_beta","all_inv_var_meta_sebeta","all_inv_var_meta_p"),
                   to=c("SNP","chr","pos","other_allele","effect_allele","eaf","beta","se","pval"))
covid$samplesize<-as.numeric(covid$all_inv_var_meta_cases)+as.numeric(covid$all_inv_var_meta_controls)
covid$Units<-c("Units")
covid$Phenotype<-c("COVID19 infection")
covid<-covid[,c("SNP","chr","pos","effect_allele","other_allele","eaf","beta","se","pval","samplesize","Units","Phenotype")]
fwrite(covid,"./Data/hgi2022_covid_infection.txt",append=FALSE, sep="\t", sep2=c("\t","|","\t"), row.names=FALSE, col.names=TRUE)
covid<-NULL


############################
### DATA FORMAT FOR 2SMR ###
############################

memory.limit(35000)

### SmkInit ###

smkinit <- read_exposure_data(
  filename = "./Data/liu2019_smkinit.txt",
  clump = FALSE,
  sep="\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  units_col = "units",
  samplesize_col = "samplesize",
  id_col = "id",
  min_pval = 1e-200,
  log_pval = FALSE)

covid <- read_outcome_data(
  snps = smkinit$SNP,
  filename = "./Data/hgi2022_covid_infection.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize")

dat<-harmonise_data(smkinit, covid, action = 2)
save(dat,file="./Data/smkinit_covid_infection.RData")
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="./two_sample_MR/snps_smkinit_2smr.csv",sep=";",col.names=FALSE,row.names=FALSE)


### CigDay ###

cigday <- read_exposure_data(
  filename = "./Data/liu2019_cigday.txt",
  clump = FALSE,
  sep="\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  units_col = "units",
  samplesize_col = "samplesize",
  id_col = "id",
  min_pval = 1e-200,
  log_pval = FALSE)

covid <- read_outcome_data(
  snps = cigday$SNP,
  filename = "./Data/hgi2022_covid_infection.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize")

dat<-harmonise_data(cigday, covid, action = 2)
save(dat,file="./Data/cigday_covid_infection.RData")

xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="./two_sample_MR/snps_cigday_2smr.csv",sep=";",col.names=FALSE,row.names=FALSE)


### AgeSmk ###

agesmk <- read_exposure_data(
  filename = "./Data/liu2019_agesmk.txt",
  clump = FALSE,
  sep="\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  units_col = "units",
  samplesize_col = "samplesize",
  id_col = "id",
  min_pval = 1e-200,
  log_pval = FALSE)

covid <- read_outcome_data(
  snps = agesmk$SNP,
  filename = "./Data/hgi2022_covid_infection.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize")

dat<-harmonise_data(agesmk, covid, action = 2)
save(dat,file="./Data/agesmk_covid_infection.RData")

xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="./two_sample_MR/snps_agesmk_2smr.csv",sep=";",col.names=FALSE,row.names=FALSE)


### SmkCes ###

smkces <- read_exposure_data(
  filename = "./Data/liu2019_smkces.txt",
  clump = FALSE,
  sep="\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  units_col = "units",
  samplesize_col = "samplesize",
  id_col = "id",
  min_pval = 1e-200,
  log_pval = FALSE)

covid <- read_outcome_data(
  snps = smkces$SNP,
  filename = "./Data/hgi2022_covid_infection.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize")

dat<-harmonise_data(smkces, covid, action = 2)
save(dat,file="./Data/smkces_covid_infection.RData")

xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="./two_sample_MR/snps_smkces_2smr.csv",sep=";",col.names=FALSE,row.names=FALSE)


##############################
### TWO-SAMPLE MR ANALYSES ###
##############################

setwd("C:/Users/ALHE/OneDrive - Folkehelseinstituttet/Articles/MR_smoking_COVID/Results/")

z<-qnorm(1-0.05/2)
Isq = function(y,s){
  k          = length(y)
  w          = 1/s^2; sum.w  = sum(w)
  mu.hat     = sum(y*w)/sum.w  
  Q          = sum(w*(y-mu.hat)^2)
  Isq        = (Q - (k-1))/Q
  Isq        = max(0,Isq)
  return(Isq)
}

# F>10: good IVW performance
# Isq close to 1: good MR-Egger performance

vars01<-c("smkinit_",
          "cigday_",
          "agesmk_",
          "smkces_")
vars02<-c("covid_infection"
          "covid_infection"
          "covid_infection"
          "covid_infection")

tab<-NULL
for(i in 1:length(vars01))
  
{
  namedat<-paste("C:/Users/ALHE/OneDrive - Folkehelseinstituttet/Articles/MR_smoking_COVID/Results/Data/",vars01[i],vars02[i],".RData",sep="")
  load(namedat)
  
  mr_results<-mr(dat,method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))
  mr_results$pval<-pval_guapa(mr_results$pval)
  mr_results$beta<-beta_se_ic_guapa3(mr_results$b,mr_results$se)
  mr_results$or<-risk_se_ic_guapa3(mr_results$b,mr_results$se)
  mr_ivw_mult<-mr_ivw_mre(b_exp=dat$beta.exposure,b_out=dat$beta.outcome,se_exp=dat$se.exposure,se_out=dat$se.outcome)
  mr_ivw_mult_beta<-beta_se_ic_guapa3(mr_ivw_mult$b,mr_ivw_mult$se)
  mr_ivw_mult_or<-risk_se_ic_guapa3(mr_ivw_mult$b,mr_ivw_mult$se)
  mr_ivw_mult_pval<-pval_guapa(mr_ivw_mult$pval)
  ivw_mult_b<-round(mr_ivw_mult$b,6)
  ivw_mult_b_lo<-round(mr_ivw_mult$b-(z*mr_ivw_mult$se),6)
  ivw_mult_b_hi<-round(mr_ivw_mult$b+(z*mr_ivw_mult$se),6)
  ivw_mult_or<-round(exp(mr_ivw_mult$b),6)
  ivw_mult_or_lo<-round(exp(mr_ivw_mult$b-(z*mr_ivw_mult$se)),6)
  ivw_mult_or_hi<-round(exp(mr_ivw_mult$b+(z*mr_ivw_mult$se)),6)
  mr_raps<-mr.raps(b_exp=dat$beta.exposure,b_out=dat$beta.outcome,se_exp=dat$se.exposure,se_out=dat$se.outcome,diagnosis=FALSE)
  mr_raps_beta<-beta_se_ic_guapa3(mr_raps$beta.hat,mr_raps$beta.se)
  mr_raps_or<-risk_se_ic_guapa3(mr_raps$beta.hat,mr_raps$beta.se)
  mr_raps_pval<-pval_guapa(mr_raps$beta.p.value)
  #mr_presso_res<-run_mr_presso(dat,NbDistribution = 1000)[[1]]$`Main MR results`
  #mr_presso_raw_beta<-beta_se_ic_guapa(mr_presso_res$`Causal Estimate`[1],mr_presso_res$Sd[1])
  #mr_presso_raw_or<-risk_se_ic_guapa(mr_presso_res$`Causal Estimate`[1],mr_presso_res$Sd[1])
  #mr_presso_raw_pval<-pval_guapa(mr_presso_res$`P-value`[1])
  #mr_presso_corr_beta<-beta_se_ic_guapa(mr_presso_res$`Causal Estimate`[2],mr_presso_res$Sd[2])
  #mr_presso_corr_or<-risk_se_ic_guapa(mr_presso_res$`Causal Estimate`[2],mr_presso_res$Sd[2])
  #mr_presso_corr_pval<-pval_guapa(mr_presso_res$`P-value`[2])
  q_cochran<-paste(guapa(mr_heterogeneity(dat)$Q[2])," (P=",pval_guapa(mr_heterogeneity(dat)$Q_pval[2]),")",sep="")
  q_rucker<-paste(guapa(mr_heterogeneity(dat)$Q[1])," (P=",pval_guapa(mr_heterogeneity(dat)$Q_pval[1]),")",sep="")
  q_cochran_ivw_mult<-paste(guapa(mr_ivw_mult$Q)," (P=",pval_guapa(mr_ivw_mult$Q_pval),")",sep="")
  p_leaveoneout<-pval_guapa(max(head(mr_leaveoneout(dat)$p, -1),na.rm=TRUE))
  plei<-pval_guapa(mr_pleiotropy_test(dat)$pval)
  fstat<-guapa(mean((abs(dat$beta.exposure))^2/dat$se.exposure^2,na.rm=TRUE))
  unIsq<-guapa(Isq(abs(dat$beta.exposure),dat$se.exposure))
  
  tab<-rbind(tab,cbind(mr_results$nsnp[1],
                       mr_results$beta[1],mr_results$or[1],mr_results$pval[1],
                       mr_ivw_mult_beta,mr_ivw_mult_or,mr_ivw_mult_pval,
                       mr_results$beta[2],mr_results$or[2],mr_results$pval[2],
                       mr_results$beta[3],mr_results$or[3],mr_results$pval[3],
                       mr_results$beta[4],mr_results$or[4],mr_results$pval[4],
                       #mr_presso_raw_beta,mr_presso_raw_or,mr_presso_raw_pval,mr_presso_corr_beta,mr_presso_corr_or,mr_presso_corr_pval,
                       mr_raps_beta,mr_raps_or,mr_raps_pval,
                       plei,q_cochran,q_rucker,q_cochran_ivw_mult,p_leaveoneout,
                       ivw_mult_b,ivw_mult_b_lo,ivw_mult_b_hi,ivw_mult_or,ivw_mult_or_lo,ivw_mult_or_hi,fstat,unIsq))
}

colnames(tab)<-c("SNPs","IVW_random_b","IVW_random_or","IVW_random_p","IVW_mult_b","IVW_mult_or","IVW_mult_p",
                 "Egger_b","Egger_or","Egger_p","WMe_b","WMe_or","Wme_p","WMo_b","WMo_or","WMo_p","RAPS_b","RAPS_or","RAPS_p",
                 #"PRESSO_raw_b","PRESSO_raw_or","PRESSO_raw_p","PRESSO_corr_b","PRESSO_corr_or","PRESSO_corr_p",
                 "Egger_pleio","Cochran_Q","Rucker_Q","Cochran_Q_IVW_mult","p_max_leaveoneout",
                 "ivw_b","ivw_b_lo","ivw_b_hi","ivw_or","ivw_or_lo","ivw_or_hi","F-stat","I2")
rownames(tab)<-paste(vars01,vars02,sep="")
write.table(tab,file="./two_sample_MR/results_2smr.csv",sep=";",col.names=NA)


