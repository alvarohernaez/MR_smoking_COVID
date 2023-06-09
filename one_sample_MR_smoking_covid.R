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


#####################
### MAIN ANALYSES ###
#####################

setwd("N:/durable/Projects/ALHE_smoking_covid/")


### POPULATION DESCRIPTION ###
##############################

load("./Data/mr_smoking_covid_all.RData")
dat$obesity<-with(dat,ifelse(bmi<25,0,
                             ifelse(bmi<30,1,
                                    ifelse(bmi>=30,2,NA))))
#plot(compareGroups(~cohabs,data=dat))

xxx<-dat[,c("role","age","edu","work","obesity","cohabs","region","smkinit","cigday","smkces","agesmk","covid")]
xxx$sel<-1
xmom<-xxx[xxx$role==0,]
xdad<-xxx[xxx$role==1,]

all01<-createTable(compareGroups(sel~.
                                 -covid,
                                 xxx, method=c("role"=3,"edu"=3,"work"=3,"obesity"=3,"region"=3,
                                               "smkinit"=3,"cigday"=3,"smkces"=3)),
                   show.n=TRUE, show.p.overall=FALSE, show.p.trend=FALSE)
comp01<-createTable(compareGroups(covid~.
                                  -sel,
                                  xxx, method=c("role"=3,"edu"=3,"work"=3,"obesity"=3,"region"=3,
                                                "smkinit"=3,"cigday"=3,"smkces"=3)),
                    show.n=TRUE, show.p.overall=TRUE, show.p.trend=FALSE)
all02<-createTable(compareGroups(sel~.
                                 -covid,
                                 xmom, method=c("role"=3,"edu"=3,"work"=3,"obesity"=3,"region"=3,
                                               "smkinit"=3,"cigday"=3,"smkces"=3)),
                   show.n=TRUE, show.p.overall=FALSE, show.p.trend=FALSE)
comp02<-createTable(compareGroups(covid~.
                                  -sel,
                                  xmom, method=c("role"=3,"edu"=3,"work"=3,"obesity"=3,"region"=3,
                                                "smkinit"=3,"cigday"=3,"smkces"=3)),
                    show.n=TRUE, show.p.overall=TRUE, show.p.trend=FALSE)
all03<-createTable(compareGroups(sel~.
                                 -covid,
                                 xdad, method=c("role"=3,"edu"=3,"work"=3,"obesity"=3,"region"=3,
                                               "smkinit"=3,"cigday"=3,"smkces"=3)),
                   show.n=TRUE, show.p.overall=FALSE, show.p.trend=FALSE)
comp03<-createTable(compareGroups(covid~.
                                  -sel,
                                  xdad, method=c("role"=3,"edu"=3,"work"=3,"obesity"=3,"region"=3,
                                                "smkinit"=3,"cigday"=3,"smkces"=3)),
                    show.n=TRUE, show.p.overall=TRUE, show.p.trend=FALSE)

tab<-NULL
tab<-as.data.frame(cbind(all01$descr[,1],comp01$descr,
                         c("",all02$descr[,1]),rbind("",comp02$descr),
                         c("",all03$descr[,1]),rbind("",comp03$descr)))
colnames(tab)<-c("All","All-No COVID","All-Yes COVID","All-Pvalue","All-N",
                 "Women","Women-No COVID","Women-Yes COVID","Women-Pvalue","Women-N",
                 "Men","Men-No COVID","Men-Yes COVID","Men-Pvalue","Men-N")
write.table(tab,file="./Outputs/Descriptive/population_description.csv",sep=";",col.names=NA)


### ASSOCIATION BETWEEN SMOKING TRAITS AND GRSs ###
###################################################

load("./Data/mr_smoking_covid_all.RData")
library(pROC)

# BINARY VARIABLES #

vars01<-c("smkinit","smkces")
vars02<-c("smkinit_grs","smkces_grs")
vars03<-c("351","22")
vars04<-c("all","women","men")

tab<-NULL

for(i in 1:length(vars01))
  
{
  xxx<-as.numeric(with(dat,scale(dat[,vars02[i]])))
  mod01<-glm(formula=as.factor(dat[,vars01[i]])~xxx, data=dat, family="binomial")
  coef_a<-risk_se_ic_guapa3(as.numeric(summary(mod01)$coefficients[2,1]),as.numeric(summary(mod01)$coefficients[2,2]))
  pval_a<-pval_guapa(as.numeric(summary(mod01)$coefficients[2,4]))
  llrtest_a<-round(as.numeric(anova(mod01,test=c("F"))[2,5]),0) # log-likelihood ratio test (~F in binomial GLM)
  pseudor2_a<-paste("'",guapa(NagelkerkeR2(mod01)$R2*100),sep="") # pseudoR2 by the Nagelkerke method
  n_snps_a<-vars03[i]
  dat$prob<-predict.glm(update(mod01,na.action=na.exclude),type=c("response"))
  roc_auc_a<-guapa(roc(as.factor(dat[,vars01[i]])~prob,data=dat)$auc)
  mean_grs_a<-paste(guapa(mean(dat[,vars02[i]],na.rm=TRUE))," (",guapa(sd(dat[,vars02[i]],na.rm=TRUE)),")",sep="")
  tab<-rbind(tab,cbind(n_snps_a,mean_grs_a,coef_a,pval_a,llrtest_a,roc_auc_a,pseudor2_a))
  
  datx<-dat[dat$role==0,]
  xxx<-as.numeric(with(datx,scale(datx[,vars02[i]])))
  mod01<-glm(formula=as.factor(datx[,vars01[i]])~xxx, data=datx, family="binomial")
  coef_w<-risk_se_ic_guapa3(as.numeric(summary(mod01)$coefficients[2,1]),as.numeric(summary(mod01)$coefficients[2,2]))
  pval_w<-pval_guapa(as.numeric(summary(mod01)$coefficients[2,4]))
  llrtest_w<-round(as.numeric(anova(mod01,test=c("F"))[2,5]),0) # log-likelihood ratio test (~F in binomial GLM)
  pseudor2_w<-paste("'",guapa(NagelkerkeR2(mod01)$R2*100),sep="") # pseudoR2 by the Nagelkerke method
  n_snps_w<-vars03[i]
  datx$prob<-predict.glm(update(mod01,na.action=na.exclude),type=c("response"))
  roc_auc_w<-guapa(roc(as.factor(datx[,vars01[i]])~prob,data=datx)$auc)
  mean_grs_w<-paste(guapa(mean(datx[,vars02[i]],na.rm=TRUE))," (",guapa(sd(datx[,vars02[i]],na.rm=TRUE)),")",sep="")
  tab<-rbind(tab,cbind(n_snps_w,mean_grs_w,coef_w,pval_w,llrtest_w,roc_auc_w,pseudor2_w))
  
  datx<-dat[dat$role==1,]
  xxx<-as.numeric(with(datx,scale(datx[,vars02[i]])))
  mod01<-glm(formula=as.factor(datx[,vars01[i]])~xxx, data=datx, family="binomial")
  coef_m<-risk_se_ic_guapa3(as.numeric(summary(mod01)$coefficients[2,1]),as.numeric(summary(mod01)$coefficients[2,2]))
  pval_m<-pval_guapa(as.numeric(summary(mod01)$coefficients[2,4]))
  llrtest_m<-round(as.numeric(anova(mod01,test=c("F"))[2,5]),0) # log-likelihood ratio test (~F in binomial GLM)
  pseudor2_m<-paste("'",guapa(NagelkerkeR2(mod01)$R2*100),sep="") # pseudoR2 by the Nagelkerke method
  n_snps_m<-vars03[i]
  datx$prob<-predict.glm(update(mod01,na.action=na.exclude),type=c("response"))
  roc_auc_m<-guapa(roc(as.factor(datx[,vars01[i]])~prob,data=datx)$auc)
  mean_grs_m<-paste(guapa(mean(datx[,vars02[i]],na.rm=TRUE))," (",guapa(sd(datx[,vars02[i]],na.rm=TRUE)),")",sep="")
  tab<-rbind(tab,cbind(n_snps_m,mean_grs_m,coef_m,pval_m,llrtest_m,roc_auc_m,pseudor2_m))
}

rownames(tab)<-apply(expand.grid(vars04,vars01), 1, paste, collapse="_")
write.table(tab,file="./Outputs/Descriptive/robustness_binary.csv",sep=";",col.names=NA)


# CONTINUOUS VARIABLES #

dat$cigday2<-with(dat,ifelse(cigday=="0",NA,as.numeric(as.character(cigday))))
dat$cigday2<-with(dat,ifelse(smkinit==0,NA,cigday2))
dat$agesmk2<-with(dat,ifelse(smkinit==0,NA,agesmk))

vars01<-c("cigday2","agesmk2")
vars02<-c("cigday_grs","agesmk_grs")
vars03<-c("53","9")
vars04<-c("all","women","men")

tab<-NULL

for(i in 1:length(vars01))
  
{
  xxx<-as.numeric(with(dat,scale(dat[,vars02[i]])))
  mod01<-lm(dat[,vars01[i]]~xxx, data=dat)
  coef_a<-beta_se_ic_guapa2(as.numeric(summary(mod01)$coefficients[2,1]),as.numeric(summary(mod01)$coefficients[2,2]))
  pval_a<-pval_guapa(as.numeric(summary(mod01)$coefficients[2,4]))
  fstat_a<-paste("'",guapa(as.numeric(summary(mod01)$fstat[1])),sep="")
  r2_a<-paste("'",guapa(summary(mod01)$adj.r.squared*100),sep="")
  n_snps_a<-vars03[i]
  mean_grs_a<-paste(guapa(mean(dat[,vars02[i]],na.rm=TRUE))," (",guapa(sd(dat[,vars02[i]],na.rm=TRUE)),")",sep="")
  tab<-rbind(tab,cbind(n_snps_a,mean_grs_a,coef_a,pval_a,fstat_a,r2_a))
  
  datx<-dat[dat$role==0,]
  xxx<-as.numeric(with(datx,scale(datx[,vars02[i]])))
  mod01<-lm(datx[,vars01[i]]~xxx, data=datx)
  coef_w<-beta_se_ic_guapa2(as.numeric(summary(mod01)$coefficients[2,1]),as.numeric(summary(mod01)$coefficients[2,2]))
  pval_w<-pval_guapa(as.numeric(summary(mod01)$coefficients[2,4]))
  fstat_w<-paste("'",guapa(as.numeric(summary(mod01)$fstat[1])),sep="")
  r2_w<-paste("'",guapa(summary(mod01)$adj.r.squared*100),sep="")
  n_snps_w<-vars03[i]
  mean_grs_w<-paste(guapa(mean(datx[,vars02[i]],na.rm=TRUE))," (",guapa(sd(datx[,vars02[i]],na.rm=TRUE)),")",sep="")
  tab<-rbind(tab,cbind(n_snps_w,mean_grs_w,coef_w,pval_w,fstat_w,r2_w))
  
  datx<-dat[dat$role==1,]
  xxx<-as.numeric(with(datx,scale(datx[,vars02[i]])))
  mod01<-lm(datx[,vars01[i]]~xxx, data=datx)
  coef_m<-beta_se_ic_guapa2(as.numeric(summary(mod01)$coefficients[2,1]),as.numeric(summary(mod01)$coefficients[2,2]))
  pval_m<-pval_guapa(as.numeric(summary(mod01)$coefficients[2,4]))
  fstat_m<-paste("'",guapa(as.numeric(summary(mod01)$fstat[1])),sep="")
  r2_m<-paste("'",guapa(summary(mod01)$adj.r.squared*100),sep="")
  n_snps_m<-vars03[i]
  mean_grs_m<-paste(guapa(mean(datx[,vars02[i]],na.rm=TRUE))," (",guapa(sd(datx[,vars02[i]],na.rm=TRUE)),")",sep="")
  tab<-rbind(tab,cbind(n_snps_m,mean_grs_m,coef_m,pval_m,fstat_m,r2_m))
}

rownames(tab)<-head(apply(expand.grid(vars04,vars01), 1, paste, collapse="_"),-1)
write.table(tab,file="./Outputs/Descriptive/robustness_continuous.csv",sep=";",col.names=NA)


### ASSOCIATION OF GRS WITH COVARIATES ###
##########################################

vars01<-c("smkinit_grs","cigday_grs","smkces_grs","agesmk_grs")

tab<-NULL

for(i in 1:length(vars01))
  
{
  mod01<-lm(age~dat[,vars01[i]], data=dat)
  mod02<-lm(eduyears~dat[,vars01[i]], data=dat)
  mod03<-lm(bmi~dat[,vars01[i]], data=dat)
  
  coef_age<-beta_se_ic_guapa2(as.numeric(summary(mod01)$coefficients[2,1]),as.numeric(summary(mod01)$coefficients[2,2]))
  pval_age<-pval_guapa(as.numeric(summary(mod01)$coefficients[2,4]))
  coef_edu<-beta_se_ic_guapa2(as.numeric(summary(mod02)$coefficients[2,1]),as.numeric(summary(mod02)$coefficients[2,2]))
  pval_edu<-pval_guapa(as.numeric(summary(mod02)$coefficients[2,4]))
  coef_bmi<-beta_se_ic_guapa2(as.numeric(summary(mod03)$coefficients[2,1]),as.numeric(summary(mod03)$coefficients[2,2]))
  pval_bmi<-pval_guapa(as.numeric(summary(mod03)$coefficients[2,4]))
  
  tab<-rbind(tab,cbind(coef_age,pval_age,coef_edu,pval_edu,coef_bmi,pval_bmi))
}

rownames(tab)<-vars01
write.table(tab,file="./Outputs/Descriptive/pleiotropy.csv",sep=";",col.names=NA)


##########################################################################
### LOGISTIC REGRESSION / MENDELIAN RANDOMIZATION: LINEAR ASSOCIATIONS ###
##########################################################################

### CALCULATION OF GENETICALLY-DETERMINED VARIABLES ###
#######################################################

load("./Data/mr_smoking_covid_all.RData")

# Continuous variables #

dat$cigday2<-with(dat,ifelse(cigday=="0",NA,as.numeric(as.character(cigday))))
dat$cigday2<-with(dat,ifelse(smkinit==0,NA,cigday2))
dat$agesmk2<-with(dat,ifelse(smkinit==0,NA,agesmk))

vars01<-c("cigday2","agesmk2")
vars02<-c("cigday_grs","agesmk_grs")
vars03<-c("cigday_gp","agesmk_gp")
vars04<-c("cigday2_z","agesmk2_z")
vars05<-c("cigday_grs_z","agesmk_grs_z")
vars06<-c("cigday_gp_z","agesmk_gp_z")

for(i in 1:length(vars01))
  
{
  mod01<-lm(dat[,vars01[i]]~dat[,vars02[i]], data=dat)
  dat[,vars03[i]]<-predict.lm(update(mod01,na.action=na.exclude), type="response")
  dat[,vars04[i]]<-as.numeric(with(dat,scale(dat[,vars01[i]])))
  dat[,vars05[i]]<-as.numeric(with(dat,scale(dat[,vars02[i]])))
  dat[,vars06[i]]<-as.numeric(with(dat,scale(dat[,vars03[i]])))
}


# Binary outcomes (genetically-predicted likelihood) #

dat$smkces<-with(dat,ifelse(smkinit=="0",NA,smkces))

vars01<-c("smkinit","smkces")
vars02<-c("smkinit_grs","smkces_grs")
vars03<-c("smkinit_gp","smkces_gp")
vars05<-c("smkinit_grs_z","smkces_grs_z")
vars06<-c("smkinit_gp_z","smkces_gp_z")

for(i in 1:length(vars01))
  
{
  mod01<-glm(dat[,vars01[i]]~dat[,vars02[i]], data=dat, family="binomial")
  dat[,vars03[i]]<-predict.glm(update(mod01,na.action=na.exclude), type="response")
  dat[,vars05[i]]<-as.numeric(with(dat,scale(dat[,vars02[i]])))
  dat[,vars06[i]]<-as.numeric(with(dat,scale(dat[,vars03[i]])))
}


### ANALYSES ###
################

vars01<-c("smkinit","cigday2_z","smkces","agesmk2_z")
vars02<-c("smkinit_gp_z","cigday_gp_z","smkces_gp_z","agesmk_gp_z")
vars03<-c("all","women","men")

tab<-NULL
for(i in 1:length(vars01))
  
{
  datx<-dat
  mod01<-glm(formula=as.factor(covid)~datx[,vars01[i]],
             data=datx, family="binomial")
  coef01<-risk_se_ic_guapa2(as.numeric(summary(mod01)$coefficients[2,1]),as.numeric(summary(mod01)$coefficients[2,2]))
  pval01<-pval_guapa(as.numeric(summary(mod01)$coefficients[2,4]))
  or01<-exp(as.numeric(summary(mod01)$coefficients[2,1]))
  orlo01<-exp(as.numeric(summary(mod01)$coefficients[2,1])-(z*as.numeric(summary(mod01)$coefficients[2,2])))
  orhi01<-exp(as.numeric(summary(mod01)$coefficients[2,1])+(z*as.numeric(summary(mod01)$coefficients[2,2])))
  
  mod02<-glm(formula=as.factor(covid)~datx[,vars01[i]]+age+role+eduyears+bmi+cohabs+as.factor(region)+as.factor(work),
             data=datx, family="binomial")
  coef02<-risk_se_ic_guapa2(as.numeric(summary(mod02)$coefficients[2,1]),as.numeric(summary(mod02)$coefficients[2,2]))
  pval02<-pval_guapa(as.numeric(summary(mod02)$coefficients[2,4]))
  or02<-exp(as.numeric(summary(mod02)$coefficients[2,1]))
  orlo02<-exp(as.numeric(summary(mod02)$coefficients[2,1])-(z*as.numeric(summary(mod02)$coefficients[2,2])))
  orhi02<-exp(as.numeric(summary(mod02)$coefficients[2,1])+(z*as.numeric(summary(mod02)$coefficients[2,2])))
  
  mod03<-glm(formula=as.factor(covid)~datx[,vars02[i]]
             +pc01+pc02+pc03+pc04+pc05+pc06+pc07+pc08+pc09+pc10
             +pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+as.factor(genotype_batch),
             data=datx, family="binomial")
  coef03<-risk_se_ic_guapa2(as.numeric(summary(mod03)$coefficients[2,1]),as.numeric(summary(mod03)$coefficients[2,2]))
  pval03<-pval_guapa(as.numeric(summary(mod03)$coefficients[2,4]))
  or03<-exp(as.numeric(summary(mod03)$coefficients[2,1]))
  orlo03<-exp(as.numeric(summary(mod03)$coefficients[2,1])-(z*as.numeric(summary(mod03)$coefficients[2,2])))
  orhi03<-exp(as.numeric(summary(mod03)$coefficients[2,1])+(z*as.numeric(summary(mod03)$coefficients[2,2])))
  
  tab<-rbind(tab,cbind(coef01,pval01,or01,orlo01,orhi01,coef02,pval02,or02,orlo02,orhi02,coef03,pval03,or03,orlo03,orhi03))
  
  datx<-dat[dat$role==0,]
  mod04<-glm(formula=as.factor(covid)~datx[,vars01[i]],
             data=datx, family="binomial")
  coef04<-risk_se_ic_guapa2(as.numeric(summary(mod04)$coefficients[2,1]),as.numeric(summary(mod04)$coefficients[2,2]))
  pval04<-pval_guapa(as.numeric(summary(mod04)$coefficients[2,4]))
  or04<-exp(as.numeric(summary(mod04)$coefficients[2,1]))
  orlo04<-exp(as.numeric(summary(mod04)$coefficients[2,1])-(z*as.numeric(summary(mod04)$coefficients[2,2])))
  orhi04<-exp(as.numeric(summary(mod04)$coefficients[2,1])+(z*as.numeric(summary(mod04)$coefficients[2,2])))
  
  mod05<-glm(formula=as.factor(covid)~datx[,vars01[i]]+age+eduyears+bmi+cohabs+as.factor(region)+as.factor(work),
             data=datx, family="binomial")
  coef05<-risk_se_ic_guapa2(as.numeric(summary(mod05)$coefficients[2,1]),as.numeric(summary(mod05)$coefficients[2,2]))
  pval05<-pval_guapa(as.numeric(summary(mod05)$coefficients[2,4]))
  or05<-exp(as.numeric(summary(mod05)$coefficients[2,1]))
  orlo05<-exp(as.numeric(summary(mod05)$coefficients[2,1])-(z*as.numeric(summary(mod05)$coefficients[2,2])))
  orhi05<-exp(as.numeric(summary(mod05)$coefficients[2,1])+(z*as.numeric(summary(mod05)$coefficients[2,2])))
  
  mod06<-glm(formula=as.factor(covid)~datx[,vars02[i]]
             +pc01+pc02+pc03+pc04+pc05+pc06+pc07+pc08+pc09+pc10
             +pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+as.factor(genotype_batch),
             data=datx, family="binomial")
  coef06<-risk_se_ic_guapa2(as.numeric(summary(mod06)$coefficients[2,1]),as.numeric(summary(mod06)$coefficients[2,2]))
  pval06<-pval_guapa(as.numeric(summary(mod06)$coefficients[2,4]))
  or06<-exp(as.numeric(summary(mod06)$coefficients[2,1]))
  orlo06<-exp(as.numeric(summary(mod06)$coefficients[2,1])-(z*as.numeric(summary(mod06)$coefficients[2,2])))
  orhi06<-exp(as.numeric(summary(mod06)$coefficients[2,1])+(z*as.numeric(summary(mod06)$coefficients[2,2])))
  
  tab<-rbind(tab,cbind(coef04,pval04,or04,orlo04,orhi04,coef05,pval05,or05,orlo05,orhi05,coef06,pval06,or06,orlo06,orhi06))
  
  datx<-dat[dat$role==1,]
  mod07<-glm(formula=as.factor(covid)~datx[,vars01[i]],
             data=datx, family="binomial")
  coef07<-risk_se_ic_guapa2(as.numeric(summary(mod07)$coefficients[2,1]),as.numeric(summary(mod07)$coefficients[2,2]))
  pval07<-pval_guapa(as.numeric(summary(mod07)$coefficients[2,4]))
  or07<-exp(as.numeric(summary(mod07)$coefficients[2,1]))
  orlo07<-exp(as.numeric(summary(mod07)$coefficients[2,1])-(z*as.numeric(summary(mod07)$coefficients[2,2])))
  orhi07<-exp(as.numeric(summary(mod07)$coefficients[2,1])+(z*as.numeric(summary(mod07)$coefficients[2,2])))
  
  mod08<-glm(formula=as.factor(covid)~datx[,vars01[i]]+age+eduyears+bmi+cohabs+as.factor(region)+as.factor(work),
             data=datx, family="binomial")
  coef08<-risk_se_ic_guapa2(as.numeric(summary(mod08)$coefficients[2,1]),as.numeric(summary(mod08)$coefficients[2,2]))
  pval08<-pval_guapa(as.numeric(summary(mod08)$coefficients[2,4]))
  or08<-exp(as.numeric(summary(mod08)$coefficients[2,1]))
  orlo08<-exp(as.numeric(summary(mod08)$coefficients[2,1])-(z*as.numeric(summary(mod08)$coefficients[2,2])))
  orhi08<-exp(as.numeric(summary(mod08)$coefficients[2,1])+(z*as.numeric(summary(mod08)$coefficients[2,2])))
  
  mod09<-glm(formula=as.factor(covid)~datx[,vars02[i]]
             +pc01+pc02+pc03+pc04+pc05+pc06+pc07+pc08+pc09+pc10
             +pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+as.factor(genotype_batch),
             data=datx, family="binomial")
  coef09<-risk_se_ic_guapa2(as.numeric(summary(mod09)$coefficients[2,1]),as.numeric(summary(mod09)$coefficients[2,2]))
  pval09<-pval_guapa(as.numeric(summary(mod09)$coefficients[2,4]))
  or09<-exp(as.numeric(summary(mod09)$coefficients[2,1]))
  orlo09<-exp(as.numeric(summary(mod09)$coefficients[2,1])-(z*as.numeric(summary(mod09)$coefficients[2,2])))
  orhi09<-exp(as.numeric(summary(mod09)$coefficients[2,1])+(z*as.numeric(summary(mod09)$coefficients[2,2])))
  
  tab<-rbind(tab,cbind(coef07,pval07,or07,orlo07,orhi07,coef08,pval08,or08,orlo08,orhi08,coef09,pval09,or09,orlo09,orhi09))
  
}

colnames(tab)<-c("Raw-Coef","Raw-pval","Raw-OR","Raw-lo","Raw-hi",
                 "MV-Coef","MV-pval","MV-OR","MV-lo","MV-hi",
                 "MR-Coef","MR-pval","MR-OR","MR-lo","MR-hi")
rownames(tab)<-head(apply(expand.grid(vars03,vars01), 1, paste, collapse="_"),-1)
write.table(tab,file="./Outputs/Results/main.csv",sep=";",col.names=NA)


### NO RELEVANCE ANALYSES ###
#############################

vars01<-c("cigday_grs","smkces_grs","agesmk_grs")
vars02<-c("cigday2","smkces","agesmk2")
vars03<-c("all","women","men")

tab<-NULL
for(i in 1:length(vars01))
  
{
  dat2<-dat[is.na(dat[,vars02[i]]),]
  datx<-dat2
  xxx<-as.numeric(with(datx,scale(datx[,vars01[i]])))
  mod03<-glm(formula=as.factor(covid)~xxx
             +pc01+pc02+pc03+pc04+pc05+pc06+pc07+pc08+pc09+pc10
             +pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+as.factor(genotype_batch),
             data=datx, family="binomial")
  coef03<-risk_se_ic_guapa2(as.numeric(summary(mod03)$coefficients[2,1]),as.numeric(summary(mod03)$coefficients[2,2]))
  pval03<-pval_guapa(as.numeric(summary(mod03)$coefficients[2,4]))

  tab<-rbind(tab,cbind(coef03,pval03))
  
  datx<-dat2[dat2$role==0,]
  xxx<-as.numeric(with(datx,scale(datx[,vars01[i]])))
  mod06<-glm(formula=as.factor(covid)~xxx
             +pc01+pc02+pc03+pc04+pc05+pc06+pc07+pc08+pc09+pc10
             +pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+as.factor(genotype_batch),
             data=datx, family="binomial")
  coef06<-risk_se_ic_guapa2(as.numeric(summary(mod06)$coefficients[2,1]),as.numeric(summary(mod06)$coefficients[2,2]))
  pval06<-pval_guapa(as.numeric(summary(mod06)$coefficients[2,4]))

  tab<-rbind(tab,cbind(coef06,pval06))
  
  datx<-dat2[dat2$role==1,]
  xxx<-as.numeric(with(datx,scale(datx[,vars01[i]])))
  mod09<-glm(formula=as.factor(covid)~xxx
             +pc01+pc02+pc03+pc04+pc05+pc06+pc07+pc08+pc09+pc10
             +pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+as.factor(genotype_batch),
             data=datx, family="binomial")
  coef09<-risk_se_ic_guapa2(as.numeric(summary(mod09)$coefficients[2,1]),as.numeric(summary(mod09)$coefficients[2,2]))
  pval09<-pval_guapa(as.numeric(summary(mod09)$coefficients[2,4]))

  tab<-rbind(tab,cbind(coef09,pval09))
}

colnames(tab)<-c("MR-Coef","MR-pval")
rownames(tab)<-apply(expand.grid(vars03,vars01), 1, paste, collapse="_")
write.table(tab,file="./Outputs/Results/no_relevance.csv",sep=";",col.names=NA)


##################################################################
### MULTIVARIABLE MENDELIAN RANDOMIZATION: LINEAR ASSOCIATIONS ###
##################################################################

### CALCULATION OF GENETICALLY-DETERMINED VARIABLES ###
#######################################################

load("./Data/mr_smoking_covid_all.RData")

# Continuous variables #

dat$cigday2<-with(dat,ifelse(cigday=="0",NA,as.numeric(as.character(cigday))))
dat$cigday2<-with(dat,ifelse(smkinit==0,NA,cigday2))
dat$agesmk2<-with(dat,ifelse(smkinit==0,NA,agesmk))

vars01<-c("cigday2","agesmk2")
vars02<-c("cigday_grs","agesmk_grs")
vars03<-c("cigday_gp","agesmk_gp")
vars04<-c("cigday2_z","agesmk2_z")
vars05<-c("cigday_grs_z","agesmk_grs_z")
vars06<-c("cigday_gp_z","agesmk_gp_z")

for(i in 1:length(vars01))
  
{
  mod01<-lm(dat[,vars01[i]]~dat[,vars02[i]]+bmi_grs+eduyears_grs+risktk_grs, data=dat)
  dat[,vars03[i]]<-predict.lm(update(mod01,na.action=na.exclude), type="response")
  dat[,vars04[i]]<-as.numeric(with(dat,scale(dat[,vars01[i]])))
  dat[,vars05[i]]<-as.numeric(with(dat,scale(dat[,vars02[i]])))
  dat[,vars06[i]]<-as.numeric(with(dat,scale(dat[,vars03[i]])))
}


# Binary outcomes (genetically-predicted likelihood) #

dat$smkces<-with(dat,ifelse(smkinit=="0",NA,smkces))

vars01<-c("smkinit","smkces")
vars02<-c("smkinit_grs","smkces_grs")
vars03<-c("smkinit_gp","smkces_gp")
vars05<-c("smkinit_grs_z","smkces_grs_z")
vars06<-c("smkinit_gp_z","smkces_gp_z")

for(i in 1:length(vars01))
  
{
  mod01<-glm(dat[,vars01[i]]~dat[,vars02[i]]+bmi_grs+eduyears_grs+risktk_grs, data=dat, family="binomial")
  dat[,vars03[i]]<-predict.glm(update(mod01,na.action=na.exclude), type="response")
  dat[,vars05[i]]<-as.numeric(with(dat,scale(dat[,vars02[i]])))
  dat[,vars06[i]]<-as.numeric(with(dat,scale(dat[,vars03[i]])))
}


### ANALYSES ###
################

vars01<-c("smkinit","cigday2_z","smkces","agesmk2_z")
vars02<-c("smkinit_gp_z","cigday_gp_z","smkces_gp_z","agesmk_gp_z")
vars03<-c("all","women","men")

tab<-NULL
for(i in 1:length(vars01))
  
{
  datx<-dat
  mod03<-glm(formula=as.factor(covid)~datx[,vars02[i]]+bmi_grs+eduyears_grs+risktk_grs
             +pc01+pc02+pc03+pc04+pc05+pc06+pc07+pc08+pc09+pc10
             +pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+as.factor(genotype_batch),
             data=datx, family="binomial")
  coef03<-risk_se_ic_guapa2(as.numeric(summary(mod03)$coefficients[2,1]),as.numeric(summary(mod03)$coefficients[2,2]))
  pval03<-pval_guapa(as.numeric(summary(mod03)$coefficients[2,4]))
  or03<-exp(as.numeric(summary(mod03)$coefficients[2,1]))
  orlo03<-exp(as.numeric(summary(mod03)$coefficients[2,1])-(z*as.numeric(summary(mod03)$coefficients[2,2])))
  orhi03<-exp(as.numeric(summary(mod03)$coefficients[2,1])+(z*as.numeric(summary(mod03)$coefficients[2,2])))
  
  tab<-rbind(tab,cbind(coef03,pval03,or03,orlo03,orhi03))
  
  datx<-dat[dat$role==0,]
  mod06<-glm(formula=as.factor(covid)~datx[,vars02[i]]+bmi_grs+eduyears_grs+risktk_grs
             +pc01+pc02+pc03+pc04+pc05+pc06+pc07+pc08+pc09+pc10
             +pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+as.factor(genotype_batch),
             data=datx, family="binomial")
  coef06<-risk_se_ic_guapa2(as.numeric(summary(mod06)$coefficients[2,1]),as.numeric(summary(mod06)$coefficients[2,2]))
  pval06<-pval_guapa(as.numeric(summary(mod06)$coefficients[2,4]))
  or06<-exp(as.numeric(summary(mod06)$coefficients[2,1]))
  orlo06<-exp(as.numeric(summary(mod06)$coefficients[2,1])-(z*as.numeric(summary(mod06)$coefficients[2,2])))
  orhi06<-exp(as.numeric(summary(mod06)$coefficients[2,1])+(z*as.numeric(summary(mod06)$coefficients[2,2])))
  
  tab<-rbind(tab,cbind(coef06,pval06,or06,orlo06,orhi06))
  
  datx<-dat[dat$role==1,]
  mod09<-glm(formula=as.factor(covid)~datx[,vars02[i]]+bmi_grs+eduyears_grs+risktk_grs
             +pc01+pc02+pc03+pc04+pc05+pc06+pc07+pc08+pc09+pc10
             +pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+as.factor(genotype_batch),
             data=datx, family="binomial")
  coef09<-risk_se_ic_guapa2(as.numeric(summary(mod09)$coefficients[2,1]),as.numeric(summary(mod09)$coefficients[2,2]))
  pval09<-pval_guapa(as.numeric(summary(mod09)$coefficients[2,4]))
  or09<-exp(as.numeric(summary(mod09)$coefficients[2,1]))
  orlo09<-exp(as.numeric(summary(mod09)$coefficients[2,1])-(z*as.numeric(summary(mod09)$coefficients[2,2])))
  orhi09<-exp(as.numeric(summary(mod09)$coefficients[2,1])+(z*as.numeric(summary(mod09)$coefficients[2,2])))
  
  tab<-rbind(tab,cbind(coef09,pval09,or09,orlo09,orhi09))
  
}

colnames(tab)<-c("MVMR-Coef","MVMR-pval","MVMR-OR","MVMR-lo","MVMR-hi")
rownames(tab)<-head(apply(expand.grid(vars03,vars01), 1, paste, collapse="_"),-1)
write.table(tab,file="./Outputs/Results/mvmr.csv",sep=";",col.names=NA)


