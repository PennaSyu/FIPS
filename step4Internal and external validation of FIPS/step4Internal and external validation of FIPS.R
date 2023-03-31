
#Internal validation:TCGA all----------------
rm(list=ls())
load(file = "../Original Data/step3constrution_riskmodel/step3riskmodel.rdata")
load(file = "../Original Data/step1DataPrepare_DEG/step1DEG.Rdata")
rownames(phe)<-sub(".","-",rownames(phe),fixed = TRUE)
rownames(phe)<-sub(".","-",rownames(phe),fixed = TRUE)
rownames(phe)<-sub(".","-",rownames(phe),fixed = TRUE)
riskmodel$FGFR3 <- phe[rownames(riskmodel),"mutation_FGFR3"]
colnames(riskmodel)[4] <- "time"
colnames(riskmodel)[5] <- "status"
library(survminer)
sur.cut <- surv_cutpoint(riskmodel, time = "time", event = "status",
                         variables = "score")
summary(sur.cut)#0.3526229 

riskmodel$class<-ifelse(riskmodel$score > 0.3526229  ,"high","low")
riskmodel$class<-factor(riskmodel$class,levels = c("low","high"))
save(riskmodel,sur.cut,file = "../Original Data/step4Internal and external validation of FIPS/step4internal validation.rdata")
######log-rank test
dat <- riskmodel
mySurv<-Surv(time = dat$time,event = dat$status)
fit<-surv_fit(mySurv~dat$class,data=dat)

p<-ggsurvplot(fit,data=dat,palette = "lancet",pval = TRUE,
              test.for.trend = FALSE,surv.median.line = "n")
ggpar(p, font.legend = c(6, "bold", "black"),title = "TCGA_FIPS")

#Internal validation:TCGA FGFR3 mutation----------------

riskmodel_m <- riskmodel[riskmodel$FGFR3==1,]
######log-rank test
dat <- riskmodel_m
mySurv<-Surv(time = dat$time,event = dat$status)
fit<-surv_fit(mySurv~dat$class,data=dat)

p<-ggsurvplot(fit,data=dat,palette = "lancet",pval = TRUE,
              test.for.trend = FALSE,surv.median.line = "n")
ggpar(p, font.legend = c(6, "bold", "black"),title = "TCGA_FIPS_mutation")


#Internal validation:TCGA FGFR3 wild----------------

riskmodel_w <- riskmodel[riskmodel$FGFR3==0,]
######log-rank test
dat <- riskmodel_w
mySurv<-Surv(time = dat$time,event = dat$status)
fit<-surv_fit(mySurv~dat$class,data=dat)

p<-ggsurvplot(fit,data=dat,palette = "lancet",pval = TRUE,
              test.for.trend = FALSE,surv.median.line = "n")
ggpar(p, font.legend = c(6, "bold", "black"),title = "TCGA_FIPS_wild")




#External validation in FUSCC cohort-------------------------
rm(list=ls())
test <- read.csv("../Original Data/data_storage/FUSCC.csv")
dat <- test
mySurv<-Surv(time = dat$time,event = dat$status)
fit<-surv_fit(mySurv~dat$class,data=dat)

p<-ggsurvplot(fit,data=dat,palette = "lancet",pval = TRUE,
              test.for.trend = FALSE,surv.median.line = "n")
ggpar(p, font.legend = c(6, "bold", "black"),title = "FUSCC")


#FUSCC ANXA1--------------------------
fit<-surv_fit(mySurv~dat$ANXA1,data=dat)
p<-ggsurvplot(fit,data=dat,palette = "lancet",pval = TRUE,
              test.for.trend = FALSE,surv.median.line = "n")
ggpar(p, font.legend = c(6, "bold", "black"),title = "FUSCC ANXA1")


#FUSCC CLSTN2--------------------------
fit<-surv_fit(mySurv~dat$CLSTN2,data=dat)
p<-ggsurvplot(fit,data=dat,palette = "lancet",pval = TRUE,
              test.for.trend = FALSE,surv.median.line = "n")
ggpar(p, font.legend = c(6, "bold", "black"),title = "FUSCC CLSTN2")

#FUSCC PLEKHG4B--------------------------
fit<-surv_fit(mySurv~dat$PLEKHG4B,data=dat)
p<-ggsurvplot(fit,data=dat,palette = "lancet",pval = TRUE,
              test.for.trend = FALSE,surv.median.line = "n")
ggpar(p, font.legend = c(6, "bold", "black"),title = "FUSCC PLEKHG4B")

