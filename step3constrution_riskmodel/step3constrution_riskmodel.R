rm(list = ls())
options(stringsAsFactors = FALSE)
getwd()
load(file = "../Original Data/step1DataPrepare_DEG/step1DEG.Rdata")
load(file = "../Original Data/step2GSEA/ImmuneDEG_result.rdata")
#obtain immune realated degs
immunedeg <- ImmuneDEG_union
#1: UniCox model fitted-----------------------------
library(survival)
colnames(phe)
sample_info<-phe[,c("OS.days","Os.status")]
rownames(sample_info) <- sub(".","-",rownames(sample_info),fixed = TRUE)
rownames(sample_info) <- sub(".","-",rownames(sample_info),fixed = TRUE)
rownames(sample_info) <- sub(".","-",rownames(sample_info),fixed = TRUE)
sample_info$ID <- rownames(sample_info)
survival_profile<-as.data.frame(t(RNA_seq_filter))
survival_profile<-survival_profile[rownames(sample_info),immunedeg]
sample_info<-cbind(sample_info,survival_profile)
survival_res<-data.frame(gene=immunedeg)
survival_res$unicox_pvalue<-NA
survival_res$expcoef<-NA
survival_res$expcoef_lower<-NA
survival_res$expcoef_upper<-NA
for (i in 4:613) {
  cox_test<-coxph(Surv(sample_info$OS.days, sample_info$Os.status) ~
                    sample_info[,i] , sample_info)
  cox_test<-summary(cox_test)
  survival_res$unicox_pvalue[i-3]<- as.data.frame(cox_test$coefficients)[1,5]
  survival_res[i-3,3:5]<-cox_test$conf.int[1,c(1,3,4)]
  
}
write.csv(survival_res,file = "../Original Data/step3constrution_riskmodel/unicox.csv")
#2. uni cox, 160 immune related DEGs with prognostic value------------------
survival_res<-survival_res[survival_res$unicox_pvalue<0.05,]
core_DEG<-survival_res$gene[survival_res$unicox_pvalue<0.05]
save(core_DEG,survival_profile,sample_info,phe,file = "../Original Data/step3constrution_riskmodel/step3unicox.Rdata")



#3: lassocox, 160 genes-----------------------------------
rm(list=ls())
load(file = "../Original Data/step3constrution_riskmodel/step3unicox.Rdata")
library(glmnet)
library(survival)
library(survminer)
set.seed(123)
siggene<-core_DEG #160
siggene <-na.omit(siggene)

lasso_profile<-sample_info[,siggene]
lasso_profile<-cbind(lasso_profile,
                     sample_info[rownames(lasso_profile),c("OS.days","Os.status")])

colnames(lasso_profile)[161:162]<-c("time","status")

x<-as.matrix(lasso_profile[,1:160])
y<-as.matrix(lasso_profile[,161:162])

cvfit = cv.glmnet(x, y, type.measure="deviance",
                  alpha=1,
                  family = "cox",nfolds=10)
plot(cvfit)
cvfit$lambda.min#0.0438361
cvfit$lambda.1se#0.1469212

coef.min = coef(cvfit, s ="lambda.min" ) 
active.min = which(coef.min != 0) 
gene= coef.min@Dimnames[[1]][active.min]#20

lassocox_gene<-coef.min@Dimnames[[1]][active.min];lassocox_gene
# [1] "ANXA1"    "IL4I1"    "CD247"    "KLRD1"    "PYHIN1"   "CYP1B1"   "SERPINE2" "CERCAM"  
# [9] "CIITA"    "FKBP10"   "SERPINB2" "F2RL2"    "PCSK5"    "NGF"      "PLEKHG4B" "NTNG1"   
# [17] "CLSTN2"   "DSG1"     "PCDH10"   "MAP1B"  

#4: 100 times lasso cox(20gene)------------------
model20<-lasso_profile[,c(lassocox_gene,"time","status")]
x<-as.matrix(model20[,1:20])
y<-as.matrix(model20[,21:22])
set.seed(119)
lassocox_gene<-c()
for (i in 1:100) {
  cvfit = cv.glmnet(x, y, type.measure="deviance",
                    alpha=1,
                    family = "cox",nfolds=10)
  coef.min = coef(cvfit, s ="lambda.1se" ) 
  active.min = which(coef.min != 0) 
  gene= coef.min@Dimnames[[1]][active.min]
  lassocox_gene<-c(lassocox_gene,gene)
}
lassocox_res<-as.data.frame(table(lassocox_gene));lassocox_res
#6 genes appeared more than 80 times out of 100 repetitions in LASSO Cox analysis.
gene <- c("SERPINB2","PLEKHG4B","CERCAM","CLSTN2","ANXA1","PCSK5")
#ANXA1      CERCAM         CLSTN2      PCSK5      PLEKHG4B   SERPINB2 
#5. AIC: model with the lowest AIC------------------------------
library(survival)
library(MASS)
riskmodel <- sample_info[,gene]
riskmodel[,gene]<-sapply(riskmodel[,gene],scale)

data1905.surv <- Surv(sample_info$OS.days,sample_info$Os.status)

coxph.data1905 <- coxph(data1905.surv ~ SERPINB2+PLEKHG4B+CERCAM+CLSTN2+ANXA1+PCSK5,
                        data=riskmodel)
stepAIC(coxph.data1905,direction = "both")
# coef exp(coef) se(coef)     z       p
# SERPINB2 0.13303   1.14228  0.09372 1.419 0.15577
# PLEKHG4B 0.21491   1.23975  0.07935 2.709 0.00676
# CLSTN2   0.22648   1.25418  0.08070 2.806 0.00501
# ANXA1    0.27938   1.32230  0.10251 2.725 0.00642


#only p<0.05 were included:PLEKHG4B,CLSTN2,ANXA1

#6. Final model: FIPS------------------------------------------
coxph.riskmodel <- coxph(data1905.surv ~ PLEKHG4B+CLSTN2+ANXA1,
                         data=riskmodel)
coxph.riskmodel$coefficients
# PLEKHG4B    CLSTN2     ANXA1 
# 0.2294032 0.2177420 0.3679997 
#5 significant figures reserved
FIPS_score <- "ANXA1*0.36800+CLSTN2*0.21774+PLEKHG4B*0.22940"

riskmodel<-riskmodel[,c("ANXA1","PLEKHG4B","CLSTN2")]
riskmodel<-cbind(riskmodel,sample_info[,c("OS.days","Os.status")])
riskmodel$score<-riskmodel$ANXA1* 0.36800+riskmodel$PLEKHG4B* 0.22940+riskmodel$CLSTN2*0.21774 
save(riskmodel,FIPS_score,file = "../Original Data/step3constrution_riskmodel/step3riskmodel.rdata")




