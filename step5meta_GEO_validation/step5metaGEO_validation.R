#GSE13507----------------------------------
rm(list = ls())
library(GEOquery)
gse_number = "GSE13507"
eSet <- getGEO(gse_number,
               destdir = '../Original Data/data_storage/', 
               getGPL = F)
class(eSet)
length(eSet)
eSet1 = eSet[[1]]
class(eSet1)
#1.1exprs()
exp <- exprs(eSet1)
exp <- as.data.frame(exp)
boxplot(exp[1:200,1:10])
pd <- read.csv(file = "../Original Data/data_storage/pd13507.csv",header = T,row.names = 1)
p = identical(rownames(pd),colnames(exp));p
exp <- exp[,rownames(pd)]
if(!p) exp = exp[,match(rownames(pd),colnames(exp))]
#1.2GPL
gpl_number <- eSet1@annotation
#2.ids
library(tinyarray)
find_anno(gpl_number)
ids <- AnnoProbe::idmap('GPL6102',destdir = '../Original Data/data_storage/')
exp$probe <- rownames(exp)
exp$probe <- ids$symbol[match(exp$probe,ids$probe_id)]
exp <- na.omit(exp)
library(dplyr)
exp <- exp %>% 
  group_by(probe) %>% 
  summarise_all(max)

exp <- as.data.frame(exp)
rownames(exp) <- exp$probe
exp <- exp[,-which(colnames(exp) == "probe")]
exp <- exp[,rownames(pd)]
pd$batch <- 1
exp13507 <- exp
pd13507 <- pd
##save data----------------------
save(exp13507,pd13507,file = "../Original Data/step5meta_GEO_validation/GSE13507.rdara")

#GSE31684----------------------------------
rm(list = ls())
library(GEOquery)
gse_number = "GSE31684"
eSet <- getGEO(gse_number,
               destdir = '../Original Data/data_storage/', 
               getGPL = F)
class(eSet)
length(eSet)
eSet1 = eSet[[1]]
class(eSet1)
#1.1exp
exp <- exprs(eSet1)
exp <- as.data.frame(exp)
exp[1:4,1:4]
boxplot(exp[1:200,1:10])
pd <- read.csv(file = "../Original Data/data_storage/pd31684.csv",header = T,row.names = 1)
p = identical(rownames(pd),colnames(exp));p
exp <- exp[,rownames(pd)]
p = identical(rownames(pd),colnames(exp));p
if(!p) exp = exp[,match(rownames(pd),colnames(exp))]
gpl_number <- eSet1@annotation
#2.ids
library(tinyarray)
find_anno(gpl_number)
ids <- AnnoProbe::idmap('GPL570',destdir = '../Original Data/data_storage/')
exp$probe <- rownames(exp)
exp$probe <- ids$symbol[match(exp$probe,ids$probe_id)]
exp <- na.omit(exp)
library(dplyr)
exp <- exp %>% 
  group_by(probe) %>% 
  summarise_all(max)
exp <- as.data.frame(exp)
rownames(exp) <- exp$probe
exp <- exp[,-which(names(exp) == "probe")]
exp <- exp[,rownames(pd)]
pd$batch <- 2
exp31684 <- exp
pd31684<- pd
##save data----------------------
save(exp31684,pd31684,file = "../Original Data/step5meta_GEO_validation/GSE31684.rdara")




#GSE48075----------------------------------
rm(list = ls())
library(GEOquery)
gse_number = "GSE48075"
eSet <- getGEO(gse_number,
               destdir = '../Original Data/data_storage/', 
               getGPL = F)
class(eSet)
length(eSet)
eSet1 = eSet[[1]]
class(eSet1)
#1.1exp
exp <- exprs(eSet1)
exp <- as.data.frame(exp)
exp[1:4,1:4]
boxplot(exp[1:200,1:10])
pd <- read.csv(file = "../Original Data/data_storage/pd48075.csv",header = T,row.names = 1)

p = identical(rownames(pd),colnames(exp));p
exp <- exp[,rownames(pd)]
p = identical(rownames(pd),colnames(exp));p
if(!p) exp = exp[,match(rownames(pd),colnames(exp))]
#1.2GPL
gpl_number <- eSet1@annotation
#2.ids
library(tinyarray)
find_anno(gpl_number)
ids <- AnnoProbe::idmap('GPL6947',
                        destdir = '../Original Data/data_storage/')
exp$probe <- rownames(exp)
exp$probe <- ids$symbol[match(exp$probe,ids$probe_id)]
exp <- na.omit(exp)
library(dplyr)
exp <- exp %>% 
  group_by(probe) %>% 
  summarise_all(max)
exp <- as.data.frame(exp)
rownames(exp) <- exp$probe
exp <- exp[,-which(names(exp) == "probe")]
exp <- exp[,rownames(pd)]
pd$batch <- 3
exp48075 <- exp
pd48075 <- pd
##save----------------------
save(exp48075,pd48075,file = "../Original Data/step5meta_GEO_validation/GSE48075.rdara")
#GSE48276----------------------------------
#数据集下载
rm(list = ls())
library(GEOquery)
gse_number = "GSE48276"
eSet <- getGEO(gse_number,
               destdir = '../Original Data/data_storage/', 
               getGPL = F)
class(eSet)
length(eSet)
eSet1 = eSet[[1]]
class(eSet1)
#(1)exprs()
exp <- exprs(eSet1)
exp <- as.data.frame(exp)
exp[1:4,1:4]
boxplot(exp[1:200,1:10])
#pd <- pData(eSet1)
#write.csv(pd,file = "pd48726qq.csv")
pd <- read.csv(file = "../Original Data/data_storage/pd48276.csv",header = T,row.names = 1)

p = identical(rownames(pd),colnames(exp));p
exp <- exp[,rownames(pd)]
p = identical(rownames(pd),colnames(exp));p
if(!p) exp = exp[,match(rownames(pd),colnames(exp))]
#1.2GPL
gpl_number <- eSet1@annotation
#2.ids
anno <- data.table::fread("../Original Data/data_storage/GPL14951-11332.txt")
colnames(anno)
ids  <-  data.frame(probe_id= anno$ID,
                    symbol  = anno$Symbol)

exp$probe <- rownames(exp)
exp$probe <- ids$symbol[match(exp$probe,ids$probe_id)]
exp <- na.omit(exp)
library(dplyr)
exp <- exp %>% 
  group_by(probe) %>% 
  summarise_all(max)

exp <- as.data.frame(exp)
rownames(exp) <- exp$probe
exp <- exp[,-which(names(exp) == "probe")]
exp <- exp[,rownames(pd)]
pd$batch <- 4
exp48276 <- exp
pd48276<- pd
##save data----------------------
save(exp48276,pd48276,file = "../Original Data/step5meta_GEO_validation/GSE48276.rdara")





#sva combat-----------------------
rm(list=ls())
load(file = "../Original Data/step5meta_GEO_validation/GSE13507.rdara")
load(file = "../Original Data/step5meta_GEO_validation/GSE31684.rdara")
load(file = "../Original Data/step5meta_GEO_validation/GSE48075.rdara")
load(file = "../Original Data/step5meta_GEO_validation/GSE48276.rdara")
library(sva)
library(tidyverse)

gene <- c("ANXA1","CLSTN2","PLEKHG4B")
gene %in% rownames(exp13507)
gene %in% rownames(exp31684)
gene %in% rownames(exp48075)
gene %in% rownames(exp48276)

x1 <- Reduce(intersect,list(
  rownames(exp13507),
  rownames(exp31684),
  rownames(exp48075),
  rownames(exp48276)
))
exp <- cbind(exp13507[x1,],exp31684[x1,],exp48075[x1,],exp48276[x1,])
pd <- rbind(pd13507,pd31684,pd48075,pd48276)
combat_edata <- ComBat(dat = exp, batch = pd$batch)
boxplot(combat_edata)
######################Meta GEO#########################################
#Validation in meta GEO cohort--------------------------------------
GEO1 <- combat_edata
GEO1<-as.data.frame(t(GEO1))
GEO1_gene<-GEO1[,c("ANXA1","CLSTN2","PLEKHG4B")]
GEO1_gene[,1:3]<-sapply(GEO1_gene[,1:3],scale)

GEO1_gene$score<-GEO1_gene$ANXA1*0.36800+GEO1_gene$CLSTN2*0.21774+GEO1_gene$PLEKHG4B*0.22940
#Use the same cutoff value as the TCGA database
GEO1_gene$score_class<-ifelse(GEO1_gene$score>=0.3526229,"high","low")
GEO1_gene$score_class <- factor(GEO1_gene$score_class,ordered = T,levels = c("low","high"))
GEO1_gene$OS<-pd[rownames(GEO1_gene),"os"]
GEO1_gene$OSstatus<-pd[rownames(GEO1_gene),"osstatus"]

###survival plot
library(survival)
library(survminer)
mySurv<-Surv(time = GEO1_gene$OS,
             event = GEO1_gene$OSstatus)
fit<-surv_fit(mySurv~GEO1_gene$score_class,
              data=GEO1_gene)
p<-ggsurvplot(fit,data=GEO1_gene,palette = "lancet",pval = T,
              test.for.trend = FALSE,surv.median.line = "v",xlab = 'Time(month)')
ggpar(p, font.legend = c(6, "bold", "black"),title = "meta-GEO cohort")


summary(coxph(mySurv~GEO1_gene$score_class,
              data=GEO1_gene))
#HE(95%CI)：1.29(1.01~1.64)   P=0.043