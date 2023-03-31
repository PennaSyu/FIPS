#Figure 6a-----------------
rm(list = ls())
load("../Original Data/step4Internal and external validation of FIPS/step4internal validation.rdata")
colnames(riskmodel)[6] <- "FIPS"
TCGAimmune <- read.table("../Original Data/data_storage/TCGAimmuneEstimation_allcancers.txt",header = T,row.names = 1)
TCGAimmune1 <- TCGAimmune[rownames(riskmodel),1:6]
riskmodel1 <- cbind(riskmodel,TCGAimmune1)
colnames(riskmodel1)
mydata <- riskmodel1[,c("FIPS","B_cell","CD4_Tcell","CD8_Tcell","Neutrophil","Macrophage","Dendritic")]
head(mydata)
library(psych)
res <- corr.test(mydata, method = 'spearman')#pearson,spearman
# Printing the correlation matrix
signif(res$r, 2)
# Printing the p-values of the correlations
signif(res$p,2)
#install.packages("PerformanceAnalytics")
library(PerformanceAnalytics)
chart.Correlation(mydata, histogram=TRUE, pch=19,method = "spearman")

#Figure 6b-------------------------
rm(list = ls())
library(tidyr)
library(dplyr)
library(tibble)
library(ggpubr)
load("../Original Data/step4Internal and external validation of FIPS/step4internal validation.rdata")
colnames(riskmodel)[6] <- "FIPS"
TCGAimmune <- read.table("../Original Data/data_storage/TCGAimmuneEstimation_allcancers.txt",header = T,row.names = 1)
TCGAimmune1 <- TCGAimmune[rownames(riskmodel),1:6]
riskmodel1 <- cbind(riskmodel,TCGAimmune1)
mydata <- riskmodel1[,c("class","B_cell","CD4_Tcell","CD8_Tcell","Neutrophil","Macrophage","Dendritic")]


p1 <- ggplot(data=mydata,aes(x=class,y=B_cell))+
  geom_violin(aes(fill=class))+
  geom_boxplot(width=0.3,outlier.size = 0)+
  scale_fill_manual(values = c("#3c8068", "#FF6245"))+
  ylim(0,1)+
  theme(legend.position = "none")

wilcox.test(mydata$B_cell~mydata$class)
p2 <- ggplot(data=mydata,aes(x=class,y=CD4_Tcell))+
  geom_violin(aes(fill=class))+
  geom_boxplot(width=0.3,outlier.size = 0)+
  scale_fill_manual(values = c("#3c8068", "#FF6245"))+
  ylim(0,1)+
  theme(legend.position = "none")
wilcox.test(mydata$CD4_Tcell~mydata$class)

p3 <- ggplot(data=mydata,aes(x=class,y=CD8_Tcell))+
  geom_violin(aes(fill=class))+
  geom_boxplot(width=0.3,outlier.size = 0)+
  scale_fill_manual(values = c("#3c8068", "#FF6245"))+
  ylim(0,1)+
  theme(legend.position = "none")
wilcox.test(mydata$CD8_Tcell~mydata$class)
p4 <- ggplot(data=mydata,aes(x=class,y=Neutrophil))+
  geom_violin(aes(fill=class))+
  geom_boxplot(width=0.3,outlier.size = 0)+
  scale_fill_manual(values = c("#3c8068", "#FF6245"))+
  ylim(0,1)+
  theme(legend.position = "none")
wilcox.test(mydata$Neutrophil~mydata$class)

p5 <- ggplot(data=mydata,aes(x=class,y=Macrophage))+
  geom_violin(aes(fill=class))+
  geom_boxplot(width=0.3,outlier.size = 0)+
  scale_fill_manual(values = c("#3c8068", "#FF6245"))+
  ylim(0,1)+
  theme(legend.position = "none")
wilcox.test(mydata$Macrophage~mydata$class)

p6 <- ggplot(data=mydata,aes(x=class,y=Dendritic))+
  geom_violin(aes(fill=class))+
  geom_boxplot(width=0.3,outlier.size = 0)+
  scale_fill_manual(values = c("#3c8068", "#FF6245"))+
  ylim(0,1)
wilcox.test(mydata$Dendritic~mydata$class)

library(patchwork)
pdf(file = "../Original Data/step7FIPS_immune_microenvironment/immunecell_class.pdf",width = 20,height = 4)
p1+p2+p3+p4+p5+p6+plot_layout(nrow=1)
dev.off()


#Figure 6c-------------------------
rm(list = ls())
load("../Original Data/step4Internal and external validation of FIPS/step4internal validation.rdata")
colnames(riskmodel)[6] <- "FIPS"
load("../Original Data/step1DataPrepare_DEG/step1DEG.Rdata")
checkpoint <- c("PDCD1","CD274","CTLA4",'LAG3','HAVCR2')
cpexp <- t(RNA_seq_filter)
cpexp1 <- cpexp[rownames(riskmodel),checkpoint]

riskmodel1 <- cbind(FIPS=riskmodel$FIPS,cpexp1)
library(psych)
res <- corr.test(riskmodel1, method = 'spearman')#pearson,spearman
res$r
res$p

#相关性散点图
riskmodel2 <- as.data.frame(cbind(class=as.data.frame(riskmodel$class),FIPS=riskmodel$FIPS,cpexp1))
colnames(riskmodel2)
colnames(riskmodel2)[1] <- "class"
p1 <- ggplot(data=riskmodel2,mapping=aes(x=FIPS,y=PDCD1))+
  geom_point(aes(color=class))+
  scale_color_manual(values = c("#004d9F","#FB5446"))+
  stat_smooth(method = "lm")+
  stat_cor(data=riskmodel2, method = "spearman",p.accuracy = 0.001)+
  theme(legend.position = "none");p1

p2 <- ggplot(data=riskmodel2,mapping=aes(x=FIPS,y=CD274))+
  geom_point(aes(color=class))+
  scale_color_manual(values = c("#004d9F","#FB5446"))+
  stat_smooth(method = "lm")+
  stat_cor(data=riskmodel2, method = "spearman",p.accuracy = 0.001)+
  theme(legend.position = "none")
p3 <- ggplot(data=riskmodel2,mapping=aes(x=FIPS,y=CTLA4))+
  geom_point(aes(color=class))+
  scale_color_manual(values = c("#004d9F","#FB5446"))+
  stat_smooth(method = "lm")+
  stat_cor(data=riskmodel2, method = "spearman",p.accuracy = 0.001)+
  theme(legend.position = "none")
p4 <- ggplot(data=riskmodel2,mapping=aes(x=FIPS,y=LAG3))+
  geom_point(aes(color=class))+
  scale_color_manual(values = c("#004d9F","#FB5446"))+
  stat_smooth(method = "lm")+
  stat_cor(data=riskmodel2, method = "spearman",p.accuracy = 0.001)+
  theme(legend.position = "none")
p5 <- ggplot(data=riskmodel2,mapping=aes(x=FIPS,y=HAVCR2))+
  geom_point(aes(color=class))+
  scale_color_manual(values = c("#004d9F","#FB5446"))+
  stat_smooth(method = "lm")+
  stat_cor(data=riskmodel2, method = "spearman",p.accuracy = 0.001)
library(patchwork)
pdf(file = "../Original Data/step7FIPS_immune_microenvironment/checkpointFIPS.pdf",width = 15,height = 4)
p1+p2+p3+p4+p5+plot_layout(nrow=1)
dev.off()
