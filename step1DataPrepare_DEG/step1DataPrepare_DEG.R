rm(list = ls())
options(stringsAsFactors = FALSE)
getwd()
###1. Data prepare-------------------
RNA_seq<-read.table("../Original Data/data_storage/TCGA_Express_Data403.txt",
                    header = TRUE,row.names = 1,sep = "\t")
phe <- read.table("../Original Data/data_storage/TCGA_clinical_399.txt",
                  header = TRUE,row.names = 1,sep = "\t")
RNA_seq <- RNA_seq[,colnames(RNA_seq) %in% rownames(phe)]
colnames(RNA_seq)<-sub(".","-",colnames(RNA_seq),fixed = TRUE)
colnames(RNA_seq)<-sub(".","-",colnames(RNA_seq),fixed = TRUE)
colnames(RNA_seq)<-sub(".","-",colnames(RNA_seq),fixed = TRUE)
#filter gene：
RNA_seq$mean<-apply(RNA_seq,1,mean,na.rm=TRUE)
RNA_seq_filter<-RNA_seq[RNA_seq$mean>1,]
RNA_seq_filter <- RNA_seq_filter[,-400]

#2. DEG analysis-------------------
#####399 samples:17106 genes,find DEG
mutation_sample <- rownames(phe)[phe$FGFR_status=="mutation"]
mutation_sample <- sub(".","-",mutation_sample,fixed = TRUE)
mutation_sample <- sub(".","-",mutation_sample,fixed = TRUE)
mutation_sample <- sub(".","-",mutation_sample,fixed = TRUE)


RNA_seq_filter<-RNA_seq_filter[,c(mutation_sample,
                                  colnames(RNA_seq_filter)[!colnames(RNA_seq_filter) %in% mutation_sample])]
library(limma)
###limma
fpkm<-as.matrix(RNA_seq_filter)
#Group
sample_lable<-factor(c(rep("mutation",55),rep("wild",344)))
design<-model.matrix(~0+sample_lable)
colnames(design)<-c("mutation","wild")
dim(fpkm)
#fit
fit<-lmFit(fpkm,design,
           weights=c(rep(6,55),rep(1,344)))
#contrast model
contrast.matrix<-makeContrasts(mutation-wild,levels = design)
fit2<-contrasts.fit(fit,contrast.matrix)
#eBayes
fit2<-eBayes(fit2,trend = TRUE)
limmaresult<-topTable(fit2,coef = 1,number = dim(fpkm)[1],
                      adjust.method = "fdr", sort.by="M",
                      resort.by = "P")   
write.table(limmaresult,"../Original Data/step1DataPrepare_DEG/limmaresult.txt",sep = "\t")


#DEG:P.Value<0.05 & abs(limmaresult$logFC)>1
###################################################################
DEG<-limmaresult[limmaresult$adj.P.Val<0.05 & 
                   abs(limmaresult$logFC)>1,]
DEG$abs<-abs(DEG$logFC)
DEG$DEG_type[DEG$logFC>0]<-"UP"
DEG$DEG_type[DEG$logFC<0]<-"Down"
write.table(DEG,"../Original Data/step1DataPrepare_DEG/DEG_limmaresult.txt",sep = "\t")
#updata:1811DEG (adj.P.Val<0.05 & abs(DEG$logFC)>1)
DEG<-rownames(DEG)

save(RNA_seq_filter,phe,DEG,mutation_sample,limmaresult,file = "../Original Data/step1DataPrepare_DEG/step1DEG.Rdata")


#3. Draw volcano plot----------------
library(ggVolcano)
limmaresult1 <- limmaresult
limmaresult1 <- cbind(row=rownames(limmaresult1),limmaresult1)
limmaresult1_test<- add_regulate(limmaresult1, log2FC_name = "logFC",
                     fdr_name = "adj.P.Val",log2FC = 1, fdr = 0.05)
# plot
ggvolcano(limmaresult1_test, x = "log2FoldChange",y = "padj",output = T,
          add_label = F,
          y_lab = "-log(adj.P.Val)",
          x_lab = "logFC",legend_position = "DR",
          legend_title = "",
          fills = c("#2166ac", "#4d4d4d", "#d6604d"),
          colors = c("#2166ac", "#4d4d4d", "#d6604d"),
          pointSize = 2,
          pointShape = 20,filename = "volcano_plot")
