rm(list=ls())
#Search for FIPS related genes---------------------
load(file = "../Original Data/step4Internal and external validation of FIPS/step4internal validation.rdata")
load(file = "../Original Data/step1DataPrepare_DEG/step1DEG.Rdata")
RNA_seq_filter <- RNA_seq_filter[,rownames(riskmodel)]

gene_name1<-c()
cor_r<-c()
pvalue<-c()
for (i in 1:nrow(RNA_seq_filter)){
  g1=rownames(RNA_seq_filter)[i]
  c_r=cor(as.numeric(RNA_seq_filter[i,]),riskmodel$score,method="pearson")
  p=cor.test(as.numeric(RNA_seq_filter[i,]),riskmodel$score,method ="pearson")[[3]]
  gene_name1=c(gene_name1,g1)
  cor_r=c(cor_r,c_r)
  pvalue=c(pvalue,p)
}

data_cor<-data.frame(gene_name1,cor_r,pvalue)
save(data_cor,file = "../Original Data/step8FIPSnote/data_cor.rdata")

rm(list = ls())
load(file = "../Original Data/step8FIPSnote/data_cor.rdata")
data_cor_sig <- data_cor[data_cor$cor_r > 0.4 | data_cor$cor_r < -0.4, ]
data_cor_sig <- data_cor_sig[data_cor_sig$pvalue<0.05,]
library(org.Hs.eg.db)
gene_ids <- mapIds(org.Hs.eg.db, keys = data_cor_sig$gene_name1, keytype = "SYMBOL", column = "ENTREZID")
library(clusterProfiler)
kk.diff <- enrichKEGG(gene= gene_ids,
                        organism = 'hsa',
                      pvalueCutoff = 0.05)
kegg_result <- kk.diff@result
library(enrichplot)
#Figure7a---------------------------
dotplot(kk.diff,
        x = "Count",
        color = "pvalue",
        showCategory = 20,
        size = 'GeneRatio',
        split = NULL,
        font.size = 12,
        title = "KEGG: FIPS related genes",
        orderBy = "x",
        label_format = 30)


#Figure7b draw heatmap----------------------
load(file = "../Original Data/step4Internal and external validation of FIPS/step4internal validation.rdata")
load(file = "../Original Data/step1DataPrepare_DEG/step1DEG.Rdata")
data_cor_heat <- data_cor[data_cor$cor_r > 0.5 | data_cor$cor_r < -0.5, ]
data_cor_heat <- data_cor_heat[data_cor_heat$pvalue<0.05,]
exp <- RNA_seq_filter[data_cor_heat$gene_name1,]
rownames(phe)<-sub(".","-",rownames(phe),fixed = TRUE)
rownames(phe)<-sub(".","-",rownames(phe),fixed = TRUE)
rownames(phe)<-sub(".","-",rownames(phe),fixed = TRUE)
exp <- exp[,rownames(phe)]
Group <- phe$score_class


library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), c("#65b5b5", "white", "#CC4C44"))
top_annotation = HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = c("#d74b12","#74bfda")),
                       labels = c("FIPShigh","FIPSlow"),
                       labels_gp = gpar(col = "white", fontsize = 12)))
Heatmap(t(scale(t(exp))),name = " ",
            col = col_fun,
            top_annotation = top_annotation,
            column_split = Group,
            show_heatmap_legend = T,
            border = F,
            show_column_names = F,
            show_row_names = F,
            column_title = NULL,
        cluster_columns = F )


