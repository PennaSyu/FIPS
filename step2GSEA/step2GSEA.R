rm(list = ls())
options(stringsAsFactors = FALSE)
getwd()
load(file = "../Original Data/step1DataPrepare_DEG/step1DEG.Rdata")
library(clusterProfiler)
library(enrichplot)
library(stringr)
#1:GSEA--------------------------

GSEA_input<-limmaresult
GSEA_input<-GSEA_input[order(GSEA_input$logFC,decreasing = TRUE),]

GSEA_gene<-bitr(rownames(GSEA_input),
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb="org.Hs.eg.db")
GSEA_gene$lfc<-GSEA_input[GSEA_gene$SYMBOL,"logFC"]

GSEA_lfc<-GSEA_gene$lfc
attr(GSEA_lfc,"names")<-GSEA_gene$ENTREZID
#####c5.go.bp.v2023.1.Hs.entrez.gmt
C5BP<-read.gmt("../Original Data/data_storage/c5.go.bp.v2023.1.Hs.entrez.gmt")
set.seed(123)
GSEA_res<-GSEA(GSEA_lfc,TERM2GENE=C5BP, pvalueCutoff = 0.05,
                verbose=TRUE)
GSEA_res2<-GSEA_res@result
write.table(GSEA_res2,"../Original Data/step2GSEA/GSEA_result.txt",sep = "\t")

#draw picture
p1 <- gseaplot2(GSEA_res,geneSetID = "GOBP_NEUTROPHIL_MIGRATION",
          title = "GOBP_NEUTROPHIL_MIGRATION",
          pvalue_table = F,
          rel_heights = c(0.8,0.2,0.5),
          ES_geom = "line",
          );p1
p2 <- gseaplot2(GSEA_res,geneSetID = "GOBP_GRANULOCYTE_MIGRATION",
                title = "GOBP_GRANULOCYTE_MIGRATION",
                pvalue_table = F,
                rel_heights = c(0.8,0.2,0.5),
                ES_geom = "line",
);p2
p3 <- gseaplot2(GSEA_res,geneSetID = "GOBP_MONOCYTE_CHEMOTAXIS",
                title = "GOBP_MONOCYTE_CHEMOTAXIS",
                pvalue_table = F,
                rel_heights = c(0.8,0.2,0.5),
                ES_geom = "line",
);p3
p4 <- gseaplot2(GSEA_res,geneSetID = "GOBP_MYELOID_LEUKOCYTE_MIGRATION",
                title = "GOBP_MYELOID_LEUKOCYTE_MIGRATION",
                pvalue_table = F,
                rel_heights = c(0.8,0.2,0.5),
                ES_geom = "line",
);p4
p5 <- gseaplot2(GSEA_res,geneSetID = "GOBP_LYMPHOCYTE_MEDIATED_IMMUNITY",
                title = "GOBP_LYMPHOCYTE_MEDIATED_IMMUNITY",
                pvalue_table = F,
                rel_heights = c(0.8,0.2,0.5),
                ES_geom = "line",
);p5

p6 <- gseaplot2(GSEA_res,geneSetID = "GOBP_TYPE_II_INTERFERON_PRODUCTION",
                title = "GOBP_TYPE_II_INTERFERON_PRODUCTION",
                pvalue_table = F,
                rel_heights = c(0.8,0.2,0.5),
                ES_geom = "line",
);p6
pdf("GSEA1.pdf",width=6,height=6,)
p1;p2;p3
p4;p5;p6
dev.off()


 #2:359 immune related pathways selected in the folowing analysis#################
GSEA_res2 <- read.csv(file="../Original Data/data_storage/immune_GSEA_result.csv",row.names = "ID")
#3: immune related DEGs-----------------
##3.1 immune related core genes---------------
immu_core_gene<-c()
for (i in 1:359) {
  immu_core_gene<-c(immu_core_gene,strsplit(as.character(GSEA_res2$core_enrichment[i]),"/")[[1]])
}
immu_core_gene<-immu_core_gene[!duplicated(immu_core_gene)]
immu_core_gene<-bitr(immu_core_gene,
                fromType = "ENTREZID",
                toType = "SYMBOL",
                OrgDb="org.Hs.eg.db")

write.table(immu_core_gene,"../Original Data/step2GSEA/359immu_core_gene.txt",sep = "\t")

##3.2 common genes between DEGs and immune related genes------------
ImmuneDEG_union <- intersect(immu_core_gene$SYMBOL,DEG)
save(immu_core_gene,ImmuneDEG_union,file = "../Original Data/step2GSEA/ImmuneDEG_result.rdata")

