rm(list = ls())
library(ComplexHeatmap)
library(ggsci)

tcgaphe <- read.table("../Original Data/step8clinical_pathologcal_information/TCGAclinical.txt",header = T,row.names = 1)
tcgaphe=tcgaphe[order(tcgaphe$Immune_Score,decreasing = T),]
cli.colors=list()
for(i in 2:ncol(tcgaphe)){
  mycolor=ggsci::pal_jco()(9)
  var=tcgaphe[i]
  var.clas=names(table(tcgaphe[,i]))
  var.color=mycolor[1:length(var.clas)]
  names(var.color)=var.clas
  cli.colors=c(cli.colors,list(var.color))
}
names(cli.colors)=colnames(tcgaphe)[2:ncol(tcgaphe)]
tcga.risk.barplot=columnAnnotation(ImmuneScore = anno_barplot(tcgaphe$Immune_Score,baseline =0,bar_width=1,gp=gpar(fill= '#3F6793',border="#3F6793"),
                                                              border = F,annotation_name_side ='left',height =unit(1,'inches')))


# draw(tcga.risk.barplot)
tcga.cli.heatmap=columnAnnotation(df = tcgaphe[,2:10]
                                  ,annotation_name_side='left'
                                  ,annotation_height =unit(4,'inches')
                                  ,col = cli.colors)
ht_list=tcga.risk.barplot %v% tcga.cli.heatmap
ht_list