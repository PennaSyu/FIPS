rm(list = ls())
load("../Original Data/step4Internal and external validation of FIPS/step4internal validation.rdata")
#Figure5a--------------------------
library(ggplot2)
library(ggplot2)
riskmodel1 <- riskmodel[order(riskmodel$score,decreasing = T),]
riskmodel1$number <- 1:nrow(riskmodel1)

riskmodel1$FGFR3 <- factor(riskmodel1$FGFR3,levels = c(1,0))
ggplot(riskmodel1,aes(x = number,y = score,fill = factor(FGFR3)))+
  geom_bar(stat = "identity",width = 1)+
  theme_bw()+
  theme(legend.position = c(0.5,0.8),
        panel.grid = element_blank())

#Figure5b-------------------------
pval <- round(t.test(score ~ FGFR3, riskmodel1)$p.value, 10)
ggplot(riskmodel1, aes(x = FGFR3, y = score, fill = FGFR3)) +
  geom_boxplot()+
  annotate("text", x = c(1.6, 2.6), y = tapply(riskmodel1$score, riskmodel1$FGFR3,median),
           label = tapply(riskmodel1$score, riskmodel1$FGFR3, median), size = 3.5)+
  annotate("text", x = 1.5, y = max(riskmodel1$score) + 1, 
           label = paste0("P = ", pval))

