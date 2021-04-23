rm(list = ls())
options(stringsAsFactors = F)
getwd()
setwd(dir = '/home/data/ylchen/RA-DRUG/code')

gene <- read.table('../data/DMP.Sig.Gene.txt')
head(gene)

# change ID
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
df <- bitr(unique(gene$V1), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)
head(df)

DEG <- gene
head(DEG)
DEG=merge(DEG,df,by.y='SYMBOL',by.x='V1')
head(DEG)
dim(DEG)
gene_all <- as.character(DEG[ ,'ENTREZID'] )

# enrich
#source('kegg_and_go_up_and_down.R')
kegg <- enrichKEGG(gene_all, pvalueCutoff = 0.05,organism   = 'hsa',qvalueCutoff  = 0.05)

#dotplot(kegg,showCategory = 10) don'know why not work

library(DOSE)
r <- kegg@result
head(r)
r_0.05 <- r[r$pvalue<0.05,]
# 将分数generatio转成小数
r_0.05$GeneRatio <- parse_ratio(r_0.05$GeneRatio)
jpeg("../pic/kegg_dotplot.jpg",units="in", width=7, height=5,res=650)
## 绘制气泡图
p <- ggplot(r_0.05, aes(x = GeneRatio, y = Description, size = Count, color=pvalue)) + 
  geom_point() +xlab("GeneRatio") +ylab(" ")
## 修改气泡颜色
p + scale_color_gradient(low='red',high='blue') +
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank())
dev.off()
