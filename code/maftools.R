#if (!require("BiocManager"))
#    install.packages("BiocManager")
#BiocManager::install("maftools")


library(maftools)       #引用包
setwd("D:\\melanoma\\anoikis\\11.maftools")      #设置工作目录

#读取突变基因文件
geneRT=read.table("gene.txt", header=T, sep="\t", check.names=F, row.names=1)
gene=row.names(geneRT)

#绘制瀑布图
pdf(file="oncoplot.pdf", width=8, height=7.5)
maf=read.maf(maf="input.maf")
oncoplot(maf=maf, genes=gene, fontSize=0.75, draw_titv=T)
dev.off()



