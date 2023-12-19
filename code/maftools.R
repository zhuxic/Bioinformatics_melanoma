#if (!require("BiocManager"))
#    install.packages("BiocManager")
#BiocManager::install("maftools")


library(maftools)       #���ð�
setwd("D:\\melanoma\\anoikis\\11.maftools")      #���ù���Ŀ¼

#��ȡͻ������ļ�
geneRT=read.table("gene.txt", header=T, sep="\t", check.names=F, row.names=1)
gene=row.names(geneRT)

#�����ٲ�ͼ
pdf(file="oncoplot.pdf", width=8, height=7.5)
maf=read.maf(maf="input.maf")
oncoplot(maf=maf, genes=gene, fontSize=0.75, draw_titv=T)
dev.off()


