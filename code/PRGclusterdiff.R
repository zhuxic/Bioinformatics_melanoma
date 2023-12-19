

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("reshape2")
#install.packages("ggpubr")


#引用包
library(limma)
library(reshape2)
library(ggpubr)

expFile="prgGeneExp.txt"          #表达数据文件
geneCluFile="geneCluster.txt"     #基因分型的结果文件
setwd("D:\\melanoma\\anoikis\\35.prgClusterDiff")       #设置工作目录

#读取表达输入文件
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=t(data)

#读取基因分型文件
geneClu=read.table(geneCluFile, header=T, sep="\t", check.names=F, row.names=1)

#合并数据
sameSample=intersect(row.names(data), row.names(geneClu))
expClu=cbind(data[sameSample,,drop=F], geneClu[sameSample,,drop=F])

#提取差异显著的基因
sigGene=c()
for(i in colnames(expClu)[1:(ncol(expClu)-1)]){
	if(sd(expClu[,i])<0.001){next}
	if(length(levels(factor(expClu[,"geneCluster"])))>2){
		test=kruskal.test(expClu[,i] ~ expClu[,"geneCluster"])
	}else{
		test=wilcox.test(expClu[,i] ~ expClu[,"geneCluster"])
	}
	pvalue=test$p.value
	if(pvalue<0.05){
		sigGene=c(sigGene, i)
	}
}
sigGene=c(sigGene, "geneCluster")
expClu=expClu[,sigGene]

#把数据转换成ggplot2输入文件
data=melt(expClu, id.vars=c("geneCluster"))
colnames(data)=c("geneCluster", "Gene", "Expression")

#设置颜色
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data[,"geneCluster"])))]

#绘制箱线图
p=ggboxplot(data, x="Gene", y="Expression", color = "geneCluster",
	     xlab="",
	     ylab="Gene expression",
	     legend.title="geneCluster",
	     palette = bioCol,
	     width=1)
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=geneCluster),
	      symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
	      label = "p.signif")

#输出箱线图
pdf(file="boxplot.pdf", width=9, height=6)
print(p1)
dev.off()



