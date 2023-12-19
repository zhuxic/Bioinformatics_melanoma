
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggpubr")


#引用包
library(limma)
library(ggpubr)
prgCluFile="PRGcluster.txt"        #细胞焦亡分型文件
geneCluFile="geneCluster.txt"      #基因分型文件
scoreFile="risk.all.txt"           #风险文件
setwd("D:\\melanoma\\anoikis\\38.clusterRisk")     #设置工作目录

#读取输入文件
prgClu=read.table(prgCluFile, header=T, sep="\t", check.names=F, row.names=1)
geneClu=read.table(geneCluFile, header=T, sep="\t", check.names=F, row.names=1)
score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)

#合并数据
twoCluster=cbind(prgClu, geneClu)
rownames(twoCluster)=gsub("(.*?)\\_(.*?)", "\\2", rownames(twoCluster))
sameSample=intersect(row.names(twoCluster), row.names(score))
data=cbind(score[sameSample,,drop=F], twoCluster[sameSample,,drop=F])


#######细胞焦亡分型的箱线图########
#设置比较组
data$PRGcluster=factor(data$PRGcluster, levels=levels(factor(data$PRGcluster)))
group=levels(factor(data$PRGcluster))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#定义颜色
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data$PRGcluster)))]
	
#绘制boxplot
boxplot=ggboxplot(data, x="PRGcluster", y="riskScore", color="PRGcluster",
			      xlab="PRGcluster",
			      ylab="Risk score",
			      legend.title="PRGcluster",
			      palette=bioCol,
			      add = "jitter")+ 
	stat_compare_means(comparisons = my_comparisons)
	
#输出图片
pdf(file="PRGcluster.pdf", width=5, height=4.5)
print(boxplot)
dev.off()
#######细胞焦亡分型的箱线图########


#######基因分型的箱线图########
#设置比较组
data$geneCluster=factor(data$geneCluster, levels=levels(factor(data$geneCluster)))
group=levels(factor(data$geneCluster))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#定义颜色
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data$geneCluster)))]
	
#绘制boxplot
boxplot=ggboxplot(data, x="geneCluster", y="riskScore", color="geneCluster",
			      xlab="geneCluster",
			      ylab="Risk score",
			      legend.title="geneCluster",
			      palette=bioCol,
			      add = "jitter")+ 
	stat_compare_means(comparisons = my_comparisons)
	
#输出图片
pdf(file="geneCluster.pdf", width=5, height=4.5)
print(boxplot)
dev.off()
#######基因分型的箱线图########




