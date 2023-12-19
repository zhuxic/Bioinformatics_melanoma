
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggpubr")


#���ð�
library(limma)
library(ggpubr)
prgCluFile="PRGcluster.txt"        #ϸ�����������ļ�
geneCluFile="geneCluster.txt"      #��������ļ�
scoreFile="risk.all.txt"           #�����ļ�
setwd("D:\\melanoma\\anoikis\\38.clusterRisk")     #���ù���Ŀ¼

#��ȡ�����ļ�
prgClu=read.table(prgCluFile, header=T, sep="\t", check.names=F, row.names=1)
geneClu=read.table(geneCluFile, header=T, sep="\t", check.names=F, row.names=1)
score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)

#�ϲ�����
twoCluster=cbind(prgClu, geneClu)
rownames(twoCluster)=gsub("(.*?)\\_(.*?)", "\\2", rownames(twoCluster))
sameSample=intersect(row.names(twoCluster), row.names(score))
data=cbind(score[sameSample,,drop=F], twoCluster[sameSample,,drop=F])


#######ϸ���������͵�����ͼ########
#���ñȽ���
data$PRGcluster=factor(data$PRGcluster, levels=levels(factor(data$PRGcluster)))
group=levels(factor(data$PRGcluster))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#������ɫ
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data$PRGcluster)))]
	
#����boxplot
boxplot=ggboxplot(data, x="PRGcluster", y="riskScore", color="PRGcluster",
			      xlab="PRGcluster",
			      ylab="Risk score",
			      legend.title="PRGcluster",
			      palette=bioCol,
			      add = "jitter")+ 
	stat_compare_means(comparisons = my_comparisons)
	
#���ͼƬ
pdf(file="PRGcluster.pdf", width=5, height=4.5)
print(boxplot)
dev.off()
#######ϸ���������͵�����ͼ########


#######������͵�����ͼ########
#���ñȽ���
data$geneCluster=factor(data$geneCluster, levels=levels(factor(data$geneCluster)))
group=levels(factor(data$geneCluster))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#������ɫ
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data$geneCluster)))]
	
#����boxplot
boxplot=ggboxplot(data, x="geneCluster", y="riskScore", color="geneCluster",
			      xlab="geneCluster",
			      ylab="Risk score",
			      legend.title="geneCluster",
			      palette=bioCol,
			      add = "jitter")+ 
	stat_compare_means(comparisons = my_comparisons)
	
#���ͼƬ
pdf(file="geneCluster.pdf", width=5, height=4.5)
print(boxplot)
dev.off()
#######������͵�����ͼ########



