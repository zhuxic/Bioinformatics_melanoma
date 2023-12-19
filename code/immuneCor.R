

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("tidyverse")
#install.packages("ggplot2")
#install.packages("ggExtra")
#install.packages("ggpubr")


#引用包
library(limma)
library(reshape2)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggExtra)

immFile="CIBERSORT-Results.txt"     #免疫细胞浸润的结果文件
riskFile="risk.all.txt"             #风险文件
setwd("D:\\melanoma\\anoikis\\44.CIBERSORT")     #设置工作目录

#读取免疫细胞结果文件，并对数据进行整理
immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)
immune=immune[immune[,"P-value"]<0.05,]
data=as.matrix(immune[,1:(ncol(immune)-3)])
rownames(data)=gsub("(.*?)\\_(.*?)", "\\2", rownames(data))

#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(data), row.names(risk))
data=data[sameSample,,drop=F]
risk=risk[sameSample,,drop=F]

#对所有免疫细胞进行循环，得到风险打分与免疫细胞的相关性
for(i in colnames(data)[1:ncol(data)]){
	x=as.numeric(risk[,"riskScore"])
	x[x>quantile(x,0.99)]=quantile(x,0.99)
	y=as.numeric(data[,i])
	if(sd(y)<0.01){next}
	cor=cor.test(x, y, method="spearma")
	#绘制相关性散点图
	if(cor$p.value<0.05){
		outFile=paste0("cor.", i, ".pdf")
		df1=as.data.frame(cbind(x,y))
		p1=ggplot(df1, aes(x, y)) + 
				  xlab("Risk score") + ylab(i)+
				  geom_point() + geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
				  stat_cor(method = 'spearman', aes(x =x, y =y))
		p2=ggMarginal(p1, type="density", xparams=list(fill = "orange"), yparams=list(fill = "blue"))
		#相关性图形
		pdf(file=outFile, width=5.2, height=5)
		print(p2)
		dev.off()
	}
}

#基因与免疫相关性分析
outTab=data.frame()
risk=risk[,3:(ncol(risk)-2),drop=F]
for(immune in colnames(data)){
	for(gene in colnames(risk)){
		x=as.numeric(data[,immune])
		y=as.numeric(risk[,gene])
		corT=cor.test(x,y,method="spearman")
		cor=corT$estimate
		pvalue=corT$p.value
		text=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
		outTab=rbind(outTab,cbind(Gene=gene, Immune=immune, cor, text, pvalue))
	}
}

#绘制相关性热图
outTab$cor=as.numeric(outTab$cor)
pdf(file="geneImmuneCor.pdf", width=7, height=6)
ggplot(outTab, aes(Gene, Immune)) + 
	geom_tile(aes(fill = cor), colour = "grey", size = 1)+
	scale_fill_gradient2(low = "#5C5DAF", mid = "white", high = "#EA2E2D") + 
	geom_text(aes(label=text),col ="black",size = 3) +
	theme_minimal() + 
	theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(),
	      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),   #x轴字体
	      axis.text.y = element_text(size = 10, face = "bold")) +       #y轴字体
	labs(fill =paste0("***  p<0.001","\n", "**  p<0.01","\n", " *  p<0.05","\n", "\n","Correlation")) +   #设置图例
	scale_x_discrete(position = "bottom")      #X轴名称显示位置
dev.off()


