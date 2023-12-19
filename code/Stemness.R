

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggplot2")
#install.packages("ggpubr")
#install.packages("ggExtra")


#引用包
library(limma)
library(ggplot2)
library(ggpubr)
library(ggExtra)
riskFile="risk.all.txt"       #风险文件
RNAssFile="StemnessScores_RNAexp_20170127.2.tsv"           #干细胞打分文件
setwd("D:\\melanoma\\anoikis\\52.Stemness")     #设置工作目录

#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#读取RNA干细胞打分文件
RNAss=read.table(RNAssFile, header=T, sep="\t",check.names=F, row.names=1)
RNAss=t(RNAss[1,,drop=F])
rownames(RNAss)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(RNAss))
RNAss=avereps(RNAss)

#样品取交集
sameSample=intersect(row.names(risk), row.names(RNAss))
risk=risk[sameSample,"riskScore",drop=F]
RNAss=RNAss[sameSample,,drop=F]
data=cbind(RNAss, risk)

#绘制相关性散点图
xlab="riskScore"
ylab="RNAss"
outFile="RNAss.cor.pdf"
x=as.numeric(data[,xlab])
x[x>quantile(x,0.99)]=quantile(x,0.99)
y=as.numeric(data[,ylab])
df1=as.data.frame(cbind(x,y))
p1=ggplot(df1, aes(x, y)) + 
		  xlab("Risk score") + ylab(ylab)+ ylim(0,0.7)+
		  geom_point() + geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
		  stat_cor(method = 'spearman', aes(x =x, y =y))
p2=ggMarginal(p1, type="density", xparams=list(fill = "orange"), yparams=list(fill = "blue"))
#相关性图形
pdf(file=outFile, width=5.2, height=5)
print(p2)
dev.off()



