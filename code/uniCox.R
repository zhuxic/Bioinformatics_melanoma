

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages('survival')


#引用包
library(limma)
library(survival)

expFile="interGeneExp.txt"     #差异基因的表达文件
cliFile="time.txt"             #生存数据文件
setwd("D:\\melanoma\\anoikis\\31.uniCox")     #设置工作目录

#读取表达文件，并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
data1=t(data)
rownames(data1)=gsub("(.*?)\\_(.*?)", "\\2", rownames(data1))

#读取生存数据
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)     #读取临床文件
cli$futime=cli$futime/365

#数据合并
sameSample=intersect(row.names(data1), row.names(cli))
data1=data1[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
rt=cbind(cli, data1)

#对基因进行循环，找出预后相关的基因
outTab=data.frame()
sigGenes=c()
for(i in colnames(rt[,3:ncol(rt)])){
	#cox分析
	cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
	coxSummary = summary(cox)
	coxP=coxSummary$coefficients[,"Pr(>|z|)"]
	if(coxP<0.05){
		sigGenes=c(sigGenes,i)
		outTab=rbind(outTab,
				         cbind(id=i,
				         HR=coxSummary$conf.int[,"exp(coef)"],
				         HR.95L=coxSummary$conf.int[,"lower .95"],
				         HR.95H=coxSummary$conf.int[,"upper .95"],
				         pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
				        )
	}
}

#输出单因素的结果
write.table(outTab,file="uniCox.txt",sep="\t",row.names=F,quote=F)

#保存单因素显著基因的表达量
sigGeneExp=data[sigGenes,]
sigGeneExp=rbind(id=colnames(sigGeneExp), sigGeneExp)
write.table(sigGeneExp, file="uniSigGeneExp.txt", sep="\t", quote=F, col.names=F)

#保存单因素显著基因表达和生存合并的文件
sigExpTime=rt[,c("futime", "fustat", sigGenes)]
sigExpTime=rbind(id=colnames(sigExpTime), sigExpTime)
write.table(sigExpTime, file="uniSigExpTime.txt", sep="\t", quote=F, col.names=F)



