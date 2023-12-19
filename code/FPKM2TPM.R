

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


library(limma)             #引用包
inputFile="symbol.txt"     #输入文件
setwd("D:\\melanoma\\anoikis\\07.FPKM2TPM")     #设置工作目录

#读取输入文件,并对输入文件整理
outTab=data.frame()
rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#FPKM转换为TPM
fpkmToTpm=function(fpkm){
    exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
tpm=apply(data, 2, fpkmToTpm)

#输出转换结果
tpmOut=rbind(ID=colnames(tpm), tpm)
write.table(tpmOut, file="TCGA.TPM.txt", sep="\t", col.names=F, quote=F)


