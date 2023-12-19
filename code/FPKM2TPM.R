

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


library(limma)             #���ð�
inputFile="symbol.txt"     #�����ļ�
setwd("D:\\melanoma\\anoikis\\07.FPKM2TPM")     #���ù���Ŀ¼

#��ȡ�����ļ�,���������ļ�����
outTab=data.frame()
rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#FPKMת��ΪTPM
fpkmToTpm=function(fpkm){
    exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
tpm=apply(data, 2, fpkmToTpm)

#���ת�����
tpmOut=rbind(ID=colnames(tpm), tpm)
write.table(tpmOut, file="TCGA.TPM.txt", sep="\t", col.names=F, quote=F)

