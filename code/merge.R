

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("sva")


#引用包
library(limma)
library(sva)

setwd("D:\\melanoma\\anoikis\\11.merge")     #设置工作目录
files=c("symbol.txt", "GSE19234.txt","GSE65904.txt","GSE66839.txt")      #输入文件名称

#对所有文件进行循环,获取交集基因的名称
geneList=list()
for(i in 1:length(files)){
    inputFile=files[i]
    rt=read.table(inputFile, header=T, sep="\t",check.names=F)
    header=unlist(strsplit(inputFile, "\\.|\\-"))
    header[1]=gsub("symbol", "TCGA", header[1])
    geneList[[header[1]]]=as.vector(rt[,1])
}
intersectGenes=Reduce(intersect, geneList)

#数据合并
allTab=data.frame()
batchType=c()
for(i in 1:length(files)){
    inputFile=files[i]
    header=unlist(strsplit(inputFile, "\\.|\\-"))
    header[1]=gsub("symbol", "TCGA", header[1])
    #读取输入文件，并对输入文件进行整理
    rt=read.table(inputFile, header=T, sep="\t", check.names=F)
    rt=as.matrix(rt)
    rownames(rt)=rt[,1]
    exp=rt[,2:ncol(rt)]
    dimnames=list(rownames(exp),colnames(exp))
    data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
    rt=avereps(data)
    colnames(rt)=paste0(header[1], "_", colnames(rt))
    #如果是TCGA的数据,删除正常样品
    if(header[1] == "TCGA"){
		group=sapply(strsplit(colnames(rt),"\\-"), "[", 4)
		group=sapply(strsplit(group,""), "[", 1)
		rt=rt[,group==0]
		rt=t(rt)
		row.names(rt)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(rt))
		rt=avereps(rt)
		rt=t(rt)
    }
    #判断数据是否取了log2,如果没有取log2,对数值自动取log2
    qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
    if(LogC){
    	rt[rt<0]=0
        rt=log2(rt+1)}
    if(header[1] != "TCGA"){
    	rt=normalizeBetweenArrays(rt)
    }
    #数据合并
    if(i==1){
    	allTab=rt[intersectGenes,]
    }else{
    	allTab=cbind(allTab, rt[intersectGenes,])
    }
    batchType=c(batchType, rep(i,ncol(rt)))
}

#对数据进行批次矫正，输出矫正后的表达数据
outTab=ComBat(allTab, batchType, par.prior=TRUE)
outTab=rbind(geneNames=colnames(outTab), outTab)
write.table(outTab, file="merge.txt", sep="\t", quote=F, col.names=F)

