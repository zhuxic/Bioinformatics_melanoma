

#install.packages("pheatmap")


library(pheatmap)        #���ð�
expFile="uniSigGeneExp.txt"        #���������ļ�
geneCluFile="geneCluster.txt"      #������͵Ľ���ļ�
prgCluFile="PRGcluster.txt"        #ϸ���������͵Ľ���ļ�
cliFile="clinical.txt"             #�ٴ������ļ�
setwd("D:\\melanoma\\anoikis\\34.geneHeatmap")     #���ù���Ŀ¼

#��ȡ�����ļ�
exp=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
prgClu=read.table(prgCluFile, header=T, sep="\t", check.names=F, row.names=1)
geneClu=read.table(geneCluFile, header=T, sep="\t", check.names=F, row.names=1)

#�ϲ�����
exp=as.data.frame(t(exp))
sameSample=intersect(row.names(exp), row.names(prgClu))
exp=exp[sameSample,,drop=F]
expData=cbind(exp, geneCluster=geneClu[sameSample,], PRGcluster=prgClu[sameSample,])
Project=gsub("(.*?)\\_.*", "\\1", rownames(expData))
rownames(expData)=gsub("(.*?)\\_(.*?)", "\\2", rownames(expData))
expData=cbind(expData, Project)

#�ϲ��ٴ�����
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli[,"Age"]=ifelse(cli[,"Age"]=="unknow", "unknow", ifelse(cli[,"Age"]>65,">65","<=65"))
sameSample=intersect(row.names(expData), row.names(cli))
expData=expData[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
data=cbind(expData, cli)

#��ȡ��ͼ����
data=data[order(data$geneCluster),]
Type=data[,((ncol(data)-2-ncol(cli)):ncol(data))]
data=t(data[,1:(ncol(expData)-3)])

#������ɫ
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
ann_colors=list()
PRGcol=bioCol[1:length(levels(factor(Type$PRGcluster)))]
names(PRGcol)=levels(factor(Type$PRGcluster))
ann_colors[["PRGcluster"]]=PRGcol
GENEcol=bioCol[1:length(levels(factor(Type$geneCluster)))]
names(GENEcol)=levels(factor(Type$geneCluster))
ann_colors[["geneCluster"]]=GENEcol

#��ͼ���ӻ�
pdf("heatmap.pdf", height=6, width=8)
pheatmap(data,
         annotation=Type,
         annotation_colors = ann_colors,
         color = colorRampPalette(c(rep("blue",5), "white", rep("red",5)))(50),
         cluster_cols =F,
         cluster_rows =F,
         scale="row",
         show_colnames=F,
         show_rownames=F,
         fontsize=6,
         fontsize_row=8,
         fontsize_col=8)
dev.off()


