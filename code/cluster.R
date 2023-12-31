

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("ConsensusClusterPlus")


library(ConsensusClusterPlus)       #引用包
expFile="prgGeneExp.txt"            #表达数据文件
workDir="D:\\melanoma\\anoikis\\22.prgCluster"     #工作目录
setwd(workDir)       #设置工作目录

#读取输入文件
data=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
data=as.matrix(data)

#对样品进行聚类
maxK=9
results=ConsensusClusterPlus(data,
              maxK=maxK,
              reps=50,
              pItem=0.8,
              pFeature=1,
              title=workDir,
              clusterAlg="km",
              distance="euclidean",
              seed=123456,
              plot="png")


#输出分型结果
clusterNum=3        #分成几个亚型
cluster=results[[clusterNum]][["consensusClass"]]
cluster=as.data.frame(cluster)
colnames(cluster)=c("PRGcluster")
letter=c("A","B","C","D","E","F","G")
uniqClu=levels(factor(cluster$PRGcluster))
cluster$PRGcluster=letter[match(cluster$PRGcluster, uniqClu)]
clusterOut=rbind(ID=colnames(cluster), cluster)
write.table(clusterOut, file="PRGcluster.txt", sep="\t", quote=F, col.names=F)


