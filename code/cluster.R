

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("ConsensusClusterPlus")


library(ConsensusClusterPlus)       #���ð�
expFile="prgGeneExp.txt"            #���������ļ�
workDir="D:\\melanoma\\anoikis\\22.prgCluster"     #����Ŀ¼
setwd(workDir)       #���ù���Ŀ¼

#��ȡ�����ļ�
data=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
data=as.matrix(data)

#����Ʒ���о���
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


#������ͽ��
clusterNum=3        #�ֳɼ�������
cluster=results[[clusterNum]][["consensusClass"]]
cluster=as.data.frame(cluster)
colnames(cluster)=c("PRGcluster")
letter=c("A","B","C","D","E","F","G")
uniqClu=levels(factor(cluster$PRGcluster))
cluster$PRGcluster=letter[match(cluster$PRGcluster, uniqClu)]
clusterOut=rbind(ID=colnames(cluster), cluster)
write.table(clusterOut, file="PRGcluster.txt", sep="\t", quote=F, col.names=F)

