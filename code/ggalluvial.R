

#install.packages("ggplot2")
#install.packages("ggalluvial")


#引用包
library(ggalluvial)
library(ggplot2)
library(dplyr)

prgCluFile="PRGcluster.txt"        #细胞焦亡分型的结果文件
geneCluFile="geneCluster.txt"      #基因分型的结果文件
riskFile="risk.all.txt"            #风险文件
setwd("D:\\melanoma\\anoikis\\37.ggalluvial")     #设置工作目录

#读取输入文件
prgClu=read.table(prgCluFile, header=T, sep="\t", check.names=F, row.names=1)
geneClu=read.table(geneCluFile, header=T, sep="\t", check.names=F, row.names=1)
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#合并数据
twoCluster=cbind(prgClu, geneClu)
rownames(twoCluster)=gsub("(.*?)\\_(.*?)", "\\2", rownames(twoCluster))
sameSample=intersect(row.names(twoCluster), row.names(risk))
rt=cbind(risk[sameSample,,drop=F], twoCluster[sameSample,,drop=F])

#准备桑基图输入文件
rt=rt[,c("PRGcluster", "geneCluster", "risk", "fustat")]
colnames(rt)=c("PRGcluster", "geneCluster", "Risk", "Fustat")
rt[,"Fustat"]=ifelse(rt[,"Fustat"]==0, "Alive", "Dead")
corLodes=to_lodes_form(rt, axes = 1:ncol(rt), id = "Cohort")

#得到输出文件
pdf(file="ggalluvial.pdf", width=6, height=5.5)
mycol=rep(c("#0066FF","#FF9900","#FF0000","#029149","#6E568C","#E0367A","#D8D155","#223D6C","#D20A13","#431A3D","#91612D","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#64495D","#7CC767"),15)
ggplot(corLodes, aes(x = x, stratum = stratum, alluvium = Cohort,fill = stratum, label = stratum)) +
  	 scale_x_discrete(expand = c(0, 0)) +  
  	 #用aes.flow控制线条颜色，forward说明颜色和前面的柱状图一致，backward说明和后面的柱状图一致。
  	 geom_flow(width = 2/10,aes.flow = "forward") + 
	 geom_stratum(alpha = .9,width = 2/10) +
	 scale_fill_manual(values = mycol) +
	 #size=3代表字体大小
	 geom_text(stat = "stratum", size = 4,color="black") +
	 xlab("") + ylab("") + theme_bw() + 
	 theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text.y = element_blank()) + #去掉坐标轴
	 theme(panel.grid =element_blank()) + 
	 theme(panel.border = element_blank()) + 
	 ggtitle("") + guides(fill = FALSE)                            
dev.off()


