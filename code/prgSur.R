
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages('survival')
#install.packages("survminer")


#���ð�
library(limma)
library(survival)
library(survminer)

expFile="prgGeneExp.txt"     #���������ļ�
cliFile="time.txt"           #���������ļ�
setwd("D:\\melanoma\\anoikis\\20.prgSur")     #���ù���Ŀ¼

#��ȡ�����ļ������������ļ�����
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
data=t(data)
rownames(data)=gsub("(.*?)\\_(.*?)", "\\2", rownames(data))

#��ȡ��������
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli$futime=cli$futime/365

#���ݺϲ�
sameSample=intersect(row.names(data), row.names(cli))
data=data[sameSample,]
cli=cli[sameSample,]
rt=cbind(cli, data)

#�Ի������ѭ�����ҳ�Ԥ����صĻ���
outTab=data.frame()
km=c()
for(i in colnames(rt[,3:ncol(rt)])){
	#cox����
	cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
	coxSummary = summary(cox)
	coxP=coxSummary$coefficients[,"Pr(>|z|)"]
	outTab=rbind(outTab,
			         cbind(id=i,
			         HR=coxSummary$conf.int[,"exp(coef)"],
			         HR.95L=coxSummary$conf.int[,"lower .95"],
			         HR.95H=coxSummary$conf.int[,"upper .95"],
			         pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
			        )
	#km����
	data=rt[,c("futime", "fustat", i)]
	colnames(data)=c("futime", "fustat", "gene")
	#��ȡ����cutoff
	res.cut=surv_cutpoint(data, time = "futime", event = "fustat", variables =c("gene"))
	res.cat=surv_categorize(res.cut)
	fit=survfit(Surv(futime, fustat) ~gene, data = res.cat)
	#print(paste0(i, " ", res.cut$cutpoint[1]))
	#�Ƚϸߵͱ����������
	diff=survdiff(Surv(futime, fustat) ~gene,data =res.cat)
	pValue=1-pchisq(diff$chisq, df=1)
	km=c(km, pValue)
	#��pvalue<0.05�Ļ��������������
	if(pValue<0.05){
		if(pValue<0.001){
			pValue="p<0.001"
		}else{
			pValue=paste0("p=",sprintf("%.03f",pValue))
		}
		
		#������������
		surPlot=ggsurvplot(fit,
						   data=res.cat,
						   pval=pValue,
						   pval.size=6,
						   legend.title=i,
						   legend.labs=c("high","low"),
						   xlab="Time(years)",
						   ylab="Overall survival",
						   palette=c("red", "blue"),
						   break.time.by=1,
						   conf.int=F,
						   risk.table=F,
						   risk.table.title="",
						   risk.table.height=.25)
		pdf(file=paste0("sur.", i, ".pdf"),onefile = FALSE,
			width = 5,         #ͼƬ�Ŀ���
			height =4.5)         #ͼƬ�ĸ߶�
		print(surPlot)
		dev.off()
	}
}

#��������صĽ��
outTab=cbind(outTab, km)
write.table(outTab,file="uniCox.txt",sep="\t",row.names=F,quote=F)


