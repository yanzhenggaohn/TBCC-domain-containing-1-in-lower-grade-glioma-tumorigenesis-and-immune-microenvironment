#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggpubr")


options(stringsAsFactors=F)
#引用包
library(limma)
library(ggpubr)

inputFile="symbol.txt"          #表达输入文件
cliFile="clinical.txt"          #临床输入文件
gene="TBCCD1"                    #基因名称
setwd("C://Users//Administrator//Desktop//TBCCD1//1.TCGA数据库//2.甲基化数据//7.临床相关性分析//mRNA的表达临床相关性分析")     #修改工作目录

#读取表达文件，并对输入文件整理
rt=read.table(inputFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#删掉正常样品
#group=sapply(strsplit(colnames(data),"\\-"),"[",4)
#group=sapply(strsplit(group,""),"[",1)
#group=gsub("2","1",group)
#data=data[,group==0]

#提取目标基因表达量
data=rbind(data,gene=data[gene,])
exp=as.matrix(t(data[c("gene",gene),]))
rownames(exp)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",rownames(exp))
exp=avereps(exp)

#读取临床数据文件
cli=read.table(cliFile,sep="\t",header=T,check.names=F,row.names=1)

#合并数据
samSample=intersect(row.names(exp),row.names(cli))
exp=exp[samSample,]
cli=cli[samSample,]
rt=cbind(exp,cli)

#临床相关性分析，输出图形结果
for(clinical in colnames(rt[,3:ncol(rt)])){
	data=rt[c(gene,clinical)]
	colnames(data)=c("gene","clinical")
	data=data[(data[,"clinical"]!="unknow"),]
	#设置比较组
	group=levels(factor(data$clinical))
	data$clinical=factor(data$clinical, levels=group)
	comp=combn(group,2)
	my_comparisons=list()
    for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
	#绘制boxplot
	boxplot=ggboxplot(data, x="clinical", y="gene", color="clinical",
	          xlab=clinical,
	          ylab=paste(gene,"expression"),
	          legend.title=clinical,
	          add = "jitter")+ 
	stat_compare_means(comparisons = my_comparisons)
	pdf(file=paste0(clinical,".pdf"),width=5.5,height=5)
	print(boxplot)
	dev.off()
}


