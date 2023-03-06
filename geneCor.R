

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggpubr")
#install.packages("ggExtra")


#引用包
library(limma)
library(reshape2)
library(ggpubr)
library(ggExtra)

gene="TBCCD1"                    #基因名称
expFile="symbol.txt"      #基因表达文件
geneFile="gene.txt"            #基因列表文件
setwd("C://Users//刘震东//Desktop//3. TBCCD1免疫细胞//7. checkpoint")     #设置工作目录

#读取表达数据文件,并对数据进行整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

#删除正常样品
#group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
#group=sapply(strsplit(group,""), "[", 1)
#group=gsub("2", "1", group)
#data=data[,group==0]

#提取基因列表文件
geneRT=read.table(geneFile, header=T, sep="\t", check.names=F)

#相关性检验
outTab=data.frame()
for(k in 1:nrow(geneRT)){
	i=geneRT[k,2]      #获取免疫细胞的marker基因
	if(i %in% row.names(data)){
		if(sd(data[i,])>0.01){
			#相关性分析
			x=as.numeric(data[gene,])
			y=as.numeric(data[i,])
			corT=cor.test(x, y, method="spearma")
			cor=corT$estimate
			pvalue=corT$p.value
			outTab=rbind(outTab,cbind(immuneCell=geneRT[k,1], gene=i, cor, pvalue))
				
			#绘制相关性散点图
			df1=as.data.frame(cbind(x,y))
			p1=ggplot(df1, aes(x, y)) + 
					  xlab(paste0(gene, " expression")) + 
					  ylab(paste0(i, " expression")) +
					  geom_point() + geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
					  stat_cor(method = 'spearman', aes(x =x, y =y))
			p2=ggMarginal(p1, type="density", xparams=list(fill = "orange"), yparams=list(fill = "blue"))
			
			#保存相关性图形
			pdf(file=paste0("cor.",gene,"_",i,".pdf"), width=5.2, height=5)
			print(p2)
			dev.off()
		}
	}
}

#输出相关性结果
write.table(file="geneCor.result.txt",outTab,sep="\t",quote=F,row.names=F)


