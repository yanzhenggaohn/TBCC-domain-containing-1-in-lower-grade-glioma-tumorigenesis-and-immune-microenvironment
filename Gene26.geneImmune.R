

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggpubr")
#install.packages("vioplot")
#install.packages("ggExtra")


#???ð?
library(limma)
library(reshape2)
library(ggpubr)
library(vioplot)
library(ggExtra)

expFile="geneExp.txt"              #?????????ļ?
immFile="CIBERSORT-Results.txt"    #????ϸ???????Ľ????ļ?
pFilter=0.05            #????ϸ???????????Ĺ???????
setwd("C://Users//����//Desktop//3. TBCCD1����ϸ��//5. immuneCor")     #???ù???Ŀ¼

#??ȡ?????????ļ?
rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
gene=colnames(rt)[1]

#ɾ????????Ʒ
tumorData=rt[rt$Type=="Tumor",1,drop=F]
tumorData=as.matrix(tumorData)
rownames(tumorData)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(tumorData))
data=avereps(tumorData)

#????Ŀ??????????��????Ʒ???з???
data=as.data.frame(data)
data$gene=ifelse(data[,gene]>median(data[,gene]), "High", "Low")

#??ȡ????ϸ???????ļ??????????ݽ???????
immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)
immune=immune[immune[,"P-value"]<pFilter,]
immune=as.matrix(immune[,1:(ncol(immune)-3)])

#ɾ????????Ʒ
group=sapply(strsplit(row.names(immune),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
immune=immune[group==0,]
row.names(immune)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(immune))
immune=avereps(immune)

#???ݺϲ?
sameSample=intersect(row.names(immune), row.names(data))
rt=cbind(immune[sameSample,,drop=F], data[sameSample,,drop=F])


##################????????ͼ##################
#??????ת????ggplot2?????ļ?
data=rt[,-(ncol(rt)-1)]
data=melt(data,id.vars=c("gene"))
colnames(data)=c("gene", "Immune", "Expression")
#????????ͼ
group=levels(factor(data$gene))
data$gene=factor(data$gene, levels=c("Low","High"))
bioCol=c("#0066FF","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(group)]
boxplot=ggboxplot(data, x="Immune", y="Expression", fill="gene",
				  xlab="",
				  ylab="Fraction",
				  legend.title=gene,
				  width=0.8,
				  palette=bioCol)+
				  rotate_x_text(50)+
	stat_compare_means(aes(group=gene),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "")), label="p.signif")
#????ͼƬ
pdf(file="immune.diff.pdf", width=7, height=6)
print(boxplot)
dev.off()


##########??????????ɢ??ͼ##########
outTab=data.frame()
for(i in colnames(rt)[1:(ncol(rt)-2)]){
	x=as.numeric(rt[,gene])
	y=as.numeric(rt[,i])
	if(sd(y)==0){y[1]=0.00001}
	cor=cor.test(x, y, method="spearman")
	outVector=cbind(Cell=i, cor=cor$estimate, pvalue=cor$p.value)
	outTab=rbind(outTab,outVector)
	if(cor$p.value<0.05){
		outFile=paste0("cor.", i, ".pdf")
		df1=as.data.frame(cbind(x,y))
		p1=ggplot(df1, aes(x, y)) + 
				  xlab(paste0(gene, " expression")) + ylab(i)+
				  geom_point() + geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
				  stat_cor(method = 'spearman', aes(x =x, y =y))
		p2=ggMarginal(p1, type="density", xparams=list(fill = "orange"), yparams=list(fill = "blue"))
		#??????ͼ??
		pdf(file=outFile, width=5.2, height=5)
		print(p2)
		dev.off()
	}
}
#?????????ԵĽ????ļ?
write.table(outTab,file="cor.result.txt",sep="\t",row.names=F,quote=F)


######Video source: https://ke.biowolf.cn
######??????ѧ??: https://www.biowolf.cn/
######΢?Ź??ںţ?biowolf_cn
######???????䣺biowolf@foxmail.com
######????΢??: 18520221056

