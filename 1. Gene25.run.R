

#install.packages('e1071')

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("preprocessCore")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


library("limma")         #???ð?
expFile="symbol.txt"     #?????????ļ?
setwd("C://Users//刘震东//Desktop//REEP4 免疫分组//25.CIBERSORT")      #???ù???Ŀ¼

#??ȡ?????ļ????????????ļ?????
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#????ת??
v=voom(data, plot=F, save.plot=F)
out=v$E
out=rbind(ID=colnames(out), out)
write.table(out,file="uniq.symbol.txt",sep="\t",quote=F,col.names=F)        #?????ļ?

#????CIBERSORT???õ?????ϸ???????Ľ???
source("Gene25.CIBERSORT.R")
results=CIBERSORT("ref.txt", "uniq.symbol.txt", perm=1000, QN=TRUE)


