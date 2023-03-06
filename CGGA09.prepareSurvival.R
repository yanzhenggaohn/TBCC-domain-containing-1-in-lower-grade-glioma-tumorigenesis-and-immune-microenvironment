setwd("C://Users//刘震东//Desktop//TBCCD1//1.TCGA数据库//1.表达数据//09.prepareSurvival")                         #修改工作目录
expFile="rocSigExp.txt"                                                              #表达数据文件
clinicalFile="clinicalNum.txt"                                                       #临床数据文件
gene="TBCCD1"

exp=read.table(expFile,sep="\t",header=T,check.names=F,row.names=1)                #读取表达数据文件
cli=read.table(clinicalFile,sep="\t",header=T,check.names=F,row.names=1)           #读取临床数据文件
samSample=intersect(row.names(exp),row.names(cli))
exp=exp[samSample,]
cli=cli[samSample,]
selectCol=c("futime","fustat",gene)
outTab=cbind(exp[,selectCol],cli)
outTab=cbind(id=row.names(outTab),outTab)
write.table(outTab,file="singleGeneSurData.txt",sep="\t",row.names=F,quote=F)

