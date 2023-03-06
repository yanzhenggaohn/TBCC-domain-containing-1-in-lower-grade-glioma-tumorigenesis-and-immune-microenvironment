
#install.packages("meta")


library(meta)      #引用包
setwd("C://Users//刘震东//Desktop//5. Meta分析TBCCD1//已经做完")     #设置工作目录
files=dir()        #获取目录下所有文件
files=grep(".cox.txt",files,value=T)      #提取cox文件

#读取输入文件
data=data.frame()
for(i in files){
    rt=read.table(i, header=T, sep="\t", check.names=F)
    data=rbind(data, rt)
}

#提取数据信息
study=data$id
HR=data$HR
lower.HR=data$HR.95L
upper.HR=data$HR.95H

#meta分析，注意选择效应模型
meta=metagen(log(HR),
             lower=log(lower.HR),
             upper=log(upper.HR),
             studlab = study,
             sm = "HR",
             comb.random=TRUE,      #TRUE代表选择随机效应模型
             comb.fixed=FALSE)        #TRUE代表选择固定效应模型
meta

#绘制森林图
pdf(file="forest.pdf", width=12, height=5)
forest(meta,
       col.square = "green",             #方框的颜色
       col.diamond = "red", 
       col.square.lines = "green",col.diamond.lines = "red"            #菱形的颜色
       
    
       
                                       )      #菱形外框的颜色
dev.off()



