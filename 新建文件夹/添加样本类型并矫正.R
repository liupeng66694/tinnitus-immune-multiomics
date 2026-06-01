#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggpubr")
#install.packages("pROC")
# 加载所需的R包
library(limma)      # 加载limma包用于数据标准化
library(ggpubr)     # 加载ggpubr包用于绘图
library(pROC)       # 加载pROC包用于绘制ROC曲线

# 设置文件路径
expFile="geneMatrix.txt"          # 表达矩阵文件路径
conFile="sample1.txt"       # 对照组样本信息文件路径
treatFile="sample2.txt"     # 实验组样本信息文件路径
setwd("D:\\02_生信分析狮\\tinnitus\\15疾病第一个数据添加样本类型并矫正")  # 设置工作目录

# 读取表达矩阵文件
rt=read.table(expFile, header=T, sep="\t", check.names=F)  # 读取表达矩阵文件
rt=as.matrix(rt)  # 将数据转换为矩阵格式
rownames(rt)=rt[,1]  # 设置行名为基因名
exp=rt[,2:ncol(rt)]  # 去除第一列的基因名，仅保留表达数据
dimnames=list(rownames(exp), colnames(exp))  # 设置矩阵的行列名
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)  # 将表达数据转换为数值型矩阵
rt=avereps(data)  # 平均重复值

# 对表达值进行log2转换
qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))  # 计算分位数
LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )  # 判断是否需要log2转换
if(LogC){  
  rt[rt<0]=0  # 将小于0的值设为0
  rt=log2(rt+1)}  # 进行log2转换
data=normalizeBetweenArrays(rt)  # 在数组之间进行标准化

# 读取样本信息
con=read.table(conFile, header=F, sep="\t", check.names=F)  # 读取对照组样本信息
treat=read.table(treatFile, header=F, sep="\t", check.names=F)  # 读取实验组样本信息
conData=data[,as.vector(con[,1])]  # 提取对照组表达数据
treatData=data[,as.vector(treat[,1])]  # 提取实验组表达数据
data=cbind(conData, treatData)  # 合并对照组和实验组数据
conNum=ncol(conData)  # 获取对照组样本数量
treatNum=ncol(treatData)  # 获取实验组样本数量

# 给样本类型添加标签并保存标准化数据
Type=c(rep("Control",conNum), rep("Treat",treatNum))  # 生成样本类型标签
outData=rbind(id=paste0(colnames(data),"_",Type),data)  # 合并样本信息和表达数据
write.table(outData, file="GSE28750.normalize.txt", sep="\t", quote=F, col.names=F)  # 保存标准化后的数据
