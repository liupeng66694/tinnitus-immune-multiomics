#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

# 如果没有安装BiocManager包，则先安装BiocManager包，方便安装生物信息学相关软件包。
#BiocManager::install("limma")  # 安装limma包，用于基因表达数据分析

library(limma)  # 加载limma包，用于基因差异表达分析
expFile = "merge.normalize.txt"  # 指定基因表达数据文件的路径
geneFile = "gene.txt"  # 指定基因列表文件的路径
setwd("D:\\02_生信分析狮\\tinnitus\\23提取基因的表达数据")  # 设置工作目录到指定路径

# 读取表达数据文件，并进行数据处理
rt = read.table(expFile, header = TRUE, sep = "\t", check.names = FALSE)  # 读取带有标题的表达数据文件，禁用列名检查
rt = as.matrix(rt)  # 将数据框转换为矩阵形式
rownames(rt) = rt[,1]  # 将第一列设置为行名（基因名）
exp = rt[, 2:ncol(rt)]  # 去除第一列，只保留表达值数据
dimnames = list(rownames(exp), colnames(exp))  # 设置矩阵的行名和列名
data = matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)  # 将数据转换为数值矩阵，确保数据为数值类型
data = avereps(data)  # 使用avereps函数去除重复基因，合并重复的基因
data = data[rowMeans(data) > 0, ]  # 过滤掉平均表达量为0或以下的基因

# 读取基因列表文件，并筛选出在表达数据中同时存在的基因
gene = read.table(geneFile, header = FALSE, sep = "\t", check.names = FALSE)  # 读取基因列表文件
sameGene = intersect(as.vector(gene[, 1]), rownames(data))  # 找到基因列表与表达数据中的共同基因
geneExp = data[sameGene, ]  # 提取共同基因的表达数据

# 将共同基因的表达数据输出到文件
outTab = rbind(ID = colnames(geneExp), geneExp)  # 将列名作为ID添加到结果的第一行
write.table(outTab, file = "snpExp.txt", sep = "\t", quote = FALSE, col.names = FALSE)  # 将结果保存为文件
