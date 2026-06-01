# 安装e1071包（如果需要用于支持向量机等算法）
#install.packages('e1071')

# 安装并加载BiocManager包，如果需要安装preprocessCore库
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("preprocessCore")

# 设置输入文件路径
inputFile = "merge.normalize.txt"  # 已标准化的合并基因表达文件

# 设置工作目录，确保文件路径正确
setwd("D:\\02_生信分析狮\\tinnitus\\27免疫浸润分析")  # 设置工作路径

# 加载包含CIBERSORT算法的自定义R脚本
source("GEOimmune.CIBERSORT.R")  # 通过CIBERSORT进行免疫细胞组成分析

# 使用CIBERSORT算法进行免疫细胞分布分析
# "ref.txt" 是参考基因表达矩阵，用于解卷积免疫细胞成分
# inputFile 是输入的基因表达数据
# perm = 1000 表示置换次数，QN=TRUE 表示是否进行量化标准化（Quantile Normalization）
outTab = CIBERSORT("ref.txt", inputFile, perm = 1000, QN = TRUE)

# 筛选P-value小于1的结果
outTab = outTab[outTab[,"P-value"] < 1,]  # 保留p值小于1的有效样本

# 将数据转化为矩阵形式，去除最后三列无关数据
outTab = as.matrix(outTab[, 1:(ncol(outTab) - 3)])  # 去掉最后三列（通常是统计信息列）

# 添加列名并将结果写入文件
outTab = rbind(id = colnames(outTab), outTab)  # 在结果中添加列名
write.table(outTab, file = "CIBERSORT-Results.txt", sep = "\t", quote = FALSE, col.names = FALSE)  # 将结果保存为txt文件