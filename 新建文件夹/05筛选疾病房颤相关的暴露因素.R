pvalFilter = 0.05     # 设置p值过滤阈值(0.05, 0.01, 0.001)
mrFile = "table.MRresult.csv"        # MR结果文件名
pleFile = "table.pleiotropy.csv"     # 多效性结果文件名
setwd("D:\\02_生信分析狮\\tinnitus\\08筛选有意义的暴露因素")  # 设置工作目录

# 读取MR结果文件
rt = read.csv(mrFile, header = TRUE, sep = ",", check.names = FALSE)

# 筛选p值符合阈值的逆方差加权（IVW）法的MR结果
ivwRT = rt[rt$method == "Inverse variance weighted", ]
ivwRT = ivwRT[ivwRT$pval < pvalFilter, ]

# 筛选OR值均大于1或均小于1的暴露变量
ivw = data.frame()
for (cell in unique(ivwRT$exposure)) {
  cellData = rt[rt$exposure == cell, ]
  if (sum(cellData$or > 1) == nrow(cellData) | sum(cellData$or < 1) == nrow(cellData)) {
    ivw = rbind(ivw, ivwRT[ivwRT$exposure == cell, ])
  }
}

# 读取多效性结果文件
pleRT = read.csv(pleFile, header = TRUE, sep = ",", check.names = FALSE)
# 筛选p值大于0.05的结果
pleRT = pleRT[pleRT$pval > 0.05, ]
cellLists = as.vector(pleRT$exposure)
outTab = ivw[ivw$exposure %in% cellLists, ]
write.csv(outTab, file = "immune-disease.IVWfilter.csv", row.names = FALSE)
