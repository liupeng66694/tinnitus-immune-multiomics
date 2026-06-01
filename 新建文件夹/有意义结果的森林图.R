# 安装和加载包（如果需要安装包，请取消注释以下两行）
#install.packages("grid")
#install.packages("readr")
#install.packages("forestploter")

# 引用包
library(grid)
library(readr)
library(forestploter)

# 设置展示的方法
selectMethod <- c("Inverse variance weighted")  # 设置展示的方法

# 设置工作目录
setwd("D:\\02_生信分析狮\\tinnitus\\10森林图")  # 设置工作目录

# 获取目录下所有文件
files <- dir()  # 获取目录下所有文件
files <- grep("MRresult.csv$", files, value = TRUE)  # 提取csv结尾的文件

# 读取孟德尔随机化分析的结果
data <- data.frame()  # 初始化数据框
for (i in files) {
  rt <- read.csv(i, header = TRUE, sep = ",", check.names = FALSE)  # 读取每个csv文件
  data <- rbind(data, rt)  # 将读取的数据绑定到一起
}
data <- data[(data$method %in% selectMethod), ]  # 过滤数据，只保留指定的方法
lineVec <- cumsum(c(1, table(data[,"exposure"])))  # 计算每个暴露变量的累积和

# 对数据进行整理
data$' ' <- paste(rep(" ", 10), collapse = " ")  # 添加一个空白列
data$'OR(95% CI)' <- ifelse(is.na(data$or), "", sprintf("%.3f (%.3f to %.3f)", data$or, data$or_lci95, data$or_uci95))  # 计算OR和95%CI
data$pval <- ifelse(data$pval < 0.001, "<0.001", sprintf("%.3f", data$pval))  # 格式化p值
data$exposure <- ifelse(is.na(data$exposure), "", data$exposure)  # 处理暴露变量的缺失值
data$nsnp <- ifelse(is.na(data$nsnp), "", data$nsnp)  # 处理SNP数目的缺失值

# 仅在存在重复暴露变量名时进行替换
if (any(duplicated(data$exposure))) {
  data[duplicated(data$exposure), ]$exposure <- ""  # 去除重复的暴露变量名
}

# 准备图形参数
tm <- forest_theme(base_size = 18,  # 图形整体的大小
                   ci_pch = 16, ci_lty = 1, ci_lwd = 1.5, ci_col = "black", ci_Theight = 0.2,  # 可信区间的形状、线条类型、宽度、颜色、两端竖线高度
                   refline_lty = "dashed", refline_lwd = 1, refline_col = "grey20",  # 参考线条的形状、宽度、颜色
                   xaxis_cex = 0.8,  # x轴刻度字体的大小
                   footnote_cex = 0.6, footnote_col = "blue")  # 脚注大小、颜色

# 绘制图形
plot <- forestploter::forest(data[, c("exposure", "nsnp", "method", "pval", " ", "OR(95% CI)")],
                             est = data$or,  # OR值
                             lower = data$or_lci95,  # OR的下界
                             upper = data$or_uci95,  # OR的上界
                             ci_column = 5,  # 可信区间所在的列
                             ref_line = 1,  # 参考线条的位置
                             xlim = c(0, 3),  # X轴的范围
                             theme = tm)  # 图形的参数

# 修改图形中可信区间的颜色
boxcolor <- c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F","#8491B4","#91D1C2","#DC0000","#7E6148")  # 设置颜色
boxcolor <- boxcolor[as.numeric(as.factor(data$method))]  # 根据方法分配颜色
for (i in 1:nrow(data)) {
  plot <- edit_plot(plot, col = 5, row = i, which = "ci", gp = gpar(fill = boxcolor[i], fontsize = 25))  # 修改可信区间的颜色
}

# 设置pvalue的字体
pos_bold_pval <- which(as.numeric(gsub('<', "", data$pval)) < 0.05)  # 找到p值小于0.05的行
if (length(pos_bold_pval) > 0) {
  for (i in pos_bold_pval) {
    plot <- edit_plot(plot, col = 4, row = i, which = "text", gp = gpar(fontface = "bold"))  # 设置p值的字体为粗体
  }
}

# 在图形中增加线段
plot <- add_border(plot, part = "header", row = 1, where = "top", gp = gpar(lwd = 2))  # 添加顶部边框
plot <- add_border(plot, part = "header", row = lineVec, gp = gpar(lwd = 1))  # 添加每个暴露变量的边框

# 设置字体大小，并且将文字居中
plot <- edit_plot(plot, col = 1:ncol(data), row = 1:nrow(data), which = "text", gp = gpar(fontsize = 12))  # 设置字体大小
plot <- edit_plot(plot, col = 1:ncol(data), which = "text", hjust = unit(0.5, "npc"), part = "header", x = unit(0.5, "npc"))  # 标题文字居中
plot <- edit_plot(plot, col = 1:ncol(data), which = "text", hjust = unit(0.5, "npc"), x = unit(0.5, "npc"))  # 内容文字居中

# 输出图形
pdf("forest.pdf", width = 20, height = 25)  # 创建PDF文件
print(plot)  # 输出图形到PDF
dev.off()  # 关闭PDF设备

### 结束 ###
