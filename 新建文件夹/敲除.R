# 如果 BiocManager 还未安装，请先取消注释下面的行
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# 从 Bioconductor 安装 scTenifoldKnk
BiocManager::install("scTenifoldKnk")

# 加载包
library(scTenifoldKnk)
# 请确保您的环境中已经存在 scRNA_clean 对象
# 如果没有，请先运行之前的代码加载它

# 检查 scRNA_clean 是否存在
if (!exists("scRNA_clean")) {
  stop("未找到 scRNA_clean 对象，请先运行前面的代码")
}

# 提取表达矩阵（Seurat v5 兼容）
expr_matrix <- as.matrix(LayerData(scRNA_clean, assay = "RNA", layer = "data"))

# 确认矩阵维度
cat("表达矩阵维度:", dim(expr_matrix), "\n")
cat("矩阵中零值比例:", mean(expr_matrix == 0), "\n")

# 加载或安装必要的包
if (!require("scTenifoldKnk", quietly = TRUE)) {
  # 从GitHub安装最新版
  if (!require("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("cailab-tamu/scTenifoldKnk")
}
if (!require("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!require("tidyr", quietly = TRUE)) install.packages("tidyr")

# 加载包
library(scTenifoldKnk)
library(dplyr)
library(ggplot2)
library(tidyr)

# 从之前保存的文件加载 scRNA_clean 对象
if (file.exists("GSE167078_clean.rds")) {
  scRNA_clean <- readRDS("GSE167078_clean.rds")
  message("成功加载 scRNA_clean 对象")
} else {
  stop("未找到 GSE167078_clean.rds 文件，请先运行之前的单细胞分析流程保存该对象")
}

# 提取表达矩阵（Seurat v5 兼容）
expr_matrix <- as.matrix(LayerData(scRNA_clean, assay = "RNA", layer = "data"))

# 定义3个目标基因
target_genes <- c("Slc9a1", "Crim1", "Gabbr1")

# 保留数据中存在的基因
existing_genes <- target_genes[target_genes %in% rownames(expr_matrix)]

if (length(existing_genes) == 0) {
  stop("未找到任何目标基因！请检查基因名是否与数据中匹配。")
} else {
  message("找到的目标基因: ", paste(existing_genes, collapse = ", "))
}

# 创建结果输出目录
if (!dir.exists("KO_results")) dir.create("KO_results")

# 1. 加载包
library(scTenifoldKnk)

# 2. 检查表达矩阵 (假设你已经通过 LayerData 创建了 expr_matrix)
print(dim(expr_matrix)) # 应输出基因数和细胞数

# 3. 确认目标基因
target_gene <- "Crim1"
if (!target_gene %in% rownames(expr_matrix)) {
  stop(paste("目标基因", target_gene, "不在表达矩阵中。请确认基因名。"))
} else {
  message(paste("找到目标基因:", target_gene))
}
# 运行虚拟敲除 (修正后的代码)
ko_result <- scTenifoldKnk(
  countMatrix = expr_matrix,
  gKO = target_gene
)

# 检查输出结构
names(ko_result)
head(ko_result$diffRegulation)
