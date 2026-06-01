# 加载Seurat包
library(Seurat)

# 设置主目录（包含两个子文件夹的路径）
main_dir <- "."  # 当前目录，如果文件在其他路径请修改

# 获取子文件夹名称（GSM5092462, GSM5092463）
sample_dirs <- list.dirs(main_dir, full.names = TRUE, recursive = FALSE)
sample_dirs <- sample_dirs[basename(sample_dirs) %in% c("GSM5092462", "GSM5092463")]

# 用于存储每个样本的Seurat对象
seurat_list <- list()

for (dir in sample_dirs) {
  # 列出该文件夹下的所有文件
  files <- list.files(dir, full.names = TRUE)
  
  # 找到 barcodes, features, matrix 文件（匹配可能的命名模式）
  barcode_file <- files[grepl("barcodes", files, ignore.case = TRUE)]
  feature_file <- files[grepl("features|genes", files, ignore.case = TRUE)]
  matrix_file <- files[grepl("matrix", files, ignore.case = TRUE)]
  
  # 检查文件是否齐全
  if (length(barcode_file) == 0 || length(feature_file) == 0 || length(matrix_file) == 0) {
    warning(paste("在", dir, "中未找到完整的10X文件，跳过该样本"))
    next
  }
  
  # 读取数据（使用ReadMtx直接指定文件路径，无需重命名）
  counts <- ReadMtx(
    mtx = matrix_file,
    cells = barcode_file,
    features = feature_file
  )
  
  # 创建Seurat对象
  seurat_obj <- CreateSeuratObject(counts = counts, project = basename(dir))
  seurat_list[[basename(dir)]] <- seurat_obj
}

# 检查是否成功读取了至少一个样本
if (length(seurat_list) == 0) {
  stop("未成功读取任何样本，请检查文件结构。")
}

# 合并样本（如果只有一个样本则不需要合并）
if (length(seurat_list) == 1) {
  combined_seurat <- seurat_list[[1]]
} else {
  combined_seurat <- merge(seurat_list[[1]], y = seurat_list[-1])
}

# 查看合并后的对象
print(combined_seurat)

# 可选：保存为RDS文件供后续分析
saveRDS(combined_seurat, file = "GSE167078_combined.rds")