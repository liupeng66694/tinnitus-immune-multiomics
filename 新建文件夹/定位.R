# ============================================================
# 基于 GSE167078_combined 的单细胞分析流程
# ============================================================
#BiocManager::install("SingleR")
# 安装出问题的包
#BiocManager::install("HDF5Array")
#BiocManager::install("celldex")
# 加载必要的包
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(SingleR)
library(celldex)

# 1. 读取您保存的 Seurat 对象
scRNA <- readRDS("GSE167078_combined.rds")
# 或者如果对象在环境中，直接使用

# 2. 查看基本信息
print(scRNA)
# 查看样本来源（GSM5092462 和 GSM5092463）
table(scRNA$orig.ident)

# 3. 质量控制
# 计算线粒体基因比例（小鼠线粒体基因以 mt- 开头）
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^mt-")

# 可视化质控指标
VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# 过滤细胞（根据实际情况调整阈值）
scRNA <- subset(scRNA, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)

# 4. 标准化、特征选择、缩放
scRNA <- NormalizeData(scRNA)
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
scRNA <- ScaleData(scRNA)

# 5. 线性降维（PCA）
scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA))
ElbowPlot(scRNA, ndims = 50)

# 6. 聚类和 UMAP（假设选择前 20 个主成分）
n_pcs <- 20
scRNA <- FindNeighbors(scRNA, dims = 1:n_pcs)
scRNA <- FindClusters(scRNA, resolution = 0.5)
scRNA <- RunUMAP(scRNA, dims = 1:n_pcs)

# 可视化聚类结果
DimPlot(scRNA, reduction = "umap", label = TRUE)

# 7. 细胞类型注释（小鼠数据使用 MouseRNAseqData）
# 加载包
library(Seurat)
library(SingleR)
library(celldex)

# 获取表达数据（注意使用 layer = "data"）
# 提取标准化后的表达矩阵（data 层）
scRNA_for_anno <- LayerData(scRNA, assay = "RNA", layer = "data")

# 如果上述命令报错，也可以尝试：
# scRNA_for_anno <- GetAssayData(scRNA, assay = "RNA", layer = "data")

# 加载小鼠参考数据集（新函数）
ref_data <- celldex::MouseRNAseqData()

# 运行 SingleR 注释
pred <- SingleR(test = scRNA_for_anno, ref = ref_data, labels = ref_data$label.main)

# 将注释结果添加到 Seurat 对象
scRNA$celltype <- pred$pruned.labels

# 可视化
DimPlot(scRNA, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE)

# ============================================================
# 修改后的代码：5个目标基因 + 稳健的亚群提取
# ============================================================
# 目标基因（小鼠名）
# ============================================================
# 修改为：只针对 3 个阳性因果基因（Slc9a1, Crim1, Gabbr1）
# ============================================================


scRNA_clean <- subset(scRNA, subset = !is.na(celltype))
# 如果还没有 scRNA_clean，先创建
scRNA_clean <- subset(scRNA, subset = !is.na(celltype))
cat("过滤后细胞数:", ncol(scRNA_clean), "\n")

# 1. 定义 3 个目标基因
target_genes <- c("Slc9a1", "Crim1", "Gabbr1")
existing_genes <- target_genes[target_genes %in% rownames(scRNA_clean)]
cat("存在的阳性基因:", paste(existing_genes, collapse = ", "), "\n")

# 2. 特征图（UMAP）- 可以使用原始 scRNA 或 scRNA_clean，通常 UMAP 坐标相同
FeaturePlot(scRNA_clean, features = existing_genes, ncol = 3)

# 3. 点图：按细胞类型查看表达（使用 scRNA_clean）
DotPlot(scRNA_clean, features = existing_genes, group.by = "celltype") + 
  RotatedAxis() + 
  labs(title = "Causal gene expression across cell types")

# 4. 提取免疫细胞亚群（从 scRNA_clean 中提取）
# 查看清理后的细胞类型
all_celltypes <- unique(scRNA_clean$celltype)
cat("所有细胞类型:", paste(all_celltypes, collapse = ", "), "\n")

# 定义免疫细胞关键词（根据您的细胞类型调整，例如 "Monocytes", "Macrophages"）
immune_keywords <- c("Monocyte", "Macrophage", "Microglia", "Dendritic", "DC", 
                     "Neutrophil", "Granulocyte", "Myeloid", "Inflammatory")

# 筛选包含关键词的细胞类型
immune_types <- all_celltypes[grepl(paste(immune_keywords, collapse = "|"), 
                                    all_celltypes, ignore.case = TRUE)]

if (length(immune_types) > 0) {
  cat("选中的免疫细胞类型:", paste(immune_types, collapse = ", "), "\n")
  sub_cells <- subset(scRNA_clean, subset = celltype %in% immune_types)
  cat("提取到", ncol(sub_cells), "个免疫细胞\n")
} else {
  warning("未识别到免疫细胞类型，将使用全部非NA细胞")
  sub_cells <- scRNA_clean
}

# 5. 保存对象
saveRDS(scRNA_clean, file = "GSE167078_clean.rds")
if (exists("sub_cells") && ncol(sub_cells) > 100) {
  saveRDS(sub_cells, file = "GSE167078_immune_cells.rds")
  cat("已保存免疫细胞亚群到 GSE167078_immune_cells.rds\n")
}
saveRDS(scRNA_clean, file = "GSE167078_clean.rds")
if (exists("scRNA_clean")) {
  saveRDS(scRNA_clean, file = "GSE167078_clean.rds")
  cat("文件已保存至:", file.path(getwd(), "GSE167078_clean.rds"), "\n")
} else {
  stop("对象 scRNA_clean 不存在")
}

# 提取标准化后的数据矩阵（data 层）
expr_matrix <- as.matrix(LayerData(scRNA_clean, assay = "RNA", layer = "data"))
