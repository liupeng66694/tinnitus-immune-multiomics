# ==========================
# 外部验证集诊断模型的ROC曲线（增强诊断版）
# ==========================

# 加载所需包
library(glmnet)
library(pROC)

# ====== 1. 设置工作目录和文件路径 ======
setwd("D:\\02_生信分析狮\\tinnitus\\34外部验证集诊断模型的ROC曲线")
expFile <- "GSE65682.normalize.txt"
geneFile <- "gene.txt"

# ====== 2. 读取表达数据 ======
rt <- read.table(expFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

# ====== 3. 提取样本分组标签 ======
y <- gsub("(.*)\\_(.*)", "\\2", colnames(rt))
y <- ifelse(y == "Control", 0, 1)
cat("分组统计：\n")
print(table(y))
if (length(unique(y)) != 2) {
  stop("分组标签不是二分类（0/1），请检查列名中的标识是否为 Control 和 Treat")
}

# ====== 4. 读取基因列表 ======
geneRT <- read.table(geneFile, header = FALSE, sep = "\t", check.names = FALSE)
geneList <- as.vector(geneRT[, 1])

# ====== 5. 统一大小写 ======
rownames(rt) <- tolower(rownames(rt))
geneList <- tolower(geneList)

# 记录缺失的基因
missing_genes <- geneList[!geneList %in% rownames(rt)]
if (length(missing_genes) > 0) {
  warning("以下基因不在表达矩阵中，将被跳过：\n", paste(missing_genes, collapse = ", "))
  geneList <- geneList[geneList %in% rownames(rt)]
}
if (length(geneList) == 0) stop("没有有效的基因可用于分析")

# ====== 6. 绘制每个基因的ROC曲线（单基因），并记录跳过详情 ======
bioCol <- rainbow(length(geneList), s = 0.9, v = 0.9)
aucText <- c()
k <- 0
skipped_genes <- list()  # 存储跳过的基因及原因

pdf(file = "ROC.genes.pdf", width = 5, height = 4.75)

for (gene in geneList) {
  expr <- as.numeric(rt[gene, ])
  
  # 检查表达值有效性
  if (any(!is.finite(expr))) {
    reason <- "存在非有限值（NA/Inf）"
    skipped_genes[[gene]] <- reason
    warning(paste("基因", gene, reason, "跳过"))
    next
  }
  if (length(unique(expr)) < 2) {
    reason <- "表达值无变异（全相同）"
    skipped_genes[[gene]] <- reason
    warning(paste("基因", gene, reason, "跳过"))
    next
  }
  
  # 计算ROC
  roc_obj <- roc(y, expr, quiet = TRUE)
  k <- k + 1
  
  if (k == 1) {
    plot(roc_obj, print.auc = FALSE, col = bioCol[k], legacy.axes = TRUE, 
         main = "Single-gene ROC Curves")
  } else {
    plot(roc_obj, print.auc = FALSE, col = bioCol[k], legacy.axes = TRUE, 
         main = "", add = TRUE)
  }
  
  aucText <- c(aucText, paste0(gene, ", AUC = ", sprintf("%.3f", roc_obj$auc[1])))
}

# 添加图例（仅显示成功绘制的基因）
if (k > 0) {
  legend("bottomright", legend = aucText, lwd = 2, bty = "n", 
         col = bioCol[1:k], cex = 0.7)
}
dev.off()

# 输出统计信息
cat(paste("\n成功绘制", k, "个基因的ROC曲线，结果保存为 ROC.genes.pdf\n"))
cat(paste("原始基因列表中共有", length(geneRT[,1]), "个基因\n"))
cat(paste("匹配后剩余", length(geneList), "个基因（不在表达矩阵中的已剔除）\n"))
cat(paste("因数据无效而被跳过的基因数：", length(skipped_genes), "\n"))

if (length(skipped_genes) > 0) {
  cat("\n===== 被跳过的基因及原因 =====\n")
  for (gene in names(skipped_genes)) {
    cat(gene, ": ", skipped_genes[[gene]], "\n")
  }
  # 保存到文件
  skip_df <- data.frame(Gene = names(skipped_genes), Reason = unlist(skipped_genes))
  write.table(skip_df, file = "skipped_genes_detail.txt", row.names = FALSE, sep = "\t", quote = FALSE)
  cat("\n详细跳过信息已保存至 skipped_genes_detail.txt\n")
}

# ====== 7. 构建多基因逻辑回归模型（使用成功绘制的基因） ======
# 注意：这里使用成功绘制的基因（即 geneList[!geneList %in% names(skipped_genes)]）
valid_genes <- geneList[!geneList %in% names(skipped_genes)]
if (length(valid_genes) == 0) stop("没有有效的基因可用于逻辑回归建模")
rt_sub <- t(rt[valid_genes, , drop = FALSE])
rt_sub <- as.data.frame(rt_sub)

logit <- glm(y ~ ., family = binomial(link = 'logit'), data = rt_sub)
pred <- predict(logit, newdata = rt_sub, type = "response")

roc_model <- roc(y, pred, quiet = TRUE)
ci_auc <- ci.auc(roc_model, method = "bootstrap")
ci_vec <- as.numeric(ci_auc)

pdf(file = "ROC.model.pdf", width = 5, height = 4.75)
plot(roc_model, print.auc = TRUE, col = "red", legacy.axes = TRUE, 
     main = "Logistic Regression Model")
polygon(c(roc_model$specificities, 1), c(roc_model$sensitivities, 0), 
        col = rgb(1, 0, 0, 0.2), border = NA)
text(0.39, 0.43, paste0("95% CI: ", sprintf("%.3f", ci_vec[1]), "-", 
                        sprintf("%.3f", ci_vec[3])), col = "red", cex = 0.8)
dev.off()
cat("模型ROC曲线已保存为 ROC.model.pdf\n")

# ====== 8. 扰动标签验证 ======
set.seed(123)
y_noisy <- sample(y)
logit_noisy <- glm(y_noisy ~ ., family = binomial(link = 'logit'), data = rt_sub)
pred_noisy <- predict(logit_noisy, newdata = rt_sub, type = "response")
roc_noisy <- roc(y_noisy, pred_noisy, quiet = TRUE)

pdf(file = "ROC.model_noisy.pdf", width = 5, height = 4.75)
plot(roc_noisy, print.auc = TRUE, col = "blue", legacy.axes = TRUE, 
     main = "Model with Shuffled Labels (Permutation Test)")
polygon(c(roc_noisy$specificities, 1), c(roc_noisy$sensitivities, 0), 
        col = rgb(0, 0, 1, 0.2), border = NA)
dev.off()
cat("扰动标签模型的ROC曲线已保存为 ROC.model_noisy.pdf\n")

cat("\n所有分析完成！\n")




