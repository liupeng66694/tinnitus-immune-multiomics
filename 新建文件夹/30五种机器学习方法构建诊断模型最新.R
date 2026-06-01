# 加载所需的R包

install.packages("caret")
install.packages("DALEX")
install.packages("ggplot2")
install.packages("ranger")
install.packages("e1071")
install.packages("pROC")
install.packages("xgboost")

library(caret)
library(DALEX)
library(ggplot2)
library(ranger)
library(e1071)
library(pROC)
library(xgboost)

set.seed(123)

inputFile = "merge.normalize.txt"
geneFile = "gene.txt"

setwd("D:\\02_生信分析狮\\tinnitus\\30四种机器学习方法构建诊断模型")

data = read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
geneRT = read.table(geneFile, header = FALSE, sep = "\t", check.names = FALSE)
data = data[as.vector(geneRT[, 1]), ]
row.names(data) = gsub("-", "_", row.names(data))

data = t(data)
group = gsub("(.*)\\_(.*)", "\\2", row.names(data))
data = as.data.frame(data)
data$Type = factor(group)
data = data[!is.na(data$Type), ]
data$Type = relevel(data$Type, ref = "Control")

inTrain <- createDataPartition(y = data$Type, p = 0.7, list = FALSE)
train <- data[inTrain, ]
test <- data[-inTrain, ]

control = trainControl(method = "repeatedcv", number = 5, repeats = 3, classProbs = TRUE, savePredictions = TRUE)
preProcess = c("center", "scale")

mod_rf = train(Type ~ ., data = train, method = 'ranger', trControl = control, preProcess = preProcess, tuneLength = 5)
mod_svm = train(Type ~ ., data = train, method = "svmLinear", trControl = control, preProcess = preProcess, tuneLength = 5)
mod_lasso = train(Type ~ ., data = train, method = "glmnet", trControl = control, preProcess = preProcess, tuneLength = 5)
mod_knn = train(Type ~ ., data = train, method = "knn", trControl = control, preProcess = preProcess, tuneLength = 5)

nb_model <- naiveBayes(Type ~ ., data = train)
mod_nb <- nb_model

nb_pred_prob <- function(object, newdata) {
  prob <- e1071:::predict.naiveBayes(object, newdata, type = "raw")
  data.frame(no = prob[, 1], yes = prob[, 2])
}

cat("\n=== 模型评估结果 ===\n")
cat("\n随机森林 (RF):\n")
print(mod_rf$results)
cat("\n支持向量机 (SVM):\n")
print(mod_svm$results)
cat("\nLASSO 模型 (LASSO):\n")
print(mod_lasso$results)
cat("\nKNN 模型 (KNN):\n")
print(mod_lasso$results)
cat("\n朴素贝叶斯 (NB):\n")
nb_pred <- nb_pred_prob(mod_nb, test)
nb_acc <- sum(apply(nb_pred, 1, which.max) == as.numeric(test$Type)) / nrow(test)
print(data.frame(Accuracy = nb_acc))

p_fun = function(object, newdata) {
  if (inherits(object, "naiveBayes")) {
    nb_pred_prob(object, newdata)[, 2]
  } else {
    predict(object, newdata = newdata, type = "prob")[, 2]
  }
}
# 确保测试集标签与训练集一致
test$Type = factor(test$Type, levels = levels(train$Type))

# 生成数值型 y（0/1），用于模型性能评估
yTest = as.numeric(test$Type) - 1
# 准备不含目标变量的特征数据集
test_features <- test[, !names(test) %in% "Type"]

# 创建 explainer 对象（使用 test_features）
explainer_rf   <- explain(mod_rf,   label = "RF",   data = test_features, y = yTest, predict_function = p_fun, verbose = FALSE)
explainer_svm  <- explain(mod_svm,  label = "SVM",  data = test_features, y = yTest, predict_function = p_fun, verbose = FALSE)
explainer_lasso<- explain(mod_lasso,label = "LASSO",data = test_features, y = yTest, predict_function = p_fun, verbose = FALSE)
explainer_knn  <- explain(mod_knn,  label = "KNN",  data = test_features, y = yTest, predict_function = p_fun, verbose = FALSE)
explainer_nb   <- explain(mod_nb,   label = "NB",   data = test_features, y = yTest, predict_function = p_fun, verbose = FALSE)

mp_rf = model_performance(explainer_rf)
mp_svm = model_performance(explainer_svm)
mp_lasso = model_performance(explainer_lasso)
mp_knn = model_performance(explainer_knn)
mp_nb = model_performance(explainer_nb)

pdf(file = "residual.pdf", width = 6, height = 6)
p1 <- plot(mp_rf, mp_svm, mp_lasso, mp_knn, mp_nb)
print(p1)
dev.off()

pdf(file = "boxplot.pdf", width = 6, height = 6)
p2 <- plot(mp_rf, mp_svm, mp_lasso, mp_knn, mp_nb, geom = "boxplot")
print(p2)
dev.off()

pred1 = predict(mod_rf, newdata = test, type = "prob")
pred2 = predict(mod_svm, newdata = test, type = "prob")
pred3 = predict(mod_lasso, newdata = test, type = "prob")
pred4 = predict(mod_knn, newdata = test, type = "prob")
pred5 = nb_pred_prob(mod_nb, test)
roc1 = roc(yTest, as.numeric(pred1[, 2]))
roc2 = roc(yTest, as.numeric(pred2[, 2]))
roc3 = roc(yTest, as.numeric(pred3[, 2]))
roc4 = roc(yTest, as.numeric(pred4[, 2]))
roc5 = roc(yTest, as.numeric(pred5[, 2]))

pdf(file = "ROC.pdf", width = 6, height = 6)
plot(roc1, print.auc = FALSE, legacy.axes = TRUE, main = "", col = "red")
plot(roc2, print.auc = FALSE, legacy.axes = TRUE, main = "", col = "blue", add = TRUE)
plot(roc3, print.auc = FALSE, legacy.axes = TRUE, main = "", col = "green", add = TRUE)
plot(roc4, print.auc = FALSE, legacy.axes = TRUE, main = "", col = "orange", add = TRUE)
plot(roc5, print.auc = FALSE, legacy.axes = TRUE, main = "", col = "purple", add = TRUE)
legend('bottomright',
       c(paste0('RF: ', sprintf("%.03f", roc1$auc)),
         paste0('SVM: ', sprintf("%.03f", roc2$auc)),
         paste0('LASSO: ', sprintf("%.03f", roc3$auc)),
         paste0('KNN: ', sprintf("%.03f", roc4$auc)),
         paste0('NB: ', sprintf("%.03f", roc5$auc))),
       col = c("red", "blue", "green", "orange","purple"), lwd = 2, bty = 'n')
dev.off()

importance_rf <- variable_importance(explainer_rf, loss_function = loss_root_mean_square)
importance_svm <- variable_importance(explainer_svm, loss_function = loss_root_mean_square)
importance_lasso <- variable_importance(explainer_lasso, loss_function = loss_root_mean_square)
importance_knn <- variable_importance(explainer_knn, loss_function = loss_root_mean_square)
importance_nb <- variable_importance(explainer_nb, loss_function = loss_root_mean_square)

pdf(file = "importance.pdf", width = 7, height = 12)
plot(importance_rf[c(1, (ncol(data) - 8):(ncol(data) + 1)), ],
     importance_svm[c(1, (ncol(data) - 8):(ncol(data) + 1)), ],
     importance_lasso[c(1, (ncol(data) - 8):(ncol(data) + 1)), ],
     importance_knn[c(1, (ncol(data) - 8):(ncol(data) + 1)), ],
     importance_nb[c(1, (ncol(data) - 8):(ncol(data) + 1)), ])
dev.off()

# ===================== 稳定的变量重要性（10 次 permutation 平均） =====================
set.seed(123)  # 确保可重复

# 定义函数：计算模型重要性（B=10, type="raw" 获取每次结果，然后聚合）
get_stable_importance <- function(explainer, n_perm = 10) {
  # 计算原始重要性（返回每个变量多次 permutation 的结果）
  imp_raw <- variable_importance(explainer, 
                                 loss_function = loss_root_mean_square,
                                 type = "raw",      # 返回每次 permutation 的原始结果
                                 B = n_perm,        # 重复次数
                                 n_sample = NULL)   # 使用全部数据
  # 过滤掉 _baseline_, _full_model_ 等特殊行
  imp_raw <- imp_raw[!grepl("^_", imp_raw$variable), ]
  # 按变量名聚合，计算 dropout_loss 的平均值
  imp_agg <- aggregate(dropout_loss ~ variable, data = imp_raw, FUN = mean)
  # 按平均损失降序排列
  imp_agg <- imp_agg[order(imp_agg$dropout_loss, decreasing = TRUE), ]
  return(imp_agg)
}

# 计算每个模型的重要性（10 次 permutation 平均）
importance_rf   <- get_stable_importance(explainer_rf,   n_perm = 10)
importance_svm  <- get_stable_importance(explainer_svm,  n_perm = 10)
importance_lasso<- get_stable_importance(explainer_lasso, n_perm = 10)
importance_knn  <- get_stable_importance(explainer_knn,  n_perm = 10)
importance_nb   <- get_stable_importance(explainer_nb,   n_perm = 10)

# 输出每个模型的前 geneNum 个基因（例如 5 个）
geneNum <- 5
write.table(head(importance_rf, geneNum),   file = "importanceGene.RF.txt",   sep = "\t", quote = FALSE, row.names = FALSE)
write.table(head(importance_svm, geneNum),  file = "importanceGene.SVM.txt",  sep = "\t", quote = FALSE, row.names = FALSE)
write.table(head(importance_lasso, geneNum),file = "importanceGene.LASSO.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(head(importance_knn, geneNum),  file = "importanceGene.KNN.txt",  sep = "\t", quote = FALSE, row.names = FALSE)
write.table(head(importance_nb, geneNum),   file = "importanceGene.NB.txt",   sep = "\t", quote = FALSE, row.names = FALSE)

# 可选：同时输出包含所有基因的完整重要性表（用于检查）
write.table(importance_rf,   file = "importance_all_RF.txt",   sep = "\t", quote = FALSE, row.names = FALSE)
write.table(importance_svm,  file = "importance_all_SVM.txt",  sep = "\t", quote = FALSE, row.names = FALSE)
write.table(importance_lasso,file = "importance_all_LASSO.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(importance_knn,  file = "importance_all_KNN.txt",  sep = "\t", quote = FALSE, row.names = FALSE)
write.table(importance_nb,   file = "importance_all_NB.txt",   sep = "\t", quote = FALSE, row.names = FALSE)