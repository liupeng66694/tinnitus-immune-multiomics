# 1. 加载需要的包
library(dplyr)

# 2. 指定目录和文件模式
directory <- "D:\\02_生信分析狮\\tinnitus\\11SNP数据汇总"
pattern <- "\\SNP.csv$"

# 3. 获取目录下所有.csv文件的路径
files <- list.files(directory, pattern = pattern, full.names = TRUE)

# 4. 读取并合并所有文件
all_data <- lapply(files, read.csv) %>% bind_rows()

# 5. 保存合并后的数据表
write.csv(all_data, file.path(directory, "合并后的表格.csv"), row.names = FALSE)
