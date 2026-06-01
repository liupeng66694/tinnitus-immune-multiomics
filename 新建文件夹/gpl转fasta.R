
setwd("D:\\02_生信分析狮\\tinnitus\\13疾病第一个数据下载（GEO数据库）")  # 设置工作目录
# 读取文件所有行
lines <- readLines("gpl.txt")

# 打开输出文件连接
out_conn <- file("gpl.fasta", open = "w")

# 创建进度条，进度范围从0到总行数减1（跳过第一行）
total_lines <- length(lines) - 1
pb <- txtProgressBar(min = 0, max = total_lines, style = 3)

# 遍历从第二行开始的每一行
for(i in 2:length(lines)) {
  # 按制表符分割每一行
  fields <- strsplit(lines[i], "\t")[[1]]
  
  # 将第一个字段前加上 ">" 写入输出文件（FASTA格式）
  cat(">", fields[1], "\n", sep = "", file = out_conn)
  
  # 将第二个字段写入输出文件
  cat(fields[2], "\n", sep = "", file = out_conn)
  
  # 更新进度条，进度以处理的行数计
  setTxtProgressBar(pb, i - 1)
}

# 关闭输出文件连接及进度条
close(out_conn)
close(pb)
