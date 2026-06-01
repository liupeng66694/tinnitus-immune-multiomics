#devtools::install_github("mrcieu/gwasglue",force = TRUE)
#BiocManager::install("VariantAnnotation")

#install.packages("devtools")
#devtools::install_github("mrcieu/gwasglue", force = TRUE)

#install.packages("remotes")
#remotes::install_github("MRCIEU/TwoSampleMR")


#引用包
library(VariantAnnotation)
library(gwasglue)
library(TwoSampleMR)


exposureFile="eqtl-a-ENSG00000150938.vcf.gz"             #结局数据(需修改)

setwd("D:\\02_生信分析狮\\tinnitus\\36筛选强关联性SNP")     #设置工作目录



#读取结局数据的vcf文件,并对数据进行格式转换
vcfRT=readVcf(exposureFile)
exposureData=gwasvcf_to_TwoSampleMR(vcf=vcfRT, type="exposure")

write.csv(exposureData, file="exposureFile.csv")

#读取整理好的结局数据



# 保留 eaf_col 小于 5e-06 的数据行
exposure_data_filtered <- exposureData[exposureData$pval.exposure < 5e-06, ]

# 将过滤后的数据保存为 CSV 文件
write.csv(exposure_data_filtered, file = "exposure_data_filtered.csv", row.names = FALSE)


