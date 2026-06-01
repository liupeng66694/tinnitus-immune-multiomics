import scanpy as sc
import pandas as pd
import scipy.io
from scTenifold import scTenifoldKnk
import os

# ==================== 配置 ====================
# 如果数据文件不在当前目录，请修改路径
mtx_file = "expression_matrix.mtx"
genes_file = "genes.tsv"
cells_file = "cells.tsv"
metadata_file = "metadata.csv"
out_dir = "KO_results_py"

# 目标基因列表
target_genes = ["Crim1", "Slc9a1", "Gabbr1"]

# ==================== 读取数据 ====================
print("读取数据...")
matrix = scipy.io.mmread(mtx_file).T.tocsc()
genes = pd.read_csv(genes_file, header=None)[0].values
cells = pd.read_csv(cells_file, header=None)[0].values
metadata = pd.read_csv(metadata_file, index_col=0)

# 创建 AnnData 对象
adata = sc.AnnData(X=matrix, obs=metadata, var=pd.DataFrame(index=genes))
adata.obs_names = cells
adata.var_names = genes

print(f"AnnData 维度: {adata.shape}")

# ==================== 可选：过滤低表达基因 ====================
sc.pp.filter_genes(adata, min_cells=3)
print(f"过滤低表达基因后维度: {adata.shape}")

# ==================== 可选：高变基因筛选（减轻内存） ====================
# 如果内存不足，可以取消注释以下三行，将基因数降到 3000 左右
# sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor='seurat_v3')
# adata = adata[:, adata.var.highly_variable]
# print(f"高变基因筛选后维度: {adata.shape}")

# ==================== 创建输出目录 ====================
os.makedirs(out_dir, exist_ok=True)

# ==================== 循环处理每个基因 ====================
for gene in target_genes:
    print(f"\n{'='*40}")
    print(f"开始处理基因: {gene}")
    print(f"{'='*40}")
    
    try:
        # 初始化并运行虚拟敲除
        knk = scTenifoldKnk(data=adata, ko_genes=[gene])
        knk.build()
        
        # 调试：列出 knk 对象的所有非私有属性（便于查看可用字段）
        print(f"可用属性: {[attr for attr in dir(knk) if not attr.startswith('_')]}")
        
        # 尝试获取差异调控结果（不同版本属性名可能不同）
        diff_df = None
        if hasattr(knk, 'diff_regulation'):
            diff_df = knk.diff_regulation
            print("使用属性 diff_regulation")
        elif hasattr(knk, 'd_regulation'):
            diff_df = knk.d_regulation
            print("使用属性 d_regulation")
        else:
            print("警告：未找到差异调控结果属性，跳过保存。")
            continue
        
        # 保存完整差异调控表
        csv_path = os.path.join(out_dir, f"{gene}_diffRegulation.csv")
        diff_df.to_csv(csv_path, index=False)
        print(f"差异调控表已保存: {csv_path}")
        
        # 提取上下调基因列表（log2FC 阈值 ±1，adj.p < 0.05）
        # 注意列名可能是 'FC' 或 'log2FoldChange'，根据实际情况调整
        fc_col = 'FC' if 'FC' in diff_df.columns else 'log2FoldChange'
        p_col = 'p.adj' if 'p.adj' in diff_df.columns else 'padj'
        
        up_genes = diff_df[(diff_df[fc_col] > 1) & (diff_df[p_col] < 0.05)]['gene'].tolist()
        down_genes = diff_df[(diff_df[fc_col] < -1) & (diff_df[p_col] < 0.05)]['gene'].tolist()
        
        # 保存上下调列表
        with open(os.path.join(out_dir, f"{gene}_up_genes.txt"), "w") as f:
            f.write("\n".join(up_genes))
        with open(os.path.join(out_dir, f"{gene}_down_genes.txt"), "w") as f:
            f.write("\n".join(down_genes))
        
        print(f"{gene} 完成: 上调基因数 = {len(up_genes)}, 下调基因数 = {len(down_genes)}")
    
    except Exception as e:
        print(f"处理基因 {gene} 时出错: {e}")
        continue

print("\n所有基因处理完毕。结果保存在", out_dir)