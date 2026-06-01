import pandas as pd
import os

out_dir = "KO_results_py"
genes = ["Crim1", "Slc9a1", "Gabbr1"]
all_genes = set()

for g in genes:
    csv_file = os.path.join(out_dir, f"{g}_diffRegulation.csv")
    if not os.path.exists(csv_file):
        print(f"警告：{csv_file} 不存在")
        continue
    df = pd.read_csv(csv_file)
    # 注意列名中可能包含空格
    sig = df[df['adjusted p-value'] < 0.05]['Gene'].tolist()
    all_genes.update(sig)
    print(f"{g}: 显著基因数 {len(sig)}")

# 保存合并后的基因列表
output_file = os.path.join(out_dir, "combined_sig_genes.txt")
with open(output_file, "w") as f:
    f.write("\n".join(sorted(all_genes)))

print(f"总显著基因数: {len(all_genes)}")
print(f"合并列表已保存至 {output_file}")