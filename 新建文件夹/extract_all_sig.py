import pandas as pd
import os

genes = ["Crim1", "Slc9a1", "Gabbr1"]
out_dir = "KO_results_py"

for gene in genes:
    df = pd.read_csv(f"{out_dir}/{gene}_diffRegulation.csv")
    sig = df[df['adjusted p-value'] < 0.05]['Gene'].tolist()
    with open(f"{out_dir}/{gene}_all_sig_genes.txt", "w") as f:
        f.write("\n".join(sig))
    print(f"{gene}: 显著基因数 {len(sig)}")