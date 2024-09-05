import sys
import scanpy as sc

# 从命令行参数获取h5ad文件路径
h5ad_path = sys.argv[1]

# 读取h5ad文件
adata = sc.read_h5ad(h5ad_path)

# 获取数据量信息
num_cells = adata.n_obs  # 细胞数量
num_genes = adata.n_vars  # 基因数量

print(f"细胞数量: {num_cells}")
print(f"基因数量: {num_genes}")

# 获取前20个基因ID
gene_ids = adata.var_names[:20]
print("前20个基因ID:")
print(gene_ids)
