# 导入必要的库
import sys
import pandas as pd
import numpy as np
import gzip
from pandas.api.types import CategoricalDtype
import scanpy as sc
import os
import matplotlib.pyplot as plt

save_path = sys.argv[2]

def save_figure(adata, key, save_path, filename, plot_type='umap', dpi=300):
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    if plot_type == 'umap':
        sc.pl.umap(adata, color=key, return_fig=True)
    elif plot_type == 'tsne':
        sc.pl.tsne(adata, color=key, return_fig=True)
    else:
        return
    fig = plt.gcf()
    full_save_path = os.path.join(save_path, f'{filename}_{plot_type}.pdf')
    fig.savefig(full_save_path, format='pdf', dpi=dpi)
    plt.close()
    print(save_path)  # 假设这是生成的文件

# 处理Xenium数据
# 获取命令行参数指定的数据文件路径
mat = sys.argv[1]

# 读取10x格式的细胞特征矩阵
adata = sc.read_10x_h5(filename=f"{mat}/cell_feature_matrix.h5")
# 读取细胞信息表格，使用gzip打开csv文件
df = pd.read_csv(f"{mat}/cells.csv.gz")

# 将文件中的索引列设置为adata的观察名称
df.set_index(adata.obs_names, inplace=True)

# 将读取的细胞信息复制到adata的观察数据中
adata.obs = df.copy()

# 将细胞的中心坐标数据从观察数据复制到观察矩阵数据中
adata.obsm["spatial"] = adata.obs[["x_centroid", "y_centroid"]].copy().to_numpy()

# 移除具有零计数的细胞
sc.pp.filter_cells(adata, min_counts=1)

# 检查是否已有UMAP和TSNE结果
has_umap = 'X_umap' in adata.obsm
has_tsne = 'X_tsne' in adata.obsm

if not has_umap and not has_tsne:
    # 数据预处理步骤
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)

    # 将稀疏矩阵转换为密集矩阵
    adata.X = adata.X.toarray()

    sc.tl.umap(adata)
    sc.tl.tsne(adata)  # 移除 'init' 参数
    sc.tl.leiden(adata, resolution=0.5)
elif has_umap and not has_tsne:
    adata.X = adata.X.toarray()  # 确保数据为密集矩阵
    sc.tl.tsne(adata)
elif not has_umap and has_tsne:
    sc.tl.umap(adata)

save_figure(adata, 'leiden', save_path, 'clustered_data', plot_type='umap')
save_figure(adata, 'leiden', save_path, 'clustered_data', plot_type='tsne')
