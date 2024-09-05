import sys
import scanpy as sc
import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# 从命令行参数获取路径
mat = sys.argv[1]
save_path = sys.argv[2]

# 定义保存图片的函数
def save_figure(adata, save_path, filename, plot_type='umap', key=None, dpi=300,Finally=False):
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    if plot_type == 'umap':
        sc.pl.umap(adata, color=key,show=False)
    elif plot_type == 'tsne':
        sc.pl.tsne(adata, color=key,show=False)
    elif plot_type == 'spatial':
        sq.pl.spatial_scatter(adata, size=1, shape=None, edges_color="black", color=key,show=False)
    else:
        return

    fig = plt.gcf()  # 获取当前图形
    full_save_path = os.path.join(save_path, f'{filename}_{plot_type}.pdf')
    fig.savefig(full_save_path, format='pdf', dpi=dpi)
    plt.close(fig)
    if Finally:
        print(save_path)

# 读取数据
adata = sc.read_10x_mtx(path=mat)

# 过滤数据
min_genes = 100  # 每个细胞的最小基因数
min_cells = 3    # 每个基因的最小细胞数
sc.pp.filter_cells(adata, min_genes=min_genes)
sc.pp.filter_genes(adata, min_cells=min_cells)

# 检查是否已经计算了UMAP和t-SNE
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
    sc.tl.umap(adata)
    sc.tl.tsne(adata)
    sc.tl.leiden(adata, resolution=0.5)
elif has_umap and not has_tsne:
    sc.tl.tsne(adata)
elif not has_umap and has_tsne:
    sc.tl.umap(adata)

# 保存聚类图
save_figure(adata, save_path, 'clustered_data_umap', plot_type='umap')  # 无leiden
save_figure(adata, save_path, 'clustered_data_tsne', plot_type='tsne')  # 无leiden
save_figure(adata, save_path, 'clustered_data_umap_leiden', plot_type='umap', key='leiden')
save_figure(adata, save_path, 'clustered_data_tsne_leiden', plot_type='tsne', key='leiden',Finally=True)

# 输出平均表达量表
average_expression = adata.to_df().groupby(adata.obs['leiden']).mean()
average_expression.to_csv(os.path.join(save_path, 'sc_average_expression.csv'))

# 计算和保存差异表达基因
sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon', key_added='rank_genes_groups')
results = sc.get.rank_genes_groups_df(adata, group="7")
results.to_csv(os.path.join(save_path, 'st_markers.csv'), index=False)

# 保存聚类后的数据为h5ad格式
# adata.write(os.path.join(save_path, 'clustered_data.h5ad'))

