# 导入必要的库
import sys
import pandas as pd
import numpy as np
import gzip
from pandas.api.types import CategoricalDtype
import scanpy as sc
import os
import matplotlib.pyplot as plt
import squidpy as sq
save_path = sys.argv[2]

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


# 处理Xenium数据
# 获取命令行参数指定的数据文件路径
mat = sys.argv[1]

# 读取10x格式的细胞特征矩阵
xenium_test = sc.read_10x_h5(filename=f"{mat}/cell_feature_matrix.h5")
# 读取细胞信息表格，使用gzip打开csv文件
df = pd.read_csv(f"{mat}/cells.csv.gz")
#Xenium数据
df.set_index(xenium_test.obs_names, inplace=True)
xenium_test.obs = df.copy()
xenium_test.obsm["spatial"] = xenium_test.obs[["x_centroid", "y_centroid"]].copy().to_numpy()
xenium_test.layers["counts"] = xenium_test.X.copy()
sc.pp.normalize_total(xenium_test, inplace=True)
sc.pp.log1p(xenium_test)
sc.pp.pca(xenium_test)
sc.pp.neighbors(xenium_test)
sc.tl.umap(xenium_test)
sc.tl.tsne(xenium_test,use_rep='X_pca')
sc.tl.leiden(xenium_test, resolution=0.5)
#输出聚类图（五张）
save_figure(xenium_test, save_path, 'clustered_data_umap', plot_type='umap')  # 无leiden
save_figure(xenium_test, save_path, 'clustered_data_tsne', plot_type='tsne')  # 无leiden
save_figure(xenium_test, save_path, 'clustered_data_umap_leiden', plot_type='umap', key='leiden')
save_figure(xenium_test, save_path, 'clustered_data_tsne_leiden', plot_type='tsne', key='leiden',Finally=True)

save_figure(xenium_test, save_path, 'clustered_data_spatial', plot_type='spatial', key='leiden',Finally=True)
#输出两表
average_expression = xenium_test.to_df().groupby(xenium_test.obs['leiden']).mean()
average_expression.to_csv(os.path.join(save_path, 'xenium_average_expression.csv'))
markers = sc.tl.rank_genes_groups(xenium_test, groupby='leiden', method='wilcoxon', key_added='rank_genes_groups')
results = sc.get.rank_genes_groups_df(xenium_test, group="7")
pd.DataFrame(results).to_csv(os.path.join(save_path, 'xenium_markers.csv'), index=False)
