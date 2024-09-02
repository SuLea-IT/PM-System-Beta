import sys
import scanpy as sc
import pandas as pd
import numpy as np
import gzip
from pandas.api.types import CategoricalDtype
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
#百迈克空间转录组数据
mat = sys.argv[1]
coord = f"{mat}/barcodes_pos.tsv.gz"
adata = sc.read_10x_mtx(path=mat)
pos = pd.read_csv(coord, sep="\t", names=["cellID", "row", "col"])
pos2 = pos[pos["cellID"].isin(adata.to_df().index.to_list())]
cat_size_order = CategoricalDtype(adata.to_df().index.to_list(), ordered=False)
pos2['cellID'] = pos2['cellID'].astype(cat_size_order)
pos2.set_index('cellID', inplace=True)
cells_in_adata = set(adata.obs_names)
cells_in_pos2 = set(pos2.index)
missing_cells = cells_in_adata - cells_in_pos2
adata = adata[~adata.obs_names.isin(missing_cells)]
adata.obsm['spatial'] = np.array(pos2)
#过滤数据
min_genes = 100  # 每个细胞的最小基因数
min_cells = 3 # 每个基因的最小细胞数
sc.pp.filter_cells(adata, min_genes=min_genes)
sc.pp.filter_genes(adata, min_cells=min_cells)
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
save_figure(adata, 'leiden', save_path, 'clustered_data', plot_type='umap')
save_figure(adata, 'leiden', save_path, 'clustered_data', plot_type='tsne')
