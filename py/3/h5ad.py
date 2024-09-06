
import scanpy as sc
import sys
import pandas as pd
import numpy as np
import squidpy as sq
import gzip
import matplotlib.pyplot as plt
from pandas.api.types import CategoricalDtype
import os
save_path = sys.argv[2]
mat = sys.argv[1]
# 定义保存图像的函数
def save_figure(adata, genes, save_path, filename, plot_type='umap', dpi=300):
    """
    保存UMAP或dotplot图像到指定路径

    参数:
    adata : AnnData
        单细胞数据对象
    genes : list
        要显示的基因列表或颜色列表
    save_path : str
        保存图像的路径
    filename : str
        文件名（不包括扩展名）
    plot_type : str
        图像类型 ('umap' 或 'dotplot')
    dpi : int
        图像的分辨率（每英寸点数）
    """
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    # 检查数据是否存在
    valid_genes = [gene for gene in genes if gene in adata.var_names or gene in adata.obs.columns]
    if not valid_genes:
        print(f"No valid data found for genes: {genes}. Skipping plot.")
        return

    # 生成指定类型的图形并获取Figure对象
    if plot_type == 'umap':
        fig = sc.pl.umap(adata, color=valid_genes, frameon=False, ncols=3, size=5, show=False, return_fig=True)
    elif plot_type == 'dotplot':
        fig = sc.pl.dotplot(adata, valid_genes, groupby="leiden", standard_scale="var", show=False, return_fig=True)
    else:
        print("Unsupported plot type.")
        return

    # 保存图像为PDF格式
    full_save_path = os.path.join(save_path, f'{filename}_{plot_type}.pdf')
    fig.savefig(full_save_path, format='pdf', dpi=dpi)
    plt.close(fig)
    print(save_path)

# # 设置保存路径
# save_path = "D:/Desktop/fsdownload/2/result"

# 检查文件并生成相应图像
gene_id_file_path = f'{mat}/geneid.txt'
sc_test = sc.read_h5ad(f"{mat}/sc_test.h5ad")

if os.path.exists(gene_id_file_path):
    # 如果基因ID文件存在，读取基因并生成多个基因的UMAP图
    gene_ids = []
    with open(gene_id_file_path, 'r') as file:
        next(file)  # 跳过标题行
        for line in file:
            gene_ids.append(line.strip())

    # 生成多个基因的UMAP图
    for gene_id in gene_ids:
        save_figure(sc_test, [gene_id, "leiden"], save_path, f'gene_projection_{gene_id}', plot_type='umap')
else:
    # 如果基因ID文件不存在，生成单个基因的UMAP图
    save_figure(sc_test, ["Solyc02g067340.5", "leiden"], save_path, 'single_gene_projection', plot_type='umap')

# 无论基因ID文件是否存在，都生成dotplot图
save_figure(sc_test, ["Solyc02g067340.5"], save_path, 'dotplot_gene_projection', plot_type='dotplot')
