import scanpy as sc
import squidpy as sq
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

save_path = sys.argv[2]
mat = sys.argv[1]

# 定义一个保存图像的函数
def save_figure(adata, genes, save_path, filename, plot_type='umap', dpi=300):
    """
    保存UMAP、dotplot或空间散点图到指定路径

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
        图像类型 ('umap', 'dotplot' 或 'spatial')
    dpi : int
        图像的分辨率（每英寸点数）
    """
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    # 检查数据是否有效
    valid_genes = [gene for gene in genes if gene in adata.var_names or gene in adata.obs.columns]
    if not valid_genes:
        print(f"没有找到有效的基因数据: {genes}。跳过此绘图。")
        return

    # 生成指定类型的图形
    if plot_type == 'umap':
        fig = sc.pl.umap(adata, color=valid_genes, frameon=False, ncols=3, size=5, show=False, return_fig=True)
    elif plot_type == 'dotplot':
        fig = sc.pl.dotplot(adata, valid_genes, groupby="leiden", standard_scale="var", show=False, return_fig=True)
    elif plot_type == 'spatial':
        fig = sq.pl.spatial_scatter(adata, color=valid_genes, size=1, shape=None, edges_color="black", show=False, return_fig=True)
    else:
        return

    # 保存图像为PDF格式
    full_save_path = os.path.join(save_path, f'{filename}_{plot_type}.pdf')
    fig.savefig(full_save_path, format='pdf', dpi=dpi)
    plt.close(fig)
    print(full_save_path)

# 检查基因ID文件并生成相应的图像
gene_id_file_path = f'{mat}/geneid.txt'
sc_test = sc.read_h5ad(f"{mat}/sc_test.h5ad")

if os.path.exists(gene_id_file_path):
    # 如果基因ID文件存在，读取基因并生成多个基因的UMAP图
    gene_ids = []
    with open(gene_id_file_path, 'r') as file:
        next(file)  # 跳过标题行
        for line in file:
            gene_ids.append(line.strip())

    # 为scRNA-seq数据生成多个基因的UMAP图
    for gene_id in gene_ids:
        save_figure(sc_test, [gene_id, "leiden"], save_path, f'gene_projection_{gene_id}', plot_type='umap')
        save_figure(sc_test, [gene_id], save_path, f'dotplot_gene_projection_{gene_id}', plot_type='dotplot')

    # 为空间数据生成多个基因的UMAP、dotplot和空间散点图
    for gene_id in gene_ids:
        save_figure(st_test, [gene_id, "leiden"], save_path, f'spatial_gene_projection_{gene_id}', plot_type='umap')
        save_figure(st_test, [gene_id], save_path, f'spatial_dotplot_gene_projection_{gene_id}', plot_type='dotplot')
        save_figure(st_test, [gene_id], save_path, f'spatial_scatter_gene_projection_{gene_id}', plot_type='spatial')
else:
    # 如果基因ID文件不存在，生成单个基因的UMAP图
    save_figure(sc_test, ["Solyc02g067340.5", "leiden"], save_path, 'single_gene_projection', plot_type='umap')
    save_figure(sc_test, ["Solyc02g067340.5"], save_path, 'dotplot_gene_projection', plot_type='dotplot')

    # 为空间数据生成单个基因的UMAP、dotplot和空间散点图
    save_figure(st_test, ["Solyc02g067340.5", "leiden"], save_path, 'spatial_single_gene_projection', plot_type='umap')
    save_figure(st_test, ["Solyc02g067340.5"], save_path, 'spatial_dotplot_gene_projection', plot_type='dotplot')
    save_figure(st_test, ["Solyc02g067340.5"], save_path, 'spatial_scatter_gene_projection', plot_type='spatial')
