import scanpy as sc
import sys
import os
import matplotlib.pyplot as plt
import squidpy as sq
save_path = sys.argv[2]
mat = sys.argv[1]

# 定义一个保存图像的函数
def save_figure(adata, genes, save_path, filename, plot_type, dpi=300, Finally=False):
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    if plot_type == 'umap':
        fig_list = sc.pl.umap(
            adata,
            color=["my_score", "leiden"],
            frameon=False,
            ncols=3,
            size=5,
            show=False,
            return_fig=True  # Make sure it returns the figure
        )
        # Check if it's a list and loop through it to save each figure
        if isinstance(fig_list, list):
            for idx, fig in enumerate(fig_list):
                full_save_path = os.path.join(save_path, f'{filename}_{plot_type}_{idx}.pdf')
                fig.savefig(full_save_path, format='pdf', dpi=dpi)
                plt.close(fig)
        else:
            full_save_path = os.path.join(save_path, f'{filename}_{plot_type}.pdf')
            fig_list.savefig(full_save_path, format='pdf', dpi=dpi)
            plt.close(fig_list)

    elif plot_type == 'dotplot':
        # Dotplot only returns a single figure
        sc.pl.dotplot(adata, "my_score", groupby="leiden", standard_scale="var", show=False)
        fig = plt.gcf()  # Get the current figure
        full_save_path = os.path.join(save_path, f'{filename}_{plot_type}.pdf')
        fig.savefig(full_save_path, format='pdf', dpi=dpi)
        plt.close(fig)
    elif plot_type == 'spatial':
        fig = sq.pl.spatial_scatter(adata, color=genes, size=1, shape=None, edges_color="black", show=False,
                                    return_fig=True)
        full_save_path = os.path.join(save_path, f'{filename}_{plot_type}.pdf')
        fig.savefig(full_save_path, format='pdf', dpi=dpi)
        plt.close(fig)

    if Finally:
        print(save_path)



xenium_gene_ids = ["PIN2","WUS","EPFL8b"]
# 读取数据并进行预处理
gene_id_file_path = f'{mat}/geneid.txt'
sc_test = sc.read_10x_mtx(path=mat)

# 数据过滤步骤
min_genes = 100  # 每个细胞的最小基因数
min_cells = 3    # 每个基因的最小细胞数
sc.pp.filter_cells(sc_test, min_genes=min_genes)
sc.pp.filter_genes(sc_test, min_cells=min_cells)

# 创建计数层
sc_test.layers["counts"] = sc_test.X.copy()

# 数据标准化和对数转换
sc.pp.normalize_total(sc_test, inplace=True)
sc.pp.log1p(sc_test)

# 筛选高度变异基因
sc.pp.highly_variable_genes(sc_test, flavor="seurat", n_top_genes=2000)

# PCA降维
sc.pp.pca(sc_test)

# 计算邻居图
sc.pp.neighbors(sc_test, n_neighbors=10, n_pcs=30)

# 运行UMAP嵌入
sc.tl.umap(sc_test)

# 运行t-SNE嵌入
sc.tl.tsne(sc_test)

# 运行Leiden聚类
sc.tl.leiden(sc_test, resolution=0.5)


# 检查基因ID文件并生成相应的图像
if os.path.exists(gene_id_file_path):
    gene_ids = []
    with open(gene_id_file_path, 'r') as file:
        next(file)  # 跳过标题行
        for line in file:
            gene_ids.append(line.strip())
sc.tl.score_genes(sc_test, xenium_gene_ids, ctrl_size=50, n_bins=25, score_name='my_score')

save_figure(sc_test, gene_ids, save_path, f'spatial_gene_projection', plot_type='umap')
save_figure(sc_test, gene_ids, save_path, f'spatial_gene_projection', plot_type='dotplot')
save_figure(sc_test, gene_ids, save_path, f'spatial_gene_projection', plot_type='spatial', Finally=True)
