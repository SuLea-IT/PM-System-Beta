import scanpy as sc
import sys
import os
import matplotlib.pyplot as plt
import pandas as pd
import squidpy as sq
save_path = sys.argv[2]
mat = sys.argv[1]

# 定义一个保存图像的函数
def save_figure(adata, genes, save_path, filename, plot_type, dpi=300, Finally=False):
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    fig = None
    if plot_type == 'umap':
        fig = sc.pl.umap(
            adata,
            color=[genes, "leiden"],
            frameon=False,
            ncols=3,
            size=5, show=False, return_fig=True
        )
    elif plot_type == 'dotplot':
        sc.pl.dotplot(xenium_test, [genes] if isinstance(genes, str) else genes, groupby="leiden", standard_scale="var", show=False)
        fig = plt.gcf()  # 获取当前的 Matplotlib 图像
    elif plot_type == 'spatial':
        sq.pl.spatial_scatter(adata,library_id="spatial", color=genes, size=2, shape=None, edges_color="black",
                              img=False,
                              show=False, return_fig=True)
    else:
        return

    if fig is not None:
        # 保存图像并关闭
        full_save_path = os.path.join(save_path, f'{filename}_{plot_type}.pdf')
        fig.savefig(full_save_path, format='pdf', dpi=dpi)
        plt.close(fig)

    if Finally:
        print(save_path)




# 读取数据并进行预处理
gene_id_file_path = f'{mat}/geneid.txt'
xenium_test = sc.read_10x_h5(filename=f"{mat}/cell_feature_matrix.h5")
# 读取细胞信息表格，使用gzip打开csv文件
df = pd.read_csv(f"{mat}/cells.csv.gz")
# 数据过滤步骤
# 移除总表达量为0的细胞，避免t-SNE过程中出现问题
# sc.pp.filter_cells(xenium_test, min_genes=1)
# 创建计数层
xenium_test.layers["counts"] = xenium_test.X.copy()

# 数据标准化和对数转换
sc.pp.normalize_total(xenium_test, inplace=True)
sc.pp.log1p(xenium_test)

# 筛选高度变异基因
sc.pp.highly_variable_genes(xenium_test, flavor="seurat", n_top_genes=2000)

# PCA降维
sc.pp.pca(xenium_test)

# 计算邻居图
sc.pp.neighbors(xenium_test, n_neighbors=10, n_pcs=30)

# 运行UMAP嵌入
sc.tl.umap(xenium_test)

# 运行t-SNE嵌入
# 运行t-SNE嵌入，使用随机初始化避免PCA和稀疏矩阵冲突
sc.tl.tsne(xenium_test,use_rep='X_pca')

# 运行Leiden聚类
sc.tl.leiden(xenium_test, resolution=0.5)


# 检查基因ID文件并生成相应的图像
if os.path.exists(gene_id_file_path):
    gene_ids = []
    with open(gene_id_file_path, 'r') as file:
        next(file)  # 跳过标题行
        for line in file:
            gene_ids.append(line.strip())

    # 为空间数据生成多个基因的UMAP、dotplot和空间散点图
    for idx, gene_id in enumerate(gene_ids):
        is_last = idx == len(gene_ids) - 1  # 检查是否是最后一个基因ID
        save_figure(xenium_test, gene_id, save_path, f'spatial_gene_projection_{gene_id}', plot_type='umap')
        save_figure(xenium_test, gene_id, save_path, f'spatial_dotplot_gene_projection_{gene_id}', plot_type='dotplot',Finally=is_last)
