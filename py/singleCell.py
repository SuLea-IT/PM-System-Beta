import scanpy as sc
import pandas as pd
import numpy as np
import gzip
from pandas.api.types import CategoricalDtype
mat = "/home/lilab/jupyter/测试数据/单细胞数据/12d1/filtered_feature_bc_matrix"
adata = sc.read_10x_mtx(path=mat)
print(adata)
