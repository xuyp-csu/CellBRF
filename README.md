# CellBRF
## Description
CellBRF: A feature selection method for single-cell RNA-seq clustering using cell balance and random forest
CellBRF is a feature selection method that considers genes’ relevance to cell types for single-cell clustering. The key idea is to identify genes with more clustering contributions through random forests guided by predicted cell labels. Moreover, it proposes a class balancing strategy to mitigate the impact of skewed cell types distributions on feature importance measurement.

![alt text](https://github.com/xuyp-csu/CellBRF/blob/main/workflow.png)

# Usage
## Requirements:

Python --- 3.7.13

Scanpy --- 1.9.1

scipy --- 1.7.3

scikit-learn --- 1.0.2

imblearn --- 0.0

## Arguments:

data : `pandas.DataFrame` or `2-D numpy.array`, optional
        Gene expression data matrix, gene in columns and samples in rows.

dataName : string
    Name of scRNA-seq dataset.

geneName : list -> string
    The length must be the same as the number of columns in the data matrix.
    Names of all genes in the scRNA-seq dataset.

n_clusters : integer
    Number of the cell types in the scRNA-seq dataset.

normalization : boolean
    Whether the data needs to be normalized. (default: True)

label_predict : list -> integer, optional
    The length must be the same as the number of rows in the data matrix.
    Predicted labels per cell. (default: None)

seed : integer
    Random seed.

RR : boolean
    Whether redundant features need to be removed. (default: True)

corr_threshold : integer
    Linear correlation threshold to use when removing redundant genes. (default: 0.8)

save_full : boolean
    Whether the full result needs to be saved. (default: True)

save_path : string
    Path to save results.

## Files:
CF.py -- implementation of CellBRF algorithm

## Demo:
```
import h5py
import numpy as np
import CF
import os

dir = './h5data'
dataName = 'Cao'
data_mat = h5py.File(os.path.join(dir, dataName + '.h5'))
X = np.array(data_mat['X'])
y = np.array(data_mat['Y'])
data_mat.close()
n_clusters = len(np.unique(y))
X_selected = CF.CellBRF(data=X, dataName=dataName, n_clusters=n_clusters, save_path='./')
```

# Contact
If any questions, please do not hesitate to contact us at: 

Yunpei Xu, xu_yunpei@csu.edu.cn

Jianxin Wang, jxwang@csu.edu.cn

# How to cite?
If you use this tool, please cite the following work.

Yunpei Xu, Hongdong Li, Cuixiang Lin, Ruiqing Zheng, Fangxiang Wu and Jianxin Wang, CellBRF: A feature selection method for single-cell RNA-seq clustering using cell balance and random forest, 2022, submitted  
