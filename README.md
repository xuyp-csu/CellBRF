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

data : Gene expression data matrix, gene in columns and samples in rows.

dataName : Name of scRNA-seq dataset.

geneName : The length must be the same as the number of columns in the data matrix.

n_clusters : Number of the cell types in the scRNA-seq dataset.

normalization : Whether the data needs to be normalized. (default: True)

label_predict : The length must be the same as the number of rows in the data matrix.

seed : Random seed.

RR : Whether redundant features need to be removed. (default: True)

corr_threshold : Linear correlation threshold to use when removing redundant genes. (default: 0.8)

save_full : Whether the full result needs to be saved. (default: True)

save_path : Path to save results.

## Files:
h5data -- a example scRNAseq data Han

CF.py -- implementation of CellBRF algorithm

Chu_DEgenes.Rdata -- genes with different possible expression paths in time course data Chu1

Feats_demo.py -- implementation of Feats algorithm

PanglaoDB_markers_27_Mar_2020.tsv -- cell-type makers from the PanglaoDB database

cell neighborhoods reconstruction.py -- The python code used to calculate k-nearest neighbor consistency and silhouette coefficient

evalcluster.R -- implementation of NMI and ARI

test_code.R -- R code for all other analyzes and visualizations

workflow.png -- CellBRF workflow

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

Xu Y, Li H D, Lin C X, et al. CellBRF: a feature selection method for single-cell clustering using cell balance and random forest[J]. Bioinformatics, 2023, 39(Supplement_1): i368-i376.
