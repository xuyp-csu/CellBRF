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

n_clusters : Number of the cell types in the scRNA-seq dataset.

balanced : Whether the data needs to be balanced. (default: True)

RR : Whether redundant features need to be removed. (default: True)

corr_threshold : Linear correlation threshold to use when removing redundant genes. (default: 0.8)

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
