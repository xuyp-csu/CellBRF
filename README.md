# CellBRF
## Description
A Balanced Random Forest-based unsupervised feature selection algorithm for single-cell RNA-seq clustering.
CellBRF is a fast unsupervised gene selection algorithm based on balanced random forest. Importantly, CellBRF provide a two-sides balancing strategy that can be used to select gene features in scRNA-seq data with imbalanced cell types. In addition, CellBRF used a more robust ensemble learning model on balanced data to assess the importance of each gene. Finally, CellBRF combine the linear correlation of gene features to propose an efficient automatic redundant feature removal method.

![alt text](https://github.com/xuyp-csu/CellBRF/blob/main/Fig1.workflow.png?raw=True)

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

```

```
