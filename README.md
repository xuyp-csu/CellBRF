# CellBRF
## Description
A Balanced Random Forest-based unsupervised feature selection algorithm for single-cell RNA-seq clustering.
CellBRF is a fast unsupervised gene selection algorithm based on balanced random forest. Importantly, CellBRF provide a two-sides balancing strategy that can be used to select gene features in scRNA-seq data with imbalanced cell types. In addition, CellBRF used a more robust ensemble learning model on balanced data to assess the importance of each gene. Finally, CellBRF combine the linear correlation of gene features to propose an efficient automatic redundant feature removal method.

![alt text](https://github.com/xuyp-csu/CellBRF/blob/main/Fig1.workflow.png?raw=True)

# Usage
## Requirements:
Python --- 3.6.8
pytorch -- 1.5.1+cu101
Scanpy --- 1.0.4
Nvidia Tesla P100

## Files:
CF.py -- implementation of CellBRF algorithm

```

```
