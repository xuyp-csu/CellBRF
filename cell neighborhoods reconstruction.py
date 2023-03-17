
import numpy as np
import os
import scanpy as sc
from sklearn.metrics.cluster import adjusted_rand_score
import pandas as pd
import random
import math
from sklearn.neighbors import NearestNeighbors

dataList = ["Buettner","Chu_celltime", "Chung", "Darmanis", "Deng", "Engel", "Goolam", "Kim", "Koh",
            "Kolodz", "Kumar", "Leng", "Li", "Maria2", "Pollen", "Robert", "Ting", "Treutlein",
            "Usoskin", "Yan", "Yeo", "Zhou"]
medlist = ["CellBRF", "DUBStepR","Feats","FEAST","geneBasisR", "HVG", "HRG"]
seed = 2022
random.seed(seed)
k=16
for dataName in dataList:
    dir = './'
    print(dataName)
    counts_df = pd.read_csv(os.path.join(dir, dataName + '.csv'))
    lab_df = pd.read_csv(os.path.join(dir, dataName + '_label.csv'))
    var_names = counts_df.iloc[:, 0].values.tolist()
    ad = ~pd.Series(var_names).duplicated(keep='first').values
    var_names = list(np.array(var_names)[ad])
    lab = lab_df.iloc[:, 1].values.tolist()
    types = np.unique(lab)
    y = np.zeros(len(lab))
    c = 1
    for t in types:
        idx = [i for i, j in enumerate(lab) if j == t]
        y[idx] = c
        c += 1
    X = counts_df.iloc[ad, 1:].values
    n_clusters = len(np.unique(y))

    adata = sc.AnnData(X.T)
    adata.var_names = var_names
    sc.pp.normalize_per_cell(adata)
    sc.pp.log1p(adata)

    auc_mat = pd.DataFrame(index=list(range(1, k)), columns=medlist)
    for m in medlist:
        print(m)
        if m == "CellBRF":
            gs = np.loadtxt('./' + dataName + '_CellBRF_gs_res.txt', dtype=str)
            norm_X = adata[:, gs].X.copy()
        else:
            gs = np.loadtxt('./' + m + '_' + dataName + '.txt', dtype=str)
            if list(set(gs) & set(var_names)) == []:
                tmp = [i.replace('_', '-') for i in var_names]
                gs = list(set(gs) & set(tmp))
                adata.var_names = tmp
                norm_X = adata[:, gs].X.copy()
                adata.var_names = var_names
            else:
                gs = list(set(gs) & set(var_names))
                norm_X = adata[:, gs].X.copy()
        print(norm_X.shape)
        nn = NearestNeighbors(n_neighbors=k).fit(norm_X)
        indices = nn.kneighbors(norm_X, return_distance=False)

        auc_list = []
        for i in range(1, k):
            auc = 0
            for j in range(norm_X.shape[0]):
                for z in range(1, (i + 1)):
                    if y[indices[j, z]] == y[indices[j, 0]]:
                        auc += 1
            auc_list.append(auc / (norm_X.shape[0] * i))
        auc_mat[m] = auc_list
    print(auc_mat)
    np.savetxt('./neighbor_auc/' + dataName + '.txt', auc_mat, fmt='%f')

