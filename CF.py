"""
MIT License

Copyright (c) [year] [fullname]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import numpy as np
import pandas as pd
import scanpy as sc
import random
from sklearn.cluster import SpectralClustering
from sklearn.ensemble import RandomForestClassifier
from collections import Counter
from sklearn.neighbors import kneighbors_graph
from sklearn.decomposition import PCA
import math
import time
from imblearn.over_sampling import SMOTE
from imblearn.under_sampling import NearMiss
from imblearn import FunctionSampler

# def predict_k(affinity_matrix):
#
#     return k

def dopca(X, dim=10):
    pcaten = PCA(n_components=dim)
    X_10 = pcaten.fit_transform(X)
    return X_10


def get_adj(count, k=15, pca=50, mode="connectivity"):
    if pca:
        countp = dopca(count, dim=pca)
    else:
        countp = count
    A = kneighbors_graph(countp, k, mode=mode, metric="euclidean", include_self=True)
    adj = A.toarray()
    return adj


def predict_label(data, k, seed):
    adj = get_adj(count=data)
    # if k is None:
    #     k = predict_k(adj)
    label_predict = SpectralClustering(n_clusters=k,
                                       affinity="precomputed",
                                       assign_labels="discretize",
                                       random_state=seed).fit_predict(adj)
    return label_predict


def freq_cal(l, thd):
    counter = Counter(l)
    gene_idx = []
    for i in dict(counter):
        if dict(counter)[i] >= thd:
            gene_idx.extend([i])
    return gene_idx


def filter_bycorr_with_orderly_genes(df,
                                    threshold=0.8,
                                    gap=0.1):

    def _get_diff_list(a_column, a_list, removed=None):
        all_cols = a_list
        if removed is not None and len(removed) > 0:
            all_cols = [aa for aa in a_list if aa not in removed]
        return [aa for aa in all_cols if aa != a_column]

    columns = df.columns.tolist()
    removed = []
    cal_ed = []
    gene_cor = df.corr()
    for cc in columns:
        cal_ed.append(cc)
        if cc not in removed:
            tmp_cols = _get_diff_list(cc,
                                      columns,
                                      removed=removed + cal_ed)
            thred_diff = gap * 1.0 / (len(tmp_cols) + 1)
            count = len(tmp_cols)
            for tt in tmp_cols:
                count -= 1
#                 relation = df[cc].corr(df[tt])
                relation = gene_cor.loc[cc, tt]
                if abs(relation) > threshold + thred_diff * count:
                    removed.append(tt)
    return [cc for cc in columns if cc not in removed]


def USSampler(x, y, random_state):
    y = y.astype(int)
    encode = dict(zip(np.unique(y), range(len(np.unique(y)))))
    tmp_y = []
    for i in y:
        tmp_y.append(encode[i])

    size = math.ceil(x.shape[0]/len(np.unique(tmp_y)))
    resize = []
    for c in np.unique(tmp_y):
        idx = list(np.where(tmp_y == c)[0])
        if len(idx) < size:
            resize.append(len(idx))
        else:
            resize.append(size)
    ratio = dict(zip(np.unique(tmp_y), resize))
    nm_us = NearMiss(sampling_strategy=ratio, version=1)
    reshape_x, tmp_reshape_y = nm_us.fit_resample(x, tmp_y)

    decode = dict(zip(range(len(np.unique(y))), np.unique(y)))
    reshape_y = []
    for i in tmp_reshape_y:
        reshape_y.append(decode[i])
    return reshape_x, reshape_y

"""
CellBRF: A Balanced Random Forest-based unsupervised feature selection algorithm for single-cell RNA-seq clustering.

This function uses a fast two-sides data balancing strategy and random forest model to identify informative genes.
"""
def CellBRF(data,
            dataName=None,
            geneNames=None,
            n_clusters=None,
            normalization=True,
            label_predict=None,
            balanced=True,
            seed=2022,
            RR=True,
            corr_threshold=0.8,
            save_full=True,
            n_features=None,
            save_path='./'):
    """
    Parameters
    ----------
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

    balanced : boolean
        Whether the data needs to be balanced. (default: True)

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

    Returns
    ----------
    X_selected : `numpy.ndarray`
        Gene expression matrix after feature selection.
    """

    random.seed(seed)
    # 0. normalization
    if normalization:
        adata = sc.AnnData(data)
        if geneNames is not None:
            adata.var_names = geneNames
        else:
            adata.var_names = [str(i) for i in range(adata.X.shape[1])]
        sc.pp.filter_genes(adata, min_cells=3)
        sc.pp.normalize_per_cell(adata)
        sc.pp.log1p(adata)
        norm_X = adata.X.copy()
        gene_names = adata.var_names
    else:
        norm_X = data.copy()
        if geneNames is None:
            gene_names = [str(i) for i in range(norm_X.shape[1])]
    n_cells = norm_X.shape[0]

    # 1. Label prediction
    if label_predict is None:
        print(">> Spectral clustering to generate predicted labels......")
        label_predict = predict_label(data=norm_X, k=n_clusters, seed=seed)
        predict = True
    else:
        predict = False

    start = time.time()
    # 2. Feature Importance Assessment
    counter = Counter(label_predict)
    cdtb = [freq for c, freq in counter.items()]
    print('>> data size: ', norm_X.shape)
    print('>> clusters distribution: ', cdtb)

    # 3. Data Balancing
    if balanced:
        balanced_thd = math.ceil(n_cells / n_clusters)
        main_clusters = []
        rare_clusters = []
        diff = []
        clusters = []
        for i in counter.items():
            if (i[1] > balanced_thd):
                main_clusters.append(i[0])
            else:
                rare_clusters.append(i[0])
            diff.append(i[1] - balanced_thd)
            clusters.append(i[0])
        central_cluster = clusters[np.where(np.abs(diff) == min(np.abs(diff)))[0][0]]
        if central_cluster in main_clusters:
            main_clusters.remove(central_cluster)
        else:
            rare_clusters.remove(central_cluster)

        if len(main_clusters) == 0 or len(rare_clusters) == 0:
            tmp_X = norm_X.copy()
            tmp_lab = label_predict.copy()
        else:
            main_idx = []
            for i in main_clusters:
                idx = list(np.where(label_predict == i)[0])
                main_idx.extend(idx)
            rare_idx = []
            min_cluster = 10
            for i in rare_clusters:
                idx = list(np.where(label_predict == i)[0])
                rare_idx.extend(idx)
                if len(idx)<min_cluster:
                    min_cluster = len(idx)
            central_idx = list(np.where(label_predict == central_cluster)[0])

            os_X = norm_X[central_idx + rare_idx, :].copy()
            os_lab = label_predict[central_idx + rare_idx].copy()
            if min_cluster<7:
                smo_os = SMOTE(random_state=seed, k_neighbors=min_cluster-1)
            else:
                smo_os = SMOTE(random_state=seed)
            os_X, os_lab = smo_os.fit_resample(os_X, os_lab)

            us_X = norm_X[central_idx + main_idx, :].copy()
            us_lab = label_predict[central_idx + main_idx].copy()
            fs_us = FunctionSampler(func=USSampler,
                                    kw_args={'random_state': seed})
            us_X, us_lab = fs_us.fit_resample(us_X, us_lab)

            tmp_X = np.concatenate((np.delete(us_X, np.where(us_lab == central_cluster)[0], axis=0),
                                    np.delete(os_X, np.where(os_lab == central_cluster)[0], axis=0),
                                    norm_X[central_idx, :]))
            tmp_lab = np.concatenate((np.delete(us_lab, np.where(us_lab == central_cluster)[0], axis=0),
                                      np.delete(os_lab, np.where(os_lab == central_cluster)[0], axis=0),
                                      label_predict[central_idx]))

        counter = Counter(tmp_lab)
        print('>> balanced data size: ', tmp_X.shape)
        print('>> balanced clusters distribution: ', [freq for c, freq in counter.items()])
    else:
        tmp_X = norm_X.copy()
        tmp_lab = label_predict.copy()

    print(">> Running Random Forest Model......")
    rf = RandomForestClassifier(n_estimators=1000, n_jobs=-1, random_state=seed)
    rf.fit(tmp_X, tmp_lab)

    # 4. Feature selection
    print(">> selecting genes with high importance ......")
    gene_imp = rf.feature_importances_

    if n_features is None:
        std = np.std(gene_imp)
        thd_mean = np.mean(gene_imp)
        thd_mean3sd = thd_mean + 3 * std
        N = np.sum(gene_imp >= thd_mean3sd)
        if RR:
            sub_X_selected = pd.DataFrame(data=norm_X[:, np.argsort(-gene_imp)[0:N]],
                                          columns=np.argsort(-gene_imp)[0:N])
            sg = filter_bycorr_with_orderly_genes(df=sub_X_selected, threshold=corr_threshold)
        else:
            sg = np.where(gene_imp > thd_mean3sd)[0]
    else:
        sg = np.argsort(-gene_imp)[0:n_features]


    end = time.time()
    runningtime = end - start
    print(">> time used:", runningtime)
    print('>> number of selected genes: ', len(sg))
    X_selected = norm_X[:, sg]

    # 5. Output results
    print(">> saving results ......")
    if dataName is None:
        if save_full:
            if balanced:
                np.savetxt(save_path + 'ClustFS_balance_data.txt', tmp_X, fmt='%f')
                np.savetxt(save_path + 'ClustFS_balance_data_label.txt', tmp_lab, fmt='%d')
        np.savetxt(save_path + 'ClustFS_filtered_res.txt', X_selected, fmt='%f')
        if predict:
            np.savetxt(save_path + 'ClustFS_pred_label.txt', label_predict, fmt='%d')
        np.savetxt(save_path + 'ClustFS_gene_imp_res.txt', gene_imp, fmt='%f')
        np.savetxt(save_path + 'ClustFS_genenames.txt', gene_names, fmt='%s')
        np.savetxt(save_path + 'ClustFS_gs_res.txt', gene_imp, fmt='%f')
    else:
        if save_full:
            if balanced:
                np.savetxt(save_path + dataName + '_ClustFS_balance_data.txt', tmp_X, fmt='%f')
                np.savetxt(save_path + dataName + '_ClustFS_balance_data_label.txt', tmp_lab, fmt='%d')
        np.savetxt(save_path + dataName + '_ClustFS_filtered_res.txt', X_selected, fmt='%f')
        if predict:
            np.savetxt(save_path + dataName + '_ClustFS_pred_label.txt', label_predict, fmt='%d')
        np.savetxt(save_path + dataName + '_ClustFS_gene_imp_res.txt', gene_imp, fmt='%f')
        np.savetxt(save_path + dataName + '_ClustFS_genenames.txt', gene_names, fmt='%s')
        np.savetxt(save_path + dataName + '_ClustFS_gs_res.txt', gene_names[sg], fmt='%s')

    return X_selected

