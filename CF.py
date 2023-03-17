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
# from RF import RandomForestClassifier
from collections import Counter
from sklearn.neighbors import kneighbors_graph, NearestNeighbors
from sklearn.decomposition import PCA
import math
import time
from imblearn.over_sampling import SMOTE, BorderlineSMOTE
from imblearn.under_sampling import NearMiss, InstanceHardnessThreshold
from imblearn import FunctionSampler
from kneed import KneeLocator
from sklearn.metrics import silhouette_samples, calinski_harabasz_score
from scipy.spatial.distance import euclidean, cdist, pdist, squareform
from sklearn.metrics.cluster import normalized_mutual_info_score, adjusted_rand_score
from sklearn.semi_supervised import LabelPropagation, LabelSpreading
from gap_statistic import OptimalK

def dopca(X, dim):
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


def predict_label(data, k_clusters, pca, seed, k):
    adj = get_adj(count=data, pca=pca, k=k)
    label_predict = SpectralClustering(n_clusters=k_clusters,
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

    # iht_us = InstanceHardnessThreshold(sampling_strategy=ratio, random_state=random_state)
    iht_us = InstanceHardnessThreshold(estimator=RandomForestClassifier(random_state=random_state, bootstrap=False),
                                       sampling_strategy=ratio, random_state=random_state, cv=5, n_jobs=None)
    _, tmp_reshape_y = iht_us.fit_resample(x, tmp_y)
    id = iht_us.sample_indices_

    reshape_x = x[id, :].copy()

    decode = dict(zip(range(len(np.unique(y))), np.unique(y)))
    reshape_y = []
    for i in tmp_reshape_y:
        reshape_y.append(decode[i])

    return reshape_x, reshape_y


def balance(seq):
    n = len(seq)
    classes = [(clas, float(count)) for clas, count in Counter(seq).items()]
    k = len(classes)
    H = -sum([(count / n) * np.log((count / n)) for clas, count in classes])  # shannon entropy
    return 1-(H / np.log(k))

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
            k=15,
            npcs=50,
            seed=2022,
            sub_factor=0.8,
            RR=True,
            corr_threshold=0.8,
            save_full=True,
            n_features=None,
            true_lab=None,
            balanced=True,
            under_sampling=True,
            over_sampling=True,
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

    k : integer
        Number of nearest neighbors considered in the KNN graph. (default: 15)

    npcs : integer
        Number of PC considered in the PCA step. (default: 50)

    seed : integer
        Random seed.

    corr_threshold : float
        Linear correlation threshold to use when removing redundant genes. (default: 0.8)

    RR : boolean
        Whether redundant features need to be removed. (default: True)

    sub_factor : float
        Proportion of cells retained in major clusters. (default: 0.8)

    save_full : boolean
        Whether the full result needs to be saved. (default: True)

    save_path : string
        Path to save results.

    Returns
    ----------
    X_selected : `numpy.ndarray`
        Gene expression matrix after feature selection.
    """

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

    pca = min(npcs, norm_X.shape[0])
    # 1. Label prediction
    if label_predict is None:
        print(">> Spectral clustering to generate predicted labels......")
        if n_clusters is None:
            optimalK = OptimalK(n_jobs=8, parallel_backend='joblib')
            n_clusters = optimalK(norm_X, cluster_array=np.arange(2, 15))
            label_predict = predict_label(data=norm_X, k_clusters=n_clusters, seed=seed, pca=pca, k=k)
        else:
            label_predict = predict_label(data=norm_X, k_clusters=n_clusters, seed=seed, pca=pca, k=k)
        predict = True
    else:
        predict = False

    label_predict = np.array([int(i) for i in label_predict])

    if true_lab is not None:
        true_lab = np.array([int(i) for i in true_lab])
        nmi = normalized_mutual_info_score(true_lab, label_predict)
        ari = adjusted_rand_score(true_lab, label_predict)
        print('>> Pre-label NMI:', round(nmi, 4), 'ARI:', round(ari, 4))
        print('>> Raw Balance entropy: ', balance(true_lab))
        print('>> Pred Balance entropy: ', balance(label_predict))

    print('>> Data size: ', norm_X.shape)
    start = time.time()
    # 2. data balance
    sub_id = []
    if balanced:
        if norm_X.shape[0] < 150 and balance(label_predict) > 0.95:
            tmp_X = norm_X.copy()
            tmp_lab = label_predict.copy()
            balanced = False
        else:
            print(">> data balancing......")
            if under_sampling:
                # Improved label accuracy
                print(">> Improved label accuracy......")
                # sub_factor = 0.8
                for c in np.unique(label_predict):
                    tmp_idx = np.where(label_predict == c)[0]
                    if len(tmp_idx) >= (norm_X.shape[0] / n_clusters):
                        num = int(len(tmp_idx) * sub_factor)
                        sub = norm_X[tmp_idx, :]
                        center = np.mean(sub, axis=0)
                        pdist = (sub - center) ** 2
                        sub_id.extend(tmp_idx[np.argsort(np.mean(pdist, axis=1))[0:num]])
                    else:
                        sub_id.extend(tmp_idx)
                sub_norm_X = norm_X[sub_id, :].copy()
                sub_lab = label_predict[sub_id].copy()
                if true_lab is not None:
                    nmi = normalized_mutual_info_score(true_lab[sub_id], sub_lab)
                    ari = adjusted_rand_score(true_lab[sub_id], sub_lab)
                    print('>> Sub size:', len(sub_id), 'Ratio:', round(len(sub_id) / len(label_predict), 2))
                    print('>> Sub-label NMI:', round(nmi, 4), 'ARI:', round(ari, 4))
                    Shannon_entropy = balance(true_lab[sub_id])
                    print('>> Raw Balance entropy after under sampling: ', Shannon_entropy)

            else:
                sub_norm_X = norm_X.copy()
                sub_lab = label_predict.copy()

            n_cells = sub_norm_X.shape[0]
            counter = Counter(sub_lab)

            if under_sampling:
                cdtb = [freq for c, freq in counter.items()]
                print('>> Sub data size: ', sub_norm_X.shape)
                print('>> clusters distribution: ', cdtb)
                print('>> balanced: ', balanced)

            if over_sampling:
                # over sampling
                balanced_thd = math.ceil(n_cells / n_clusters)
                min_cluster = 10
                diff = []
                clusters = []
                for i in counter.items():
                    if i[1] < min_cluster:
                        min_cluster = i[1]
                    diff.append(i[1] - balanced_thd)
                    clusters.append(i[0])
                if n_clusters == 2:
                    resize = []
                    for c in np.unique(sub_lab):
                        idx = list(np.where(sub_lab == c)[0])
                        if len(idx) > balanced_thd:
                            resize.append(len(idx))
                        else:
                            resize.append(balanced_thd)
                    ratio = dict(zip(np.unique(sub_lab), resize))
                    if min_cluster < 7:
                        smo_os = SMOTE(random_state=seed, k_neighbors=min_cluster - 1, sampling_strategy=ratio)
                    else:
                        smo_os = SMOTE(random_state=seed, sampling_strategy=ratio)
                else:
                    central_cluster = clusters[np.where(np.abs(diff) == min(np.abs(diff)))[0][0]]
                    central_idx = list(np.where(sub_lab == central_cluster)[0])
                    # Over sampling
                    resize = []
                    for c in np.unique(sub_lab):
                        idx = list(np.where(sub_lab == c)[0])
                        if len(idx) > len(central_idx):
                            resize.append(len(idx))
                        else:
                            resize.append(len(central_idx))
                    ratio = dict(zip(np.unique(sub_lab), resize))
                    if min_cluster < 7:
                        smo_os = SMOTE(random_state=seed, k_neighbors=min_cluster-1, sampling_strategy=ratio)
                    else:
                        smo_os = SMOTE(random_state=seed, sampling_strategy=ratio)

                tmp_X, tmp_lab = smo_os.fit_resample(sub_norm_X, sub_lab)
                counter = Counter(tmp_lab)
                print('>> balanced data size: ', tmp_X.shape)
                print('>> balanced clusters distribution: ', [freq for c, freq in counter.items()])
            else:
                tmp_X = sub_norm_X.copy()
                tmp_lab = sub_lab.copy()

    else:
        tmp_X = norm_X.copy()
        tmp_lab = label_predict.copy()

    # 4. Feature Importance Assessment
    print(">> Running Random Forest Model......")
    rf = RandomForestClassifier(n_estimators=1000, n_jobs=-1, random_state=seed, bootstrap=False)
    rf.fit(tmp_X, tmp_lab)

    # 5. Feature selection
    print(">> selecting genes with high importance ......")
    gene_imp = rf.feature_importances_.copy()

    if n_features is None:
        std = np.std(gene_imp)
        thd_mean = np.mean(gene_imp)
        thd_mean3sd = thd_mean + 3 * std
        N = np.sum(gene_imp >= thd_mean3sd)
    else:
        N = n_features

    if RR:
        if balanced:
            sub_X_selected = pd.DataFrame(data=sub_norm_X[:, np.argsort(-gene_imp)[0:N]],
                                          columns=np.argsort(-gene_imp)[0:N])
        else:
            sub_X_selected = pd.DataFrame(data=norm_X[:, np.argsort(-gene_imp)[0:N]],
                                          columns=np.argsort(-gene_imp)[0:N])
        sg = filter_bycorr_with_orderly_genes(df=sub_X_selected, threshold=corr_threshold)
    else:
        sg = np.argsort(-gene_imp)[0:N]

    end = time.time()
    runningtime = end - start
    print(">> time used:", runningtime)
    print('>> number of selected genes: ', len(sg))
    X_selected = norm_X[:, sg]

    # Output results
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
        np.savetxt(save_path + 'ClustFS_gs_res.txt', gene_names[sg], fmt='%f')
        if under_sampling:
            np.savetxt(save_path + 'ClustFS_subid.txt', sub_id, fmt='%d')
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
        if under_sampling:
            np.savetxt(save_path + dataName + '_ClustFS_subid.txt', sub_id, fmt='%d')

    return X_selected

