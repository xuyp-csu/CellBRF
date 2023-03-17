
import numpy as np
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_samples
from sklearn.feature_selection import f_classif, SelectKBest, mutual_info_classif
import pandas as pd
import os
import scanpy as sc
import h5py
import time

# Compute Clustering Score
def ClusteringScore(X, labels):

    s_score = silhouette_samples(X, labels)
    s_average = np.mean(s_score)

    return s_average


def SelectTopNScores(clustering_score, n):

    idx = np.argsort(clustering_score, kind='mergesort')
    idx = np.flip(idx)
    #print("s_score array size: ", idx.size)
    # print(clustering_score)
    if (idx.size < n):
        return idx
    else:
        return idx[0:n]

def normalize(adata):
    sc.pp.filter_genes(adata, min_cells=int((10 / 100) * adata.shape[0]))
    sc.pp.filter_genes(adata, max_cells=int((90 / 100) * adata.shape[0]))
    sc.pp.log1p(adata)
    return adata

# data load
dataList = ["xxx"]
print(dataList)

for dataName in dataList:
    print(dataName)
    dir = './'
    # data load
    counts_df = pd.read_csv(os.path.join(dir, dataName + '.csv'),encoding='gbk')
    lab_df = pd.read_csv(os.path.join(dir, dataName + '_label.csv'))
    var_names = counts_df.iloc[:, 0].values.tolist()
    obs_names = counts_df.columns[1:].values.tolist()
    lab = lab_df.iloc[:, 1].values.tolist()
    X = counts_df.iloc[:, 1:].values

    # data_mat = h5py.File(os.path.join(dir, dataName + '.h5'))
    # X = np.array(data_mat['X'])
    # y = np.array(data_mat['Y'])
    # data_mat.close()
    # obs_names = ['cell'+str(i) for i in range(1, X.shape[0]+1)]
    # var_names = ['gene'+str(i) for i in range(1, X.shape[1]+1)]

    # adata = sc.AnnData(X.T)
    adata = sc.AnnData(X)
    print(adata)
    adata.obs['Group'] = y
    adata.obs_names = obs_names
    adata.var_names = var_names
    n_clusters = len(np.unique(adata.obs['Group']))
    print("dataset: %s" % dataName)
    print("cells: %d; genes: %d; n clusters: %d" % (len(obs_names), len(var_names), n_clusters))

    adata = normalize(adata)
    print(adata)
    X = adata.X.copy()

    sc.pp.scale(adata)
    X_nrm = adata.X.copy()

    q = np.arange(1, int((5 / 100) * adata.shape[1]) + 1)
    # q = np.arange(100, int((5 / 100) * adata.X.shape[1]) + 1, 10)

    # FEATURE SELECTION STEP
    # Using Hierarchical Clustering, compute and store temporary clusters in sc object
    print("Computing Temporary Clusters . . .")
    model_hc = AgglomerativeClustering(n_clusters=n_clusters, linkage='ward')
    temp_labels = model_hc.fit_predict(X_nrm) + 1

    start = time.time()
    # Perform ANOVA analysis to find significant genes
    print("Performing Feature Selection . . .")
    feat_sel = SelectKBest(f_classif, k="all")
    feat_sel.fit(X_nrm, temp_labels)
    feature_scores = feat_sel.scores_
    idx = np.argsort(feature_scores, kind='mergesort')
    idx = idx[::-1]  # Sort descending

    # CLUSTERING STEP
    clust_score = np.zeros(q.shape[0])  # for storing clustering score
    pred_labels = np.zeros([X_nrm.shape[0], q.shape[0]])  # for storing cluster labels

    # compute clusters and iterate
    print("Computing " + str(n_clusters) + " clusters using best features . . .")
    i = 0
    for j in q:
        X_red = X_nrm[:, idx[0:j]]
        pred_labels[:, i] = model_hc.fit_predict(X_red) + 1
        clust_score[i] = ClusteringScore(X, pred_labels[:, i])
        i = i + 1

    print(clust_score)
    print(q)
    # FINAL CLUSTERING
    # Choose clustering with the max. clustering score
    mask = SelectTopNScores(clust_score, 1)
    final_labels = pred_labels[:, mask]
    final_labels = final_labels.astype(int)

    print("Optimal number of features = ", mask+1)

    end = time.time()
    runningtime = end - start
    print(">> time used:", runningtime)
    np.savetxt('./Re_sg_res/Feats_' + dataName + '.txt', np.array(adata.var_names)[idx[range(0, int(q[mask]))]],
               fmt="%s", delimiter=" ")



