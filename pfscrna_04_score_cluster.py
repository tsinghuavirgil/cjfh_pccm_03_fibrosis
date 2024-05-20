from sklearn.cluster import KMeans
from numpy import unique
from numpy import where
##########
score_mtx = pd.read_csv('score_adata.csv',header =0, index_col=0)
col_mtx = score_mtx['ecm_02_collagen_score']
col_mtx = np.array(col_mtx)
cluster_col = col_mtx.reshape(-1,1)
kmeans_res = KMeans(n_clusters=3).fit(cluster_col)
cluster_label = kmeans_res.labels_
kmeans_pre = kmeans_res.fit_predict(cluster_col)
new_column = pd.DataFrame(kmeans_pre)
new_column.index=score_mtx.index
new_column.columns=['kmeans_res']
kmeans_mtx = score_mtx.merge(new_column,left_index=True,right_index=True)
kmeans_mtx.to_csv('kmean_adata.csv')