import h5py
import numpy as np
import scanpy.api as sc

print("Reading ICA (bone marrow) dataset")
adata = sc.read_10x_h5("ica_bone_marrow_h5.h5", genome='GRCh38')
adata.var_names_make_unique()

print("Filtering cells")
sc.pp.filter_cells(adata, min_genes=500)

print("Filtering genes")
sc.pp.filter_genes(adata, min_cells=50)

print("Normalizing data")
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)

print("Finding variable genes")
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, inplace=True)
adata_variable_genes = adata[:, adata.var['highly_variable']]

print("Scale expression matrix of variable genes")
sc.pp.scale(adata_variable_genes, max_value=10)

print("Running PCA")
sc.tl.pca(adata_variable_genes, svd_solver='arpack')

print("Computing neighborhood graph")
sc.pp.neighbors(adata_variable_genes, n_neighbors=10, n_pcs=40)

print("Finding Clusters")
sc.tl.louvain(adata_variable_genes)

print("Computing UMAP embedding")
sc.tl.umap(adata_variable_genes)

print("Computing tSNE embedding")
sc.tl.tsne(adata_variable_genes)

loom_file = h5py.File('ica_bone_marrow.loom', 'w')
# genes along rows in loom file
for name in adata.var:  # genes
    val = np.array(adata.var[name])
    if val.dtype is np.dtype('O'):
        val = val.astype(str)
loom_file['/row_attrs/' + name] = val
loom_file['/row_attrs/Gene'] = adata.var_names.values.astype(np.string_)

# cells along columns in loom file
for name in adata_variable_genes.obsm_keys():
    array = adata_variable_genes.obsm[name]
    for j in range(array.shape[1]):
        loom_file['/col_attrs/' + name + str(j + 1)] = array[:, j]

for name in adata_variable_genes.obs:
    val = np.array(adata_variable_genes.obs[name])
    if val.dtype is np.dtype('O'):
        val = val.astype(np.string_)
    loom_file['/col_attrs/' + name] = val
loom_file['/col_attrs/CellID'] = adata.obs_names.values.astype(np.string_)

# save in chunks in dense format
X = adata.X.transpose()
dset = loom_file.create_dataset('/matrix', shape=X.shape, chunks=(1000, 1000), compression='gzip',
                                compression_opts=9)

start = 0
step = min(X.shape[0], 1000)
stop = step
nchunks = int(np.ceil(max(1, X.shape[0] / step)))
for i in range(nchunks):
    stop = min(X.shape[0], stop)
    dset[start:stop] = X[start:stop].toarray()
    start += step
    stop += step

loom_file.close()
