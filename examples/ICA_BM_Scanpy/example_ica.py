import scanpy.api as sc

print("Reading ICA (bone marrow) dataset")
adata = sc.read_10x_h5("ica_bone_marrow_h5.h5", genome='GRCh38')
adata.var_names_make_unique()

print("Filtering cells")
sc.pp.filter_cells(adata, min_genes=500)

print("Filtering genes")
sc.pp.filter_genes(adata, min_cells=3)

print("Normalizing data")
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)

print("Finding variable genes")
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, inplace=True)
adata = adata[:, adata.var['highly_variable']]

print("Scale expression matrix of variable genes")
sc.pp.scale(adata, max_value=10)

print("Running PCA")
sc.tl.pca(adata, svd_solver='arpack')

print("Computing neighborhood graph")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

print("Finding Clusters")
sc.tl.louvain(adata)

print("Computing UMAP embedding")
sc.tl.umap(adata)

print("Computing tSNE embedding")
sc.tl.tsne(adata)
#adata.X = None  # don't save matrix
adata.write('ica_bone_marrow_sc.h5ad', compression='gzip')
