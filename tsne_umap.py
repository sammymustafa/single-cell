from setup import *

## Run PCA and plot
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca(adata)

## Number of PCs selection
sc.pl.pca_variance_ratio(adata, log=True)

## Compute neighborhood graph
sc.pp.neighbors(adata, n_pcs=4)

## Run clustering and plot
sc.tl.leiden(adata, key_added = 'leiden')
sc.tl.umap(adata)
sc.pl.umap(adata, color='leiden')

## Run tSNE and plot, colored by cluster id which you get in the clustering step, colored by "major.celltype" in the metadata
sc.tl.tsne(adata, n_pcs=4)
sc.pl.tsne(adata, color='leiden')
sc.pl.tsne(adata, color='major.celltype')

## Run UMAP and plot, colored by cluster id which you get in the clustering step, colored by "major.celltype" in the metadata
sc.tl.umap(adata)
sc.pl.umap(adata, color='leiden')

# Plot UMAP, colored by "major.celltype"
sc.pl.umap(adata, color='major.celltype')