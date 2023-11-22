from setup import * 

## Read data and create anndata object
# load sparse matrix :
X = io.mmread("data/counts.mtx")
# create anndata object
adata = anndata.AnnData(X=X.transpose().tocsr() )
# load cell metadata:
cell_meta = pd.read_csv("data/meta.csv")
# load gene names:
with open("data/gene_names.csv", 'r') as f:
    gene_names = f.read().splitlines()
# set anndata observations and index obs by barcodes, var by gene names
adata.obs = cell_meta
adata.obs.index = adata.obs['barcode']
adata.var.index = gene_names
print(adata)

## Since the provided data was pre-processed, we will skip the data fitering steps.
## Perform log-transformation, scale data, and identify highly variable genes
sc.pp.log1p(adata)
sc.pp.scale(adata, max_value=10)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)