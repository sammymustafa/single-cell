from setup import *

# Find marker genes
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')

# This is a dictionary to hold the cell type annotations for each cluster
cell_type_annotations = {}

# Reference marker genes for brain cell types
reference_markers = {
    'Neurons': ['SYT1', 'SNAP25', 'GRIN1'],
    'Excitatory neurons': ['NRGN', 'SLC17A7', 'CAMK2A'],
    'Inhibitory neurons': ['GAD1', 'GAD2'],
    'Astrocytes': ['AQP4', 'GFAP'],
    'Oligodendrocytes': ['MBP', 'MOBP', 'PLP1'],
    'Microglia': ['CSF1R', 'CD74', 'C3'],
    'Oligodendrocyte progenitor cells': ['VCAN', 'PDGFRA', 'CSPG4'],
    'Vascular cells': ['FLT1', 'CLDN5', 'AMBP']
}

# Determine cell type annotation based on marker genes
for cluster_id in adata.obs['leiden'].cat.categories:
    markers = adata.uns['rank_genes_groups']['names'][cluster_id]
    for cell_type, marker_genes in reference_markers.items():
        if any(marker in markers[:10] for marker in marker_genes):  # Check top 10 markers
            if cluster_id not in cell_type_annotations:
                cell_type_annotations[cluster_id] = []
            cell_type_annotations[cluster_id].append(cell_type)

# Add the cell type annotations to the adata object
for cluster_id, cell_types in cell_type_annotations.items():
    # This step assumes a cluster can only have one cell type; adjust as needed.
    cell_type = cell_types[0] if len(cell_types) == 1 else 'Mixed'
    adata.obs.loc[adata.obs['leiden'] == cluster_id, 'annotated_cell_type'] = cell_type

# plot UMAP colored with your cell type annotation and compare it with the "major.celltype" in the metadata
sc.pl.umap(adata, color='annotated_cell_type', title='UMAP colored by Annotated Cell Type')
sc.pl.umap(adata, color='major.celltype', title='UMAP colored by Major Cell Type from Metadata')