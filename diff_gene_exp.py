from setup import *

# Your code and plot
cell_types = adata.obs['major.celltype'].unique()

# Dictionary to store DEG results
degs_wilcox = {}
degs_ttest = {}

for cell_type in cell_types:
    # Subset the data for the current cell type
    adata_subset = adata[adata.obs['major.celltype'] == cell_type]

    # Wilcoxon test
    sc.tl.rank_genes_groups(adata_subset, groupby='condition', method='wilcoxon')
    degs_wilcox[cell_type] = adata_subset.uns['rank_genes_groups']

    # t-test
    sc.tl.rank_genes_groups(adata_subset, groupby='condition', method='t-test')
    degs_ttest[cell_type] = adata_subset.uns['rank_genes_groups']

# Function to plot top DEGs
def plot_top_degs(degs, cell_type, condition, method, n_genes=10):
    sc.pl.rank_genes_groups_dotplot(adata_subset, groupby='condition',
                                    groups=[condition], n_genes=n_genes,
                                    title=f"Top {n_genes} DEGs for {cell_type} - {condition} ({method})")

# Plot for each cell type and condition
for cell_type in cell_types:
    plot_top_degs(degs_wilcox, cell_type, 'AD', 'Wilcoxon')
    plot_top_degs(degs_wilcox, cell_type, 'nonAD', 'Wilcoxon')
    plot_top_degs(degs_ttest, cell_type, 'AD', 't-test')
    plot_top_degs(degs_ttest, cell_type, 'nonAD', 't-test')

# Function to extract top DEG names
def get_deg_names(degs, cell_type, condition, n_genes=10):
    return set(degs[cell_type]['names'][condition][:n_genes])

# Function to plot Venn diagram
import matplotlib.pyplot as plt
def plot_venn(cell_type, condition, degs_wilcox, degs_ttest):
    plt.figure(figsize=(8, 8))
    venn2([get_deg_names(degs_wilcox, cell_type, condition),
           get_deg_names(degs_ttest, cell_type, condition)],
          set_labels=('Wilcoxon', 't-test'))
    plt.title(f"Venn Diagram of DEGs for {cell_type} - {condition}")
    plt.show()

# Plot Venn diagram for each cell type and condition
!pip install matplotlib-venn
from matplotlib_venn import venn2
for cell_type in cell_types:
    plot_venn(cell_type, 'AD', degs_wilcox, degs_ttest)
    plot_venn(cell_type, 'nonAD', degs_wilcox, degs_ttest)