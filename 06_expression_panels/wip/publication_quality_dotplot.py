import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches

# Set publication-quality style
plt.style.use('default')
sns.set_style("whitegrid", {"axes.spines.right": False, "axes.spines.top": False})
plt.rcParams.update({
    'font.size': 12,
    'font.family': 'Arial',
    'axes.linewidth': 1.2,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'xtick.major.width': 1.2,
    'ytick.major.width': 1.2,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.grid': False
})

print("Loading single cell data...")
adata_sc = sc.read_h5ad("/Users/yashkulkarni/cellxgene_data/processed_deconvolution_sc_spleen_240713_fixed_full_data.h5ad")

# Filter for 3wk (infected) timepoint only
print("Filtering for 3wk (infected) timepoint...")
infected_mask = adata_sc.obs['timepoint'] == '3wk'
adata_infected = adata_sc[infected_mask].copy()

print(f"Original data: {adata_sc.shape}")
print(f"3wk infected data: {adata_infected.shape}")

# VERIFY DATA IS PROPERLY NORMALIZED
print("\n=== DATA NORMALIZATION CHECK ===")
print(f"Data type: {adata_infected.X.dtype}")
print(f"Data range: {adata_infected.X.min():.3f} to {adata_infected.X.max():.3f}")
print(f"Mean expression: {adata_infected.X.mean():.3f}")

# Check if data looks like log-normalized (should be mostly 0-10 range)
sample_values = adata_infected.X[0:100, 0:100]
if hasattr(sample_values, 'toarray'):
    sample_values = sample_values.toarray()

print(f"Sample values (first 5x5):")
print(sample_values[:5, :5])

# If data doesn't look normalized, normalize it
if adata_infected.X.max() > 15 or adata_infected.X.min() < 0:
    print("⚠️  Data appears to be raw counts. Applying normalization...")
    # Save raw counts
    adata_infected.layers['raw_counts'] = adata_infected.X.copy()
    # Normalize
    sc.pp.normalize_total(adata_infected, target_sum=1e4)
    sc.pp.log1p(adata_infected)
    print("✅ Normalization applied")
else:
    print("✅ Data appears already normalized")

# Genes organized by functional pathways for publication
gene_groups = {
    'Chemokines': ['Cxcl9', 'Cxcl10', 'Ccl5', 'Ccl6'],
    'Chemokine\nReceptors': ['Cxcr3', 'Ccr5', 'Ccr6'], 
    'IFN-γ\nSignaling': ['Ifng', 'Ifngr1', 'Ifngr2'],
    'IL-18\nSignaling': ['Il18', 'Il18bp', 'Il18r1', 'Il18rap'],
    'IL-27\nSignaling': ['Il27', 'Il6st', 'Il27ra']
}

# Find available genes
available_genes = []
gene_to_group = {}
group_positions = {}

for group, genes in gene_groups.items():
    group_start = len(available_genes)
    for gene in genes:
        found_gene = None
        if gene in adata_infected.var_names:
            found_gene = gene
        else:
            # Try different capitalizations
            for alt in [gene.upper(), gene.lower(), gene.capitalize()]:
                if alt in adata_infected.var_names:
                    found_gene = alt
                    break
        
        if found_gene:
            available_genes.append(found_gene)
            gene_to_group[found_gene] = group
    
    if len(available_genes) > group_start:
        group_positions[group] = (group_start, len(available_genes) - 1)

print(f"\nAvailable genes: {available_genes}")
print(f"Group positions: {group_positions}")

# Create publication-quality figure
fig, ax = plt.subplots(figsize=(16, 8))

# Calculate mean expression and detection rate per cell type
cell_types = sorted(adata_infected.obs['coarse_redo'].unique())
n_cells = len(cell_types)
n_genes = len(available_genes)

# Get gene indices
gene_indices = [list(adata_infected.var_names).index(gene) for gene in available_genes]
X_subset = adata_infected.X[:, gene_indices]
if hasattr(X_subset, 'toarray'):
    X_subset = X_subset.toarray()

# Calculate mean expression and detection rate
mean_expr = np.zeros((n_cells, n_genes))
detection_rate = np.zeros((n_cells, n_genes))

for i, cell_type in enumerate(cell_types):
    mask = adata_infected.obs['coarse_redo'] == cell_type
    cell_data = X_subset[mask]
    
    # Mean expression
    mean_expr[i, :] = cell_data.mean(axis=0)
    
    # Detection rate (percentage of cells expressing)
    detection_rate[i, :] = (cell_data > 0).mean(axis=0) * 100

# Z-score normalize across cell types for each gene (gene-wise normalization)
mean_expr_zscore = np.zeros_like(mean_expr)
for j in range(n_genes):
    gene_values = mean_expr[:, j]
    if gene_values.std() > 0:
        mean_expr_zscore[:, j] = (gene_values - gene_values.mean()) / gene_values.std()
    else:
        mean_expr_zscore[:, j] = 0

print(f"\n=== NORMALIZATION VERIFICATION ===")
print(f"Z-score range: {mean_expr_zscore.min():.3f} to {mean_expr_zscore.max():.3f}")
print(f"Detection rate range: {detection_rate.min():.1f}% to {detection_rate.max():.1f}%")

# Create the dotplot manually for full control
x_pos = np.arange(n_genes)
y_pos = np.arange(n_cells)
X, Y = np.meshgrid(x_pos, y_pos)

# Flatten for scatter plot
x_flat = X.flatten()
y_flat = Y.flatten()
colors = mean_expr_zscore.flatten()
sizes = detection_rate.flatten()

# Scale sizes appropriately (20 to 400 pixels)
sizes_scaled = 20 + (sizes / 100) * 380

# Create scatter plot
scatter = ax.scatter(x_flat, y_flat, c=colors, s=sizes_scaled, 
                    cmap='RdBu_r', vmin=-2.5, vmax=2.5, 
                    alpha=0.8, edgecolors='black', linewidth=0.3)

# Customize axes
ax.set_xlim(-0.5, n_genes - 0.5)
ax.set_ylim(-0.5, n_cells - 0.5)

# Set ticks and labels
ax.set_xticks(range(n_genes))
ax.set_xticklabels([f'${gene}$' for gene in available_genes], 
                   rotation=45, ha='right', fontsize=11)
ax.set_yticks(range(n_cells))
ax.set_yticklabels(cell_types, fontsize=12)

# Add grid for readability
ax.set_axisbelow(True)
ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)

# Labels
ax.set_xlabel('Genes', fontsize=14, fontweight='bold', labelpad=10)
ax.set_ylabel('Cell Types', fontsize=14, fontweight='bold', labelpad=10)

# Add functional group brackets and labels
bracket_y = n_cells + 0.3
for group, (start, end) in group_positions.items():
    # Draw bracket
    ax.plot([start - 0.3, start - 0.3, end + 0.3, end + 0.3], 
            [bracket_y, bracket_y + 0.15, bracket_y + 0.15, bracket_y], 
            'k-', lw=1.5)
    
    # Add group label
    mid_pos = (start + end) / 2
    ax.text(mid_pos, bracket_y + 0.25, group, ha='center', va='bottom',
            fontsize=11, fontweight='bold')

# Add colorbars
# Expression colorbar
cbar1 = plt.colorbar(scatter, ax=ax, shrink=0.6, aspect=20, pad=0.02)
cbar1.set_label('Normalized Expression\n(Z-score)', fontsize=12, fontweight='bold')
cbar1.ax.tick_params(labelsize=10)

# Size legend for detection rate
legend_sizes = [20, 40, 60, 80, 100]
legend_handles = []
for size in legend_sizes:
    size_scaled = 20 + (size / 100) * 380
    handle = plt.scatter([], [], s=size_scaled, c='gray', alpha=0.7, 
                        edgecolors='black', linewidth=0.3)
    legend_handles.append(handle)

size_legend = ax.legend(legend_handles, [f'{s}%' for s in legend_sizes],
                       title='Fraction of Cells\nExpressing Gene',
                       loc='center left', bbox_to_anchor=(1.15, 0.3),
                       frameon=True, fancybox=True, shadow=True,
                       title_fontsize=12, fontsize=10)
size_legend.get_title().set_fontweight('bold')

# Final styling
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(1.2)
ax.spines['bottom'].set_linewidth(1.2)

# Adjust layout
plt.subplots_adjust(top=0.85, bottom=0.15, left=0.12, right=0.75)

# Add title
fig.suptitle('Gene Expression Patterns in Infected Spleen Cell Types\n(Gene-wise Normalized)', 
             fontsize=16, fontweight='bold', y=0.95)

# Save high-resolution figures
plt.savefig('publication_dotplot_gene_expression_FINAL.pdf', dpi=300, bbox_inches='tight')
plt.savefig('publication_dotplot_gene_expression_FINAL.png', dpi=300, bbox_inches='tight')

plt.show()

print("\n✅ PUBLICATION-QUALITY DOTPLOT COMPLETED!")
print("✅ Data normalization verified and applied")
print("✅ Gene-wise Z-score normalization applied")
print("✅ Functional groupings added")
print("✅ High-resolution figures saved:")
print("   - publication_dotplot_gene_expression_FINAL.pdf")
print("   - publication_dotplot_gene_expression_FINAL.png")
