import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches
from scipy import sparse

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
adata_sc = sc.read_h5ad("/Users/yashkulkarni/cellxgene_data/sc_spleen_072225.h5ad")

print(f"Original data: {adata_sc.shape}")
print(f"Timepoints: {adata_sc.obs['timepoint'].unique()}")

# NORMALIZE TOGETHER for fair comparison
print("\n=== NORMALIZING DATA TOGETHER ===")
if adata_sc.X.max() > 15 or adata_sc.X.min() < 0:
    print("⚠️  Data appears to be raw counts. Applying normalization...")
    adata_sc.layers['raw_counts'] = adata_sc.X.copy()
    sc.pp.normalize_total(adata_sc, target_sum=1e4)
    sc.pp.log1p(adata_sc)
    print("✅ Normalization applied to full dataset")
else:
    print("✅ Data appears already normalized")

# Separate infected and uninfected AFTER normalization
infected_mask = adata_sc.obs['timepoint'] == '3wk'
uninfected_mask = adata_sc.obs['timepoint'] == '0wk'

adata_infected = adata_sc[infected_mask].copy()
adata_uninfected = adata_sc[uninfected_mask].copy()

print(f"3wk infected data: {adata_infected.shape}")
print(f"0wk uninfected data: {adata_uninfected.shape}")

# Explicitly exclude specific cell types
exclude_celltypes = [
    'Erythrocyte',
    'Erythrocyte_precursor', 
    'Thrombocyte',
    'Hematopoietic_stem_cell',
    'Hematopoeitic_stem_cell',
    'CD4-Tcell_naive',
    'Bcell_plasma',
    'Myeloid_pDC'
]

# Get all cell types present in both conditions (excluding unwanted ones)
all_celltypes_infected = set(adata_infected.obs['paper_names'].unique())
all_celltypes_uninfected = set(adata_uninfected.obs['paper_names'].unique())
all_celltypes = sorted(all_celltypes_infected & all_celltypes_uninfected)

# Remove excluded cell types
all_celltypes = [ct for ct in all_celltypes if ct not in exclude_celltypes]

print(f"\nCell types present in both conditions: {len(all_celltypes)}")

# Use the SAME gene list from the previous plot
module_labels = [
    ('Chemokines', ['Cxcl9', 'Cxcl10', 'Ccl5', 'Cxcl16']),
    ('Chemokine Receptors', ['Cxcr3', 'Cxcr5', 'Cxcr6', 'Ccr5']),
    ('IFN-γ Signaling', ['Ifng','Ifngr1', 'Ifngr2', 'Ifngas1']),
    ('IL-18 Signaling', ['Il18', 'Il18bp', 'Il18r1', 'Il18rap']),
    ('IL-27 Signaling', ['Il27', 'Ebi3', 'Il6st', 'Il27ra']),
    ('Activation', ['Cd80', 'Cd86', 'H2-K1'])
]

# Get all genes in order
genes_passing_filter = []
for family_name, gene_list in module_labels:
    for gene in gene_list:
        # Find the gene (case insensitive)
        found_gene = None
        for var_gene in adata_infected.var_names:
            if gene.lower() == var_gene.lower():
                found_gene = var_gene
                break
        
        if found_gene is not None and found_gene not in genes_passing_filter:
            genes_passing_filter.append(found_gene)

print(f"\nUsing same genes as previous plot: {len(genes_passing_filter)} genes")
print(f"Genes: {genes_passing_filter}")

# Get gene indices
gene_indices = [list(adata_infected.var_names).index(gene) for gene in genes_passing_filter]

# Infected data
X_infected = adata_infected.X[:, gene_indices]
if hasattr(X_infected, 'toarray'):
    X_infected = X_infected.toarray()

# Uninfected data
X_uninfected = adata_uninfected.X[:, gene_indices]
if hasattr(X_uninfected, 'toarray'):
    X_uninfected = X_uninfected.toarray()

n_cells = len(all_celltypes)
n_genes = len(genes_passing_filter)

# Calculate mean expression and detection rate for each cell type
mean_expr_infected = np.zeros((n_cells, n_genes))
detection_rate_infected = np.zeros((n_cells, n_genes))
mean_expr_uninfected = np.zeros((n_cells, n_genes))
detection_rate_uninfected = np.zeros((n_cells, n_genes))

for i, cell_type in enumerate(all_celltypes):
    # Infected
    mask_inf = adata_infected.obs['paper_names'] == cell_type
    if mask_inf.sum() > 0:
        cell_data_inf = X_infected[mask_inf]
        mean_expr_infected[i, :] = cell_data_inf.mean(axis=0)
        detection_rate_infected[i, :] = (cell_data_inf > 0).mean(axis=0) * 100
    
    # Uninfected
    mask_uninf = adata_uninfected.obs['paper_names'] == cell_type
    if mask_uninf.sum() > 0:
        cell_data_uninf = X_uninfected[mask_uninf]
        mean_expr_uninfected[i, :] = cell_data_uninf.mean(axis=0)
        detection_rate_uninfected[i, :] = (cell_data_uninf > 0).mean(axis=0) * 100

# Z-score normalize TOGETHER for fair comparison
mean_expr_infected_zscore = np.zeros_like(mean_expr_infected)
mean_expr_uninfected_zscore = np.zeros_like(mean_expr_uninfected)

for j in range(n_genes):
    # Combine infected and uninfected for normalization
    combined_values = np.concatenate([mean_expr_infected[:, j], mean_expr_uninfected[:, j]])
    
    if combined_values.std() > 0:
        mean_combined = combined_values.mean()
        std_combined = combined_values.std()
        
        mean_expr_infected_zscore[:, j] = (mean_expr_infected[:, j] - mean_combined) / std_combined
        mean_expr_uninfected_zscore[:, j] = (mean_expr_uninfected[:, j] - mean_combined) / std_combined
    else:
        mean_expr_infected_zscore[:, j] = 0
        mean_expr_uninfected_zscore[:, j] = 0

print(f"\n=== NORMALIZATION VERIFICATION ===")
print(f"Infected Z-score range: {mean_expr_infected_zscore.min():.3f} to {mean_expr_infected_zscore.max():.3f}")
print(f"Uninfected Z-score range: {mean_expr_uninfected_zscore.min():.3f} to {mean_expr_uninfected_zscore.max():.3f}")

# Create side-by-side plots
print("\n=== CREATING SIDE-BY-SIDE PLOTS ===")

fig_width = max(24, n_genes * 0.8)
fig_height = max(15, n_cells * 0.4)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(fig_width, fig_height), 
                                sharey=True, sharex=False)

# Common settings
x_pos = np.arange(n_genes)
y_pos = np.arange(n_cells)
X, Y = np.meshgrid(x_pos, y_pos)

# LEFT PANEL - INFECTED
x_flat = X.flatten()
y_flat = Y.flatten()
colors_inf = mean_expr_infected_zscore.flatten()
sizes_inf = detection_rate_infected.flatten()
sizes_scaled_inf = 20 + (sizes_inf / 100) * 380

scatter1 = ax1.scatter(x_flat, y_flat, c=colors_inf, s=sizes_scaled_inf, 
                      cmap='RdBu_r', vmin=-2.5, vmax=2.5, 
                      alpha=0.8, edgecolors='black', linewidth=0.3)

ax1.set_xlim(-0.5, n_genes - 0.5)
ax1.set_ylim(-0.5, n_cells - 0.5)
ax1.set_xticks(range(n_genes))
ax1.set_xticklabels([f'${gene}$' for gene in genes_passing_filter], 
                    rotation=45, ha='right', fontsize=9)
ax1.set_yticks(range(n_cells))
ax1.set_yticklabels(all_celltypes, fontsize=9)
ax1.set_axisbelow(True)
ax1.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
ax1.set_xlabel('Genes', fontsize=12, fontweight='bold', labelpad=10)
ax1.set_ylabel('Fine Cell Types', fontsize=12, fontweight='bold', labelpad=10)
ax1.set_title('Infected (3wk)', fontsize=14, fontweight='bold', pad=30)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['left'].set_linewidth(1.2)
ax1.spines['bottom'].set_linewidth(1.2)

# RIGHT PANEL - UNINFECTED
colors_uninf = mean_expr_uninfected_zscore.flatten()
sizes_uninf = detection_rate_uninfected.flatten()
sizes_scaled_uninf = 20 + (sizes_uninf / 100) * 380

scatter2 = ax2.scatter(x_flat, y_flat, c=colors_uninf, s=sizes_scaled_uninf, 
                      cmap='RdBu_r', vmin=-2.5, vmax=2.5, 
                      alpha=0.8, edgecolors='black', linewidth=0.3)

ax2.set_xlim(-0.5, n_genes - 0.5)
ax2.set_ylim(-0.5, n_cells - 0.5)
ax2.set_xticks(range(n_genes))
ax2.set_xticklabels([f'${gene}$' for gene in genes_passing_filter], 
                    rotation=45, ha='right', fontsize=9)
ax2.set_axisbelow(True)
ax2.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
ax2.set_xlabel('Genes', fontsize=12, fontweight='bold', labelpad=10)
ax2.set_title('Uninfected (0wk)', fontsize=14, fontweight='bold', pad=30)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['left'].set_linewidth(1.2)
ax2.spines['bottom'].set_linewidth(1.2)

# Add gene group brackets ABOVE both panels
gene_to_group = {}
group_positions = {}

for group, genes in module_labels:
    group_genes_in_order = []
    
    # Find genes in the order they appear in genes_passing_filter
    for gene in genes_passing_filter:
        for target_gene in genes:
            if target_gene.lower() == gene.lower():
                if gene not in group_genes_in_order:
                    group_genes_in_order.append(gene)
                    gene_to_group[gene] = group
                break
    
    if group_genes_in_order:
        # Use actual positions in genes_passing_filter list
        start_idx = genes_passing_filter.index(group_genes_in_order[0])
        end_idx = genes_passing_filter.index(group_genes_in_order[-1])
        group_positions[group] = (start_idx, end_idx)

# Add brackets to BOTH panels
bracket_y = n_cells 
bracket_height = 0.15
text_offset = 0.25

for ax in [ax1, ax2]:
    for group, (start, end) in group_positions.items():
        # Draw horizontal bracket line
        ax.plot([start - 0.2, end + 0.2], [bracket_y, bracket_y], 
                'k-', lw=2, clip_on=False)
        
        # Draw vertical bracket ends pointing down
        ax.plot([start - 0.2, start - 0.2], [bracket_y, bracket_y - bracket_height], 
                'k-', lw=2, clip_on=False)
        ax.plot([end + 0.2, end + 0.2], [bracket_y, bracket_y - bracket_height], 
                'k-', lw=2, clip_on=False)
        
        # Add group label above the bracket
        mid_pos = (start + end) / 2
        ax.text(mid_pos, bracket_y + text_offset, group, 
                ha='center', va='bottom', fontsize=9, fontweight='bold',
                clip_on=False)

# Add shared colorbar to the right of both panels
cbar = fig.colorbar(scatter1, ax=[ax1, ax2], shrink=0.6, aspect=20, pad=0.02)
cbar.set_label('Normalized Expression', fontsize=12, fontweight='bold')
cbar.ax.tick_params(labelsize=10)

# Add shared size legend
legend_sizes = [20, 40, 60, 80, 100]
legend_handles = []
for size in legend_sizes:
    size_scaled = 20 + (size / 100) * 380
    handle = plt.scatter([], [], s=size_scaled, c='gray', alpha=0.7, 
                        edgecolors='black', linewidth=0.3)
    legend_handles.append(handle)

# Position legend to the right of both panels
size_legend = fig.legend(legend_handles, [f'{s}%' for s in legend_sizes],
                        title='Fraction of Cells\nExpressing Gene',
                        loc='center right', bbox_to_anchor=(0.99, 0.5),
                        frameon=True, fancybox=True, shadow=True,
                        title_fontsize=12, fontsize=10,
                        scatterpoints=1,
                        markerscale=1,
                        handletextpad=2.0,
                        borderpad=1.5,
                        labelspacing=1.2)
size_legend.get_title().set_fontweight('bold')

# Adjust layout - give more space on right for colorbar and legend
plt.subplots_adjust(top=0.90, bottom=0.12, left=0.08, right=0.80, wspace=0.08)

# Add overall title
fig.suptitle('Gene Expression: Infected vs Uninfected Comparison', 
             fontsize=16, fontweight='bold', y=0.96)

# Save high-resolution figures
plt.savefig('infected_vs_uninfected_comparison_sidebyside.pdf', dpi=300, bbox_inches='tight')
plt.savefig('infected_vs_uninfected_comparison_sidebyside.png', dpi=300, bbox_inches='tight')

print("\n✅ Plot saved!")
print(f"Dimensions: {n_genes} genes × {n_cells} cell types (2 panels)")
print(f"Showing same genes as previous plot: {len(genes_passing_filter)} genes")
plt.show()
