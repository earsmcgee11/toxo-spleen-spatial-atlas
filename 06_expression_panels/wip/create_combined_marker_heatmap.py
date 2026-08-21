import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob

print("Creating heatmaps from combined marker analysis...")

# Read all the marker results files
marker_files = glob.glob("*_top10_markers_combined.csv")
condition_files = glob.glob("*_condition_differences.csv")

print(f"Found {len(marker_files)} marker files")
print(f"Found {len(condition_files)} condition difference files")

# === HEATMAP 1: Cluster Markers (Combined Analysis) ===
print("\n1. Creating cluster marker fold change heatmap...")

marker_data = []
all_genes = set()

for file in marker_files:
    # Extract cluster name from filename  
    cluster_name = file.replace("_top10_markers_combined.csv", "")
    cluster_name = cluster_name.replace("_", " ").replace(" ", ": ", 1)
    
    # Read data
    df = pd.read_csv(file)
    df['cluster_clean'] = cluster_name
    
    # Add genes to set
    all_genes.update(df['gene'].tolist())
    
    marker_data.append(df)
    print(f"  {cluster_name}: {len(df)} marker genes")

# Combine marker data
combined_markers = pd.concat(marker_data, ignore_index=True)
clusters = sorted(combined_markers['cluster_clean'].unique())

# Order genes by cluster (like before)
ordered_genes_list = []
for cluster in clusters:
    cluster_data = combined_markers[combined_markers['cluster_clean'] == cluster].copy()
    cluster_data = cluster_data.sort_values('fold_change', ascending=False)
    cluster_genes = cluster_data['gene'].tolist()
    ordered_genes_list.extend(cluster_genes)

# Remove duplicates while preserving order
seen = set()
all_genes_list = []
for gene in ordered_genes_list:
    if gene not in seen:
        all_genes_list.append(gene)
        seen.add(gene)

# Create fold change matrix
marker_fc_matrix = pd.DataFrame(index=clusters, columns=all_genes_list, dtype=float)

for _, row in combined_markers.iterrows():
    cluster = row['cluster_clean']
    gene = row['gene']
    fold_change = row['fold_change']
    marker_fc_matrix.loc[cluster, gene] = fold_change

# Fill NaN with 1 (no change)
marker_fc_matrix = marker_fc_matrix.fillna(1.0)
marker_log2_matrix = np.log2(marker_fc_matrix)

# Plot marker heatmap
plt.figure(figsize=(max(20, len(all_genes_list) * 0.3), max(10, len(clusters) * 0.5)))

sns.heatmap(marker_log2_matrix, 
            cmap='RdBu_r',
            center=0,
            annot=False,
            cbar_kws={'label': 'Log2 Fold Change (Cluster vs All Others)'},
            xticklabels=True,
            yticklabels=True)

# Add gene group separators
gene_boundaries = []
current_pos = 0
for cluster in clusters:
    cluster_data = combined_markers[combined_markers['cluster_clean'] == cluster]
    n_genes = len(cluster_data)
    current_pos += n_genes
    if current_pos < len(all_genes_list):
        gene_boundaries.append(current_pos)

for boundary in gene_boundaries:
    plt.axvline(x=boundary, color='white', linewidth=2, alpha=0.8)

plt.title('Cluster Marker Genes (Combined Analysis)\nEach cluster vs all others', 
          fontsize=16, pad=20)
plt.xlabel('Marker Genes (grouped by cluster)', fontsize=14)
plt.ylabel('Clusters', fontsize=14)
plt.xticks(rotation=45, ha='right')
plt.yticks(rotation=0)

plt.tight_layout()
plt.savefig('cluster_markers_heatmap_combined.png', dpi=300, bbox_inches='tight')
plt.savefig('cluster_markers_heatmap_combined.pdf', bbox_inches='tight')
plt.close()

# === HEATMAP 2: Condition Differences for Marker Genes ===
print("\n2. Creating condition difference heatmap...")

condition_data = []
for file in condition_files:
    cluster_name = file.replace("_condition_differences.csv", "")
    cluster_name = cluster_name.replace("_", " ").replace(" ", ": ", 1)
    
    df = pd.read_csv(file)
    df['cluster_clean'] = cluster_name
    condition_data.append(df)

# Combine condition data
combined_conditions = pd.concat(condition_data, ignore_index=True)

# Create condition difference matrix using same gene order
condition_fc_matrix = pd.DataFrame(index=clusters, columns=all_genes_list, dtype=float)

for _, row in combined_conditions.iterrows():
    cluster = row['cluster_clean']
    gene = row['gene']
    condition_fc = row['condition_fold_change']
    condition_fc_matrix.loc[cluster, gene] = condition_fc

# Fill NaN with 1 (no difference between conditions)
condition_fc_matrix = condition_fc_matrix.fillna(1.0)
condition_log2_matrix = np.log2(condition_fc_matrix)

# Plot condition difference heatmap
plt.figure(figsize=(max(20, len(all_genes_list) * 0.3), max(10, len(clusters) * 0.5)))

sns.heatmap(condition_log2_matrix, 
            cmap='RdBu_r',
            center=0,
            annot=False,
            cbar_kws={'label': 'Log2 Fold Change (Infected/Uninfected)'},
            xticklabels=True,
            yticklabels=True)

# Add gene group separators
for boundary in gene_boundaries:
    plt.axvline(x=boundary, color='white', linewidth=2, alpha=0.8)

plt.title('Condition Differences for Cluster Marker Genes\nInfected vs Uninfected expression', 
          fontsize=16, pad=20)
plt.xlabel('Marker Genes (grouped by cluster)', fontsize=14)
plt.ylabel('Clusters', fontsize=14)
plt.xticks(rotation=45, ha='right')
plt.yticks(rotation=0)

plt.tight_layout()
plt.savefig('condition_differences_heatmap.png', dpi=300, bbox_inches='tight')
plt.savefig('condition_differences_heatmap.pdf', bbox_inches='tight')
plt.close()

# === COMBINED HEATMAP: Side by side ===
print("\n3. Creating combined side-by-side heatmap...")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(max(40, len(all_genes_list) * 0.6), max(10, len(clusters) * 0.5)))

# Left: Marker fold changes
sns.heatmap(marker_log2_matrix, 
            cmap='RdBu_r',
            center=0,
            annot=False,
            cbar_kws={'label': 'Log2 FC (Cluster vs Others)'},
            xticklabels=True,
            yticklabels=True,
            ax=ax1)

for boundary in gene_boundaries:
    ax1.axvline(x=boundary, color='white', linewidth=2, alpha=0.8)

ax1.set_title('A) Cluster Markers\n(Combined Analysis)', fontsize=14)
ax1.set_xlabel('Marker Genes', fontsize=12)
ax1.set_ylabel('Clusters', fontsize=12)
ax1.tick_params(axis='x', rotation=45)

# Right: Condition differences
sns.heatmap(condition_log2_matrix, 
            cmap='RdBu_r',
            center=0,
            annot=False,
            cbar_kws={'label': 'Log2 FC (Inf/Uninf)'},
            xticklabels=True,
            yticklabels=True,
            ax=ax2)

for boundary in gene_boundaries:
    ax2.axvline(x=boundary, color='white', linewidth=2, alpha=0.8)

ax2.set_title('B) Condition Differences\n(Same Marker Genes)', fontsize=14)
ax2.set_xlabel('Marker Genes', fontsize=12)
ax2.set_ylabel('')
ax2.tick_params(axis='x', rotation=45)

plt.tight_layout()
plt.savefig('combined_marker_and_condition_heatmap.png', dpi=300, bbox_inches='tight')
plt.savefig('combined_marker_and_condition_heatmap.pdf', bbox_inches='tight')
plt.close()

# Save matrices
marker_fc_matrix.to_csv('cluster_marker_fold_changes.csv')
marker_log2_matrix.to_csv('cluster_marker_log2_fold_changes.csv')
condition_fc_matrix.to_csv('condition_difference_fold_changes.csv')
condition_log2_matrix.to_csv('condition_difference_log2_fold_changes.csv')

# Summary statistics
print("\nSummary Statistics:")
print(f"Marker fold change range: {marker_fc_matrix.min().min():.2f} to {marker_fc_matrix.max().max():.2f}")
print(f"Condition fold change range: {condition_fc_matrix.min().min():.2f} to {condition_fc_matrix.max().max():.2f}")

# Most infection-responsive marker genes
if len(combined_conditions) > 0:
    top_upregulated = combined_conditions.nlargest(10, 'condition_fold_change')[['gene', 'cluster', 'condition_fold_change']]
    top_downregulated = combined_conditions.nsmallest(10, 'condition_fold_change')[['gene', 'cluster', 'condition_fold_change']]
    
    print(f"\nTop 10 marker genes upregulated in infection:")
    print(top_upregulated.to_string(index=False))
    
    print(f"\nTop 10 marker genes downregulated in infection:")
    print(top_downregulated.to_string(index=False))

print(f"\nFiles generated:")
print("- cluster_markers_heatmap_combined.png/pdf")
print("- condition_differences_heatmap.png/pdf") 
print("- combined_marker_and_condition_heatmap.png/pdf")
print("- cluster_marker_fold_changes.csv")
print("- condition_difference_fold_changes.csv") 