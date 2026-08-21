import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob
from pathlib import Path

# Read all the top10 genes files
print("Reading fold change results...")
top10_files = glob.glob("INFECTED_*_vs_ALL_top10_genes.csv")
print(f"Found {len(top10_files)} cluster result files")

# Collect all data
all_data = []
all_genes = set()

for file in top10_files:
    # Extract cluster name from filename
    cluster_name = file.replace("INFECTED_", "").replace("_vs_ALL_top10_genes.csv", "")
    cluster_name = cluster_name.replace("_", " ").replace(" ", ": ", 1)  # Convert back to original format
    
    # Read the data
    df = pd.read_csv(file)
    df['cluster_clean'] = cluster_name
    
    # Add genes to our set
    all_genes.update(df['gene'].tolist())
    
    all_data.append(df)
    print(f"  {cluster_name}: {len(df)} genes")

# Combine all data
combined_df = pd.concat(all_data, ignore_index=True)
print(f"\nTotal unique genes: {len(all_genes)}")
print(f"Total cluster-gene combinations: {len(combined_df)}")

# Create a matrix for the heatmap
print("\nCreating fold change matrix...")

# Get all unique clusters (sorted for consistent ordering)
clusters = sorted(combined_df['cluster_clean'].unique())
print(f"Clusters: {len(clusters)}")

# Get genes ordered by cluster (top genes from each cluster in sequence)
ordered_genes_list = []
for cluster in clusters:
    cluster_data = combined_df[combined_df['cluster_clean'] == cluster].copy()
    # Sort by fold change to get top genes in order
    cluster_data = cluster_data.sort_values('fold_change', ascending=False)
    cluster_genes = cluster_data['gene'].tolist()
    ordered_genes_list.extend(cluster_genes)
    print(f"  {cluster}: {cluster_genes}")

print(f"Total genes (ordered by cluster): {len(ordered_genes_list)}")

# Remove duplicates while preserving order
seen = set()
all_genes_list = []
for gene in ordered_genes_list:
    if gene not in seen:
        all_genes_list.append(gene)
        seen.add(gene)

print(f"Unique genes (deduplicated): {len(all_genes_list)}")

# Create the matrix
fold_change_matrix = pd.DataFrame(index=clusters, columns=all_genes_list, dtype=float)

# Fill the matrix
for _, row in combined_df.iterrows():
    cluster = row['cluster_clean']
    gene = row['gene']
    fold_change = row['fold_change']
    fold_change_matrix.loc[cluster, gene] = fold_change

# Fill NaN values with 1 (no change) for visualization
fold_change_matrix = fold_change_matrix.fillna(1.0)

print(f"Matrix shape: {fold_change_matrix.shape}")

# Convert to log2 for better visualization
log2_fold_change_matrix = np.log2(fold_change_matrix)

# Create the heatmap
print("\nCreating heatmap...")
plt.figure(figsize=(max(20, len(all_genes_list) * 0.3), max(10, len(clusters) * 0.5)))

# Create heatmap
sns.heatmap(log2_fold_change_matrix, 
            cmap='RdBu_r',  # Red for high, blue for low
            center=0,  # Center colormap at 0 (no change)
            annot=False,  # Don't annotate with values (too many)
            cbar_kws={'label': 'Log2 Fold Change'},
            xticklabels=True,
            yticklabels=True)

# Add vertical lines to separate gene groups by cluster
gene_boundaries = []
current_pos = 0
for cluster in clusters:
    cluster_data = combined_df[combined_df['cluster_clean'] == cluster]
    n_genes = len(cluster_data)
    current_pos += n_genes
    if current_pos < len(all_genes_list):  # Don't add line after last group
        gene_boundaries.append(current_pos)

for boundary in gene_boundaries:
    plt.axvline(x=boundary, color='white', linewidth=2, alpha=0.8)

plt.title('Fold Change Heatmap: Infected Clusters vs All Others\n(Genes grouped by cluster, top 10 per cluster)', 
          fontsize=16, pad=20)
plt.xlabel('Genes (grouped by cluster)', fontsize=14)
plt.ylabel('Infected Clusters', fontsize=14)

# Rotate x-axis labels for readability
plt.xticks(rotation=45, ha='right')
plt.yticks(rotation=0)

plt.tight_layout()
plt.savefig('cluster_fold_change_heatmap.png', dpi=300, bbox_inches='tight')
plt.savefig('cluster_fold_change_heatmap.pdf', bbox_inches='tight')
plt.show()

# Save the matrix for reference
fold_change_matrix.to_csv('fold_change_matrix.csv')
log2_fold_change_matrix.to_csv('log2_fold_change_matrix.csv')

# Print some summary statistics
print("\nSummary Statistics:")
print(f"Fold change range: {fold_change_matrix.min().min():.2f} to {fold_change_matrix.max().max():.2f}")
print(f"Log2 fold change range: {log2_fold_change_matrix.min().min():.2f} to {log2_fold_change_matrix.max().max():.2f}")

# Show which genes appear most frequently across clusters
gene_counts = combined_df['gene'].value_counts()
print(f"\nMost frequent genes across clusters:")
print(gene_counts.head(10))

# Show cluster with highest fold changes
max_fc_per_cluster = combined_df.groupby('cluster_clean')['fold_change'].max().sort_values(ascending=False)
print(f"\nHighest fold changes per cluster:")
print(max_fc_per_cluster)

print(f"\nFiles generated:")
print("- cluster_fold_change_heatmap.png")
print("- cluster_fold_change_heatmap.pdf") 
print("- fold_change_matrix.csv")
print("- log2_fold_change_matrix.csv") 