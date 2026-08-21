import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob

print("Creating heatmaps from complete marker analysis...")

# Read all marker result files
shared_marker_files = glob.glob("*_top10_markers_combined.csv")
infected_only_marker_files = glob.glob("*_top10_markers_infected_only.csv")
condition_files = glob.glob("*_condition_differences.csv")

print(f"Found {len(shared_marker_files)} shared cluster marker files")
print(f"Found {len(infected_only_marker_files)} infected-only cluster marker files")
print(f"Found {len(condition_files)} condition difference files")

# === COLLECT ALL MARKER DATA ===
all_marker_data = []
all_genes = set()

# Process shared clusters
for file in shared_marker_files:
    cluster_name = file.replace("_top10_markers_combined.csv", "")
    cluster_name = cluster_name.replace("_", " ").replace(" ", ": ", 1)
    
    df = pd.read_csv(file)
    df['cluster_clean'] = cluster_name
    df['cluster_type'] = 'shared'
    
    all_genes.update(df['gene'].tolist())
    all_marker_data.append(df)
    print(f"  {cluster_name}: {len(df)} marker genes (shared)")

# Process infected-only clusters
for file in infected_only_marker_files:
    cluster_name = file.replace("_top10_markers_infected_only.csv", "")
    cluster_name = cluster_name.replace("_", " ").replace(" ", ": ", 1)
    
    df = pd.read_csv(file)
    df['cluster_clean'] = cluster_name
    df['cluster_type'] = 'infected_only'
    
    all_genes.update(df['gene'].tolist())
    all_marker_data.append(df)
    print(f"  {cluster_name}: {len(df)} marker genes (infected-only)")

# Combine all marker data
combined_markers = pd.concat(all_marker_data, ignore_index=True)
all_clusters = sorted(combined_markers['cluster_clean'].unique())

# Order genes by cluster (shared first, then infected-only)
shared_clusters = sorted([c for c in all_clusters if (combined_markers[combined_markers['cluster_clean'] == c]['cluster_type'] == 'shared').any()])
infected_only_clusters = sorted([c for c in all_clusters if (combined_markers[combined_markers['cluster_clean'] == c]['cluster_type'] == 'infected_only').any()])

ordered_clusters = shared_clusters + infected_only_clusters

print(f"\nCluster order: {len(shared_clusters)} shared + {len(infected_only_clusters)} infected-only")

# Order genes by cluster
ordered_genes_list = []
gene_boundaries = []
current_pos = 0

for cluster in ordered_clusters:
    cluster_data = combined_markers[combined_markers['cluster_clean'] == cluster].copy()
    cluster_data = cluster_data.sort_values('fold_change', ascending=False)
    cluster_genes = cluster_data['gene'].tolist()
    ordered_genes_list.extend(cluster_genes)
    
    current_pos += len(cluster_genes)
    if current_pos < len(ordered_genes_list):  # Will be updated after deduplication
        gene_boundaries.append(current_pos)

# Remove duplicates while preserving order
seen = set()
all_genes_list = []
for gene in ordered_genes_list:
    if gene not in seen:
        all_genes_list.append(gene)
        seen.add(gene)

# Recalculate gene boundaries after deduplication
gene_boundaries = []
current_pos = 0
for cluster in ordered_clusters:
    cluster_data = combined_markers[combined_markers['cluster_clean'] == cluster].copy()
    n_unique_genes = len([g for g in cluster_data['gene'].tolist() if g in all_genes_list])
    current_pos += n_unique_genes
    if current_pos < len(all_genes_list):
        gene_boundaries.append(current_pos)

print(f"Total unique genes: {len(all_genes_list)}")

# === HEATMAP 1: Cluster Markers ===
print("\n1. Creating cluster marker heatmap...")

marker_fc_matrix = pd.DataFrame(index=ordered_clusters, columns=all_genes_list, dtype=float)

for _, row in combined_markers.iterrows():
    cluster = row['cluster_clean']
    gene = row['gene']
    fold_change = row['fold_change']
    marker_fc_matrix.loc[cluster, gene] = fold_change

marker_fc_matrix = marker_fc_matrix.fillna(1.0)
marker_log2_matrix = np.log2(marker_fc_matrix)

# Plot marker heatmap
plt.figure(figsize=(max(25, len(all_genes_list) * 0.4), max(12, len(ordered_clusters) * 0.8)))

sns.heatmap(marker_log2_matrix, 
            cmap='RdBu_r',
            center=0,
            annot=False,
            cbar_kws={'label': 'Log2 Fold Change (Cluster vs All Others)'},
            xticklabels=True,
            yticklabels=True)

# Add separators
for boundary in gene_boundaries:
    plt.axvline(x=boundary, color='white', linewidth=2, alpha=0.8)

# Add horizontal line to separate shared from infected-only clusters
if len(infected_only_clusters) > 0:
    plt.axhline(y=len(shared_clusters), color='black', linewidth=3, alpha=0.8)

plt.title('Cluster Marker Genes: All Clusters\n(Shared clusters + Infected-only clusters)', 
          fontsize=18, pad=20)
plt.xlabel('Marker Genes (grouped by cluster)', fontsize=14)
plt.ylabel('Clusters', fontsize=14)
plt.xticks(rotation=45, ha='right')
plt.yticks(rotation=0)

# Add text annotations for cluster types
plt.text(-0.1, len(shared_clusters)/2, 'SHARED\nCLUSTERS', 
         rotation=90, va='center', ha='center', fontsize=12, weight='bold',
         transform=plt.gca().transAxes)

if len(infected_only_clusters) > 0:
    plt.text(-0.1, len(shared_clusters) + len(infected_only_clusters)/2, 'INFECTED-ONLY\nCLUSTERS', 
             rotation=90, va='center', ha='center', fontsize=12, weight='bold', color='red',
             transform=plt.gca().transAxes)

plt.tight_layout()
plt.savefig('complete_cluster_markers_heatmap.png', dpi=300, bbox_inches='tight')
plt.savefig('complete_cluster_markers_heatmap.pdf', bbox_inches='tight')
plt.close()

# === HEATMAP 2: Condition Differences (Shared Clusters Only) ===
print("\n2. Creating condition difference heatmap...")

condition_data = []
for file in condition_files:
    cluster_name = file.replace("_condition_differences.csv", "")
    cluster_name = cluster_name.replace("_", " ").replace(" ", ": ", 1)
    
    df = pd.read_csv(file)
    df['cluster_clean'] = cluster_name
    condition_data.append(df)

if condition_data:
    combined_conditions = pd.concat(condition_data, ignore_index=True)
    
    # Create condition difference matrix (only for shared clusters)
    condition_fc_matrix = pd.DataFrame(index=shared_clusters, columns=all_genes_list, dtype=float)
    
    for _, row in combined_conditions.iterrows():
        cluster = row['cluster_clean']
        gene = row['gene']
        condition_fc = row['condition_fold_change']
        condition_fc_matrix.loc[cluster, gene] = condition_fc
    
    condition_fc_matrix = condition_fc_matrix.fillna(1.0)
    condition_log2_matrix = np.log2(condition_fc_matrix)
    
    # Plot condition difference heatmap
    plt.figure(figsize=(max(25, len(all_genes_list) * 0.4), max(10, len(shared_clusters) * 0.8)))
    
    sns.heatmap(condition_log2_matrix, 
                cmap='RdBu_r',
                center=0,
                annot=False,
                cbar_kws={'label': 'Log2 Fold Change (Infected/Uninfected)'},
                xticklabels=True,
                yticklabels=True)
    
    # Add gene separators (only for shared clusters)
    shared_gene_boundaries = []
    current_pos = 0
    for cluster in shared_clusters:
        cluster_data = combined_markers[combined_markers['cluster_clean'] == cluster]
        n_unique_genes = len([g for g in cluster_data['gene'].tolist() if g in all_genes_list])
        current_pos += n_unique_genes
        if current_pos < len(all_genes_list):
            shared_gene_boundaries.append(current_pos)
    
    for boundary in shared_gene_boundaries:
        plt.axvline(x=boundary, color='white', linewidth=2, alpha=0.8)
    
    plt.title('Condition Differences for Cluster Marker Genes\n(Shared clusters only: Infected vs Uninfected)', 
              fontsize=18, pad=20)
    plt.xlabel('Marker Genes (grouped by cluster)', fontsize=14)
    plt.ylabel('Shared Clusters', fontsize=14)
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    
    plt.tight_layout()
    plt.savefig('condition_differences_heatmap_shared_only.png', dpi=300, bbox_inches='tight')
    plt.savefig('condition_differences_heatmap_shared_only.pdf', bbox_inches='tight')
    plt.close()

# === COMBINED HEATMAP: Side by side ===
print("\n3. Creating combined side-by-side heatmap...")

if condition_data:
    fig = plt.figure(figsize=(max(50, len(all_genes_list) * 0.8), max(12, len(ordered_clusters) * 0.8)))
    
    # Left: All cluster markers
    ax1 = plt.subplot(1, 2, 1)
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
    
    if len(infected_only_clusters) > 0:
        ax1.axhline(y=len(shared_clusters), color='black', linewidth=3, alpha=0.8)
    
    ax1.set_title('A) Cluster Markers\n(All clusters)', fontsize=16)
    ax1.set_xlabel('Marker Genes', fontsize=12)
    ax1.set_ylabel('Clusters', fontsize=12)
    ax1.tick_params(axis='x', rotation=45)
    
    # Right: Condition differences (shared only)
    ax2 = plt.subplot(1, 2, 2)
    sns.heatmap(condition_log2_matrix, 
                cmap='RdBu_r',
                center=0,
                annot=False,
                cbar_kws={'label': 'Log2 FC (Inf/Uninf)'},
                xticklabels=True,
                yticklabels=True,
                ax=ax2)
    
    for boundary in shared_gene_boundaries:
        ax2.axvline(x=boundary, color='white', linewidth=2, alpha=0.8)
    
    ax2.set_title('B) Condition Differences\n(Shared clusters only)', fontsize=16)
    ax2.set_xlabel('Marker Genes', fontsize=12)
    ax2.set_ylabel('Shared Clusters', fontsize=12)
    ax2.tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig('complete_combined_marker_and_condition_heatmap.png', dpi=300, bbox_inches='tight')
    plt.savefig('complete_combined_marker_and_condition_heatmap.pdf', bbox_inches='tight')
    plt.close()

# Save matrices
marker_fc_matrix.to_csv('complete_cluster_marker_fold_changes.csv')
marker_log2_matrix.to_csv('complete_cluster_marker_log2_fold_changes.csv')

if condition_data:
    condition_fc_matrix.to_csv('shared_clusters_condition_differences.csv')
    condition_log2_matrix.to_csv('shared_clusters_condition_log2_differences.csv')

# Summary statistics
print("\nSummary Statistics:")
print(f"Marker fold change range: {marker_fc_matrix.min().min():.2f} to {marker_fc_matrix.max().max():.2f}")

if condition_data:
    print(f"Condition fold change range: {condition_fc_matrix.min().min():.2f} to {condition_fc_matrix.max().max():.2f}")
    
    # Most infection-responsive marker genes (shared clusters)
    top_upregulated = combined_conditions.nlargest(10, 'condition_fold_change')[['gene', 'cluster', 'condition_fold_change']]
    top_downregulated = combined_conditions.nsmallest(10, 'condition_fold_change')[['gene', 'cluster', 'condition_fold_change']]
    
    print(f"\nTop 10 marker genes upregulated in infection (shared clusters):")
    print(top_upregulated.to_string(index=False))
    
    print(f"\nTop 10 marker genes downregulated in infection (shared clusters):")
    print(top_downregulated.to_string(index=False))

# Show infected-only cluster markers
if infected_only_clusters:
    print(f"\nInfected-only cluster markers:")
    for cluster in infected_only_clusters:
        cluster_data = combined_markers[combined_markers['cluster_clean'] == cluster]
        top_marker = cluster_data.nlargest(1, 'fold_change')
        if len(top_marker) > 0:
            print(f"  {cluster}: {top_marker.iloc[0]['gene']} (FC={top_marker.iloc[0]['fold_change']:.2f})")

print(f"\nFiles generated:")
print("- complete_cluster_markers_heatmap.png/pdf")
if condition_data:
    print("- condition_differences_heatmap_shared_only.png/pdf")
    print("- complete_combined_marker_and_condition_heatmap.png/pdf")
print("- complete_cluster_marker_fold_changes.csv")
if condition_data:
    print("- shared_clusters_condition_differences.csv") 