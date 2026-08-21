import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import glob

print("Creating expression heatmap with STATISTICALLY SIGNIFICANT cluster-condition combinations...")

# Load the data
print("Loading data...")
adata_infected = sc.read_h5ad('cell2location_infected_030525_yk.h5ad')
adata_uninfected = sc.read_h5ad('cell2location_uninfected_030525_yk.h5ad')

# Get cluster names directly from data
infected_clusters = set(adata_infected.obs['paper_clusters'].unique())
uninfected_clusters = set(adata_uninfected.obs['paper_clusters'].unique())

# Identify shared and infected-only clusters  
shared_clusters = sorted(list(infected_clusters & uninfected_clusters))
infected_only_clusters = sorted(list(infected_clusters - uninfected_clusters))

print(f"Shared clusters: {len(shared_clusters)}")
print(f"Infected-only clusters: {len(infected_only_clusters)}")

# Read STATISTICAL marker files to get genes grouped by cluster
shared_marker_files = glob.glob("statistical_*_top10_markers_combined.csv")
infected_only_marker_files = glob.glob("statistical_*_top10_markers_infected_only.csv")
all_marker_files = shared_marker_files + infected_only_marker_files

print(f"Found {len(shared_marker_files)} shared cluster marker files")
print(f"Found {len(infected_only_marker_files)} infected-only cluster marker files")
print(f"Total marker files: {len(all_marker_files)}")

# Create mapping from filename to actual cluster names
filename_to_cluster = {
    'statistical_RP_A__GZMK': 'RP-A: GZMK',
    'statistical_RP_B__H_MK': 'RP-B: H_MK',
    'statistical_RP_C__NGP': 'RP-C: NGP', 
    'statistical_RP_D__F480': 'RP-D: F480',
    'statistical_RP_E__Rhag': 'RP-E: Rhag',
    'statistical_WP_A__TZ': 'WP-A: TZ',
    'statistical_WP_B__BZ': 'WP-B: BZ',
    'statistical_WP_C__MZ': 'WP-C: MZ',
    'statistical_WP_D__GC': 'WP-D: GC'
}

# Collect genes grouped by cluster and handle duplicates
cluster_genes = {}
all_gene_info = {}  # To track highest fold change for each gene

for file in all_marker_files:
    df = pd.read_csv(file)
    
    # Extract cluster name from filename
    base_filename = file.replace("_top10_markers_combined.csv", "").replace("_top10_markers_infected_only.csv", "")
    
    if base_filename in filename_to_cluster:
        cluster_name = filename_to_cluster[base_filename]
        print(f"  Processing {file} -> {cluster_name} ({len(df)} genes)")
        
        # Store genes for this cluster
        cluster_genes[cluster_name] = df['gene'].tolist()
        
        # Track the best fold change for each gene
        for _, row in df.iterrows():
            gene = row['gene']
            fc = row['fold_change']
            
            if gene not in all_gene_info or fc > all_gene_info[gene]['fold_change']:
                all_gene_info[gene] = {
                    'fold_change': fc,
                    'primary_cluster': cluster_name,
                    'padj': row.get('padj', 0.0)
                }
    else:
        print(f"  Warning: Could not map filename {file} to cluster name")

print(f"Successfully mapped {len(cluster_genes)} clusters")

# Create ordered gene list grouped by cluster, removing duplicates
final_genes = []
gene_boundaries = []  # For visual separators

# Process clusters in order
cluster_order = sorted(shared_clusters) + sorted(infected_only_clusters)

for cluster in cluster_order:
    cluster_specific_genes = []
    
    if cluster in cluster_genes:
        for gene in cluster_genes[cluster]:
            # Only include if this cluster is the primary cluster for this gene
            if all_gene_info[gene]['primary_cluster'] == cluster:
                if gene not in final_genes:  # Extra safety check
                    cluster_specific_genes.append(gene)
                    final_genes.append(gene)
    
    if cluster_specific_genes:
        gene_boundaries.append(len(final_genes) - 0.5)
        print(f"{cluster}: {len(cluster_specific_genes)} unique significant genes")

# Remove the last boundary (no line after last cluster)
gene_boundaries = gene_boundaries[:-1]

# Filter to genes present in both datasets
genes_in_both = set(adata_infected.var_names) & set(adata_uninfected.var_names)
final_genes = [gene for gene in final_genes if gene in genes_in_both]

print(f"Total unique statistically significant marker genes: {len(final_genes)}")
print(f"Genes available in both datasets: {len(final_genes)}")

# Create cluster-condition combinations
cluster_condition_combinations = []

# Add shared clusters (both infected and uninfected)
for cluster in shared_clusters:
    cluster_condition_combinations.append(f"{cluster} (Infected)")
    cluster_condition_combinations.append(f"{cluster} (Uninfected)")

# Add infected-only clusters  
for cluster in infected_only_clusters:
    cluster_condition_combinations.append(f"{cluster} (Infected)")

print(f"Total cluster-condition combinations: {len(cluster_condition_combinations)}")

# Calculate expression matrix
print("Calculating mean expression for each cluster-condition combination...")
expression_matrix = pd.DataFrame(index=cluster_condition_combinations, columns=final_genes, dtype=float)

for combo in cluster_condition_combinations:
    if "(Infected)" in combo:
        cluster_name = combo.replace(" (Infected)", "")
        if cluster_name in adata_infected.obs['paper_clusters'].values:
            cluster_mask = adata_infected.obs['paper_clusters'] == cluster_name
            cluster_spots = cluster_mask.sum()
            
            if cluster_spots > 0:
                cluster_expr = adata_infected[cluster_mask, final_genes].X
                if hasattr(cluster_expr, 'toarray'):
                    cluster_expr = cluster_expr.toarray()
                
                mean_expr = np.mean(cluster_expr, axis=0)
                expression_matrix.loc[combo, final_genes] = mean_expr
                print(f"  {combo}: {cluster_spots} spots")
    
    elif "(Uninfected)" in combo:
        cluster_name = combo.replace(" (Uninfected)", "")
        if cluster_name in adata_uninfected.obs['paper_clusters'].values:
            cluster_mask = adata_uninfected.obs['paper_clusters'] == cluster_name
            cluster_spots = cluster_mask.sum()
            
            if cluster_spots > 0:
                cluster_expr = adata_uninfected[cluster_mask, final_genes].X
                if hasattr(cluster_expr, 'toarray'):
                    cluster_expr = cluster_expr.toarray()
                
                mean_expr = np.mean(cluster_expr, axis=0)
                expression_matrix.loc[combo, final_genes] = mean_expr
                print(f"  {combo}: {cluster_spots} spots")

# Fill any NaN values with 0
expression_matrix = expression_matrix.fillna(0)

print(f"Final expression matrix shape: {expression_matrix.shape}")

# Normalize expression (z-score across all samples)
from scipy.stats import zscore
normalized_matrix = expression_matrix.apply(zscore, axis=0)

print("Creating heatmap...")

# Create the heatmap with much wider figure for gene readability
fig_width = max(35, len(final_genes) * 0.6)  # Even wider for statistical results
plt.figure(figsize=(fig_width, 14))

ax = sns.heatmap(normalized_matrix, 
                 cmap='RdBu_r', 
                 center=0,
                 cbar_kws={'label': 'Z-score normalized expression'},
                 xticklabels=True,
                 yticklabels=True)

# Add vertical lines to separate gene clusters
for boundary in gene_boundaries:
    plt.axvline(x=boundary, color='white', linewidth=2, alpha=0.8)

plt.title('Statistically Significant Cluster-Condition Expression Heatmap\n' + 
          '(Z-score normalized expression, genes with padj < 0.05, ranked by fold change)', 
          fontsize=16, pad=20)
plt.xlabel('Statistically Significant Marker Genes (grouped by primary cluster)', fontsize=12)
plt.ylabel('Cluster-Condition Combinations', fontsize=12)
plt.xticks(rotation=45, ha='right', fontsize=9)
plt.yticks(rotation=0, fontsize=10)

plt.tight_layout()
plt.savefig('statistical_cluster_condition_expression_heatmap.png', dpi=300, bbox_inches='tight')
plt.savefig('statistical_cluster_condition_expression_heatmap.pdf', bbox_inches='tight')

print("Heatmap saved as 'statistical_cluster_condition_expression_heatmap.png' and '.pdf'")

# Save the expression matrix
normalized_matrix.to_csv('statistical_cluster_condition_expression_matrix.csv')
print("Expression matrix saved as 'statistical_cluster_condition_expression_matrix.csv'")

# Create summary statistics
print(f"\nSummary of Statistical Analysis:")
print(f"- Total significant genes found: {len(final_genes)}")
print(f"- Average genes per cluster: {len(final_genes) / len(cluster_genes):.1f}")
print(f"- Significance threshold: padj < 0.05")
print(f"- Specificity filter: ≥10% difference in expression frequency")
print(f"- Ranking method: Fold change among significant genes")

# Show gene distribution by cluster
print(f"\nGenes per cluster:")
for cluster in cluster_order:
    if cluster in cluster_genes:
        n_genes = len([g for g in cluster_genes[cluster] if g in final_genes])
        print(f"  {cluster}: {n_genes} genes")

print("Analysis complete!")
print("\nKey differences from simple fold-change analysis:")
print("1. Statistical significance testing (Wilcoxon rank-sum)")
print("2. Multiple testing correction (FDR)")
print("3. Specificity filtering (expression frequency difference)")
print("4. Adaptive pseudocount based on data distribution")
print("5. Only significant genes included in final ranking") 