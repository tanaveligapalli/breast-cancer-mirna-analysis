"""
Breast Cancer Subtype Prediction — Stage 1 & 2
================================================
Stage 1: Unsupervised clustering of cancer patients
Stage 2: Label clusters using known subtype miRNA markers
Stage 3: Save predicted subtype labels for differential expression

Input:  final_dataset.csv
Output: subtype_labels.csv        — predicted subtype per patient
        clustering_plots.png      — PCA + silhouette + heatmap
        subtype_summary.csv       — cluster statistics
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.metrics import silhouette_score, silhouette_samples
import warnings
warnings.filterwarnings('ignore')

DATA_PATH     = "/Users/tanaveligapalli/Downloads/miRNA_results/final_dataset.csv"
OUTPUT_FOLDER = "/Users/tanaveligapalli/Downloads/miRNA_results"

# ── LOAD DATA ─────────────────────────────────────────────────
print("Loading data...")
df         = pd.read_csv(DATA_PATH, index_col=0)
mirna_cols = [c for c in df.columns if c.startswith('hsa-')]

# Keep ONLY cancer patients for subtype clustering
# Healthy patients have no subtype — including them would distort clusters
cancer_df  = df[df['diagnosis'] == 'breast cancer'].copy()
expr       = cancer_df[mirna_cols].astype(float)

print(f"  Cancer patients for clustering: {len(expr)}")
print(f"  miRNA features: {len(mirna_cols)}")

# ── STEP 1: SCALE THE DATA ────────────────────────────────────
# StandardScaler: subtract mean, divide by std for each miRNA
# This ensures miRNAs with large absolute values don't dominate
# clustering just because of their scale
scaler      = StandardScaler()
expr_scaled = scaler.fit_transform(expr)

# ── STEP 2: PCA FOR DIMENSIONALITY REDUCTION ─────────────────
# 545 miRNAs is too many dimensions for clustering to work well
# PCA reduces to the most informative components
# We keep enough components to explain 80% of variance
pca        = PCA(n_components=0.80, random_state=42)
expr_pca   = pca.fit_transform(expr_scaled)
n_pca      = expr_pca.shape[1]
var_explained = pca.explained_variance_ratio_.cumsum()[-1]

print(f"\nPCA: {n_pca} components explain {var_explained*100:.1f}% of variance")

# Also keep 2D for visualisation
pca_2d     = PCA(n_components=2, random_state=42)
coords_2d  = pca_2d.fit_transform(expr_scaled)

# ── STEP 3: FIND OPTIMAL NUMBER OF CLUSTERS ──────────────────
# Test k = 2, 3, 4, 5 clusters
# Silhouette score measures cluster quality: 
#   +1 = well separated, 0 = overlapping, -1 = wrong cluster
print("\nFinding optimal number of clusters...")
print(f"  {'k':>3}  {'silhouette':>12}  {'inertia':>12}")

silhouette_scores = {}
inertia_scores    = {}

for k in range(2, 6):
    km     = KMeans(n_clusters=k, random_state=42, n_init=20)
    labels = km.fit_predict(expr_pca)
    sil    = silhouette_score(expr_pca, labels)
    silhouette_scores[k] = sil
    inertia_scores[k]    = km.inertia_
    print(f"  {k:>3}  {sil:>12.4f}  {km.inertia_:>12.1f}")

best_k = max(silhouette_scores, key=silhouette_scores.get)
print(f"\n  Best k = {best_k} (silhouette = {silhouette_scores[best_k]:.4f})")

# Force k=4 if best_k < 4 — we expect 4 biological subtypes
# If silhouette at k=4 is within 10% of best, use k=4 for interpretability
if best_k < 4:
    sil_4 = silhouette_scores.get(4, 0)
    sil_best = silhouette_scores[best_k]
    if sil_4 >= sil_best * 0.9:
        print(f"  Using k=4 (biologically motivated, silhouette={sil_4:.4f})")
        best_k = 4
    else:
        print(f"  k={best_k} has much better silhouette than k=4")
        print(f"  Proceeding with k={best_k} but biological interpretation may differ")

# ── STEP 4: FINAL CLUSTERING ──────────────────────────────────
print(f"\nRunning final clustering with k={best_k}...")
km_final  = KMeans(n_clusters=best_k, random_state=42, n_init=50)
clusters  = km_final.fit_predict(expr_pca)
cluster_labels = pd.Series(clusters, index=expr.index)

for c in range(best_k):
    n = (clusters == c).sum()
    print(f"  Cluster {c}: {n} patients ({n/len(clusters)*100:.1f}%)")

# ── STEP 5: LABEL CLUSTERS USING KNOWN SUBTYPE MARKERS ───────
# Known serum/blood miRNA markers for each breast cancer subtype
# From published literature — these are the miRNAs most consistently
# associated with each subtype in blood-based studies

SUBTYPE_MARKERS = {
    'Luminal_A': {
        'up'  : ['hsa-miR-100-5p', 'hsa-miR-125b-5p', 'hsa-miR-145-5p',
                 'hsa-let-7a-5p', 'hsa-let-7b-5p', 'hsa-miR-99a-5p'],
        'down': ['hsa-miR-21-5p', 'hsa-miR-155-5p', 'hsa-miR-10b-5p']
    },
    'Luminal_B': {
        'up'  : ['hsa-miR-21-5p', 'hsa-miR-210-3p', 'hsa-miR-196a-5p'],
        'down': ['hsa-miR-99a-5p', 'hsa-miR-100-5p', 'hsa-miR-125b-5p',
                 'hsa-let-7c-5p']
    },
    'HER2': {
        'up'  : ['hsa-miR-4728-3p', 'hsa-miR-21-5p', 'hsa-miR-141-3p',
                 'hsa-miR-210-3p', 'hsa-miR-375'],
        'down': ['hsa-miR-125b-5p', 'hsa-miR-145-5p']
    },
    'TNBC': {
        'up'  : ['hsa-miR-155-5p', 'hsa-miR-10b-5p', 'hsa-miR-17-5p',
                 'hsa-miR-19a-3p', 'hsa-miR-93-5p'],
        'down': ['hsa-miR-145-5p', 'hsa-miR-125b-5p', 'hsa-miR-100-5p']
    }
}

print("\nScoring clusters against subtype markers...")

# For each cluster, compute a score for each subtype
# Score = mean expression of upregulated markers - mean of downregulated
# Higher score = cluster looks more like that subtype

cluster_means = {}
for c in range(best_k):
    cluster_patients = expr.index[clusters == c]
    cluster_means[c] = expr.loc[cluster_patients].mean()

# Score each cluster against each subtype
subtype_scores = pd.DataFrame(index=range(best_k),
                               columns=SUBTYPE_MARKERS.keys(),
                               dtype=float)

for c in range(best_k):
    for subtype, markers in SUBTYPE_MARKERS.items():
        # Get markers that actually exist in our dataset
        up_markers   = [m for m in markers['up']   if m in expr.columns]
        down_markers = [m for m in markers['down'] if m in expr.columns]

        up_score   = cluster_means[c][up_markers].mean()   if up_markers   else 0
        down_score = cluster_means[c][down_markers].mean() if down_markers else 0

        # High up + low down = better match for this subtype
        subtype_scores.loc[c, subtype] = up_score - down_score

print("\nSubtype scores per cluster (higher = better match):")
print(subtype_scores.round(3).to_string())

# Assign each cluster to its best-matching subtype
# Use Hungarian-style greedy assignment to avoid duplicates
assigned_subtypes = {}
available_subtypes = list(SUBTYPE_MARKERS.keys())
available_clusters = list(range(best_k))

# Sort by score confidence (difference between top 2 scores)
for _ in range(min(best_k, len(available_subtypes))):
    best_score  = -np.inf
    best_c      = None
    best_sub    = None

    for c in available_clusters:
        for sub in available_subtypes:
            score = float(subtype_scores.loc[c, sub])
            if score > best_score:
                best_score = score
                best_c     = c
                best_sub   = sub

    if best_c is not None:
        assigned_subtypes[best_c] = best_sub
        available_clusters.remove(best_c)
        available_subtypes.remove(best_sub)

# Any remaining clusters get labelled by remaining subtypes or Unknown
for c in available_clusters:
    if available_subtypes:
        assigned_subtypes[c] = available_subtypes.pop(0)
    else:
        assigned_subtypes[c] = f'Unknown_{c}'

print("\nCluster → Subtype assignments:")
for c, sub in sorted(assigned_subtypes.items()):
    n = (clusters == c).sum()
    print(f"  Cluster {c} → {sub} ({n} patients)")

# Map cluster numbers to subtype names
subtype_map     = {idx: assigned_subtypes[c]
                   for idx, c in zip(expr.index, clusters)}
predicted_labels = pd.Series(subtype_map, name='predicted_subtype')

# ── STEP 6: SAVE SUBTYPE LABELS ───────────────────────────────
# Merge predicted subtypes back into the full dataset
cancer_df['predicted_subtype'] = predicted_labels

# Save labels file
labels_df = cancer_df[['diagnosis', 'label', 'predicted_subtype']].copy()
labels_df.to_csv(f"{OUTPUT_FOLDER}/subtype_labels.csv")
print(f"\nSaved: subtype_labels.csv")

# Summary statistics
print("\n" + "="*55)
print("SUBTYPE DISTRIBUTION")
print("="*55)
print(predicted_labels.value_counts().to_string())

# ── STEP 7: PLOTS ─────────────────────────────────────────────
print("\nGenerating plots...")

SUBTYPE_COLOURS = {
    'Luminal_A': '#378ADD',
    'Luminal_B': '#639922',
    'HER2'     : '#E24B4A',
    'TNBC'     : '#7F77DD',
    'Unknown_0': '#AAAAAA',
    'Unknown_1': '#CCCCCC',
}

fig, axes = plt.subplots(1, 3, figsize=(18, 6))
fig.suptitle("Breast Cancer Subtype Clustering", fontsize=13, fontweight='bold')

# Panel 1: PCA coloured by predicted subtype
ax = axes[0]
for subtype in predicted_labels.unique():
    mask = predicted_labels == subtype
    pts  = coords_2d[mask.reindex(expr.index).fillna(False).values]
    col  = SUBTYPE_COLOURS.get(subtype, '#AAAAAA')
    ax.scatter(pts[:, 0], pts[:, 1], s=8, alpha=0.6,
               color=col, label=f"{subtype} (n={mask.sum()})")
ax.set_xlabel(f"PC1 ({pca_2d.explained_variance_ratio_[0]*100:.1f}%)")
ax.set_ylabel(f"PC2 ({pca_2d.explained_variance_ratio_[1]*100:.1f}%)")
ax.set_title("PCA — predicted subtypes\n(separation = subtype signal present)")
ax.legend(fontsize=7, markerscale=2)

# Panel 2: Silhouette scores for different k
ax2 = axes[1]
ks   = list(silhouette_scores.keys())
sils = list(silhouette_scores.values())
cols = ['#E24B4A' if k == best_k else '#378ADD' for k in ks]
ax2.bar([str(k) for k in ks], sils, color=cols, alpha=0.85)
ax2.set_xlabel("Number of clusters (k)")
ax2.set_ylabel("Silhouette score")
ax2.set_title("Optimal cluster number\n(red = chosen, higher = better separation)")
for k, s in zip(ks, sils):
    ax2.text(str(k), s + 0.005, f"{s:.3f}", ha='center', fontsize=9)

# Panel 3: Subtype marker expression heatmap
ax3 = axes[2]
# Get all markers that exist in dataset
all_markers = []
for markers in SUBTYPE_MARKERS.values():
    all_markers.extend(markers['up'] + markers['down'])
all_markers = list(dict.fromkeys(  # deduplicate preserving order
    [m for m in all_markers if m in expr.columns]
))

if len(all_markers) > 0:
    # Mean expression per subtype for marker miRNAs
    subtype_expr = pd.DataFrame(index=predicted_labels.unique(),
                                columns=all_markers)
    for sub in predicted_labels.unique():
        patients = predicted_labels[predicted_labels == sub].index
        subtype_expr.loc[sub] = expr.loc[patients, all_markers].mean()

    subtype_expr = subtype_expr.astype(float)
    # Normalise each column to 0-1 for colour comparison
    norm = (subtype_expr - subtype_expr.min()) / (subtype_expr.max() - subtype_expr.min() + 1e-10)

    im = ax3.imshow(norm.values, aspect='auto', cmap='RdBu_r', vmin=0, vmax=1)
    ax3.set_xticks(range(len(all_markers)))
    ax3.set_xticklabels([m.replace('hsa-','') for m in all_markers],
                         rotation=90, fontsize=6)
    ax3.set_yticks(range(len(subtype_expr)))
    ax3.set_yticklabels(subtype_expr.index, fontsize=8)
    ax3.set_title("Marker miRNA expression\nper predicted subtype")
    plt.colorbar(im, ax=ax3, label='Relative expression')
else:
    ax3.text(0.5, 0.5, "No marker miRNAs\nfound in dataset",
             ha='center', va='center', transform=ax3.transAxes)
    ax3.set_title("Marker expression")

plt.tight_layout()
plt.savefig(f"{OUTPUT_FOLDER}/clustering_plots.png", dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved: clustering_plots.png")

print(f"""
{'='*55}
DONE — STAGE 1 & 2 COMPLETE
{'='*55}
subtype_labels.csv      — predicted subtype per patient
clustering_plots.png    — PCA + silhouette + heatmap

NEXT STEP — Stage 3:
  Run subtype_differential_expression.py to find
  which miRNAs distinguish each subtype from healthy
  and from each other subtype.

IMPORTANT NOTE:
  These are PREDICTED subtypes from serum miRNA clustering.
  They are not confirmed by pathology.
  State this clearly in your report as exploratory analysis.
""")
