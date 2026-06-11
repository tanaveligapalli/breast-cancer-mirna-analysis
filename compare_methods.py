import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

RESULTS_FOLDER = "/Users/tanaveligapalli/R/results"
PYDESEQ_PATH   = "/Users/tanaveligapalli/Downloads/miRNA_results/pydeseq2_results.csv"
EDGER_PATH     = f"{RESULTS_FOLDER}/edgeR_results.csv"

# ── LOAD BOTH RESULTS ─────────────────────────────────────────
print("Loading results...")
edger  = pd.read_csv(EDGER_PATH)
deseq  = pd.read_csv(PYDESEQ_PATH)

# Standardise column names to log2FC and fdr
edger  = edger.rename(columns={'log2FC': 'log2FC', 'fdr': 'fdr'})
deseq  = deseq.rename(columns={'log2FoldChange': 'log2FC', 'padj': 'fdr'})

# Keep only significant in each method
edger_sig = edger[edger['significant'] == True].copy()
deseq_sig = deseq[deseq['significant'] == True].copy()

print(f"  edgeR significant:   {len(edger_sig)}")
print(f"  PyDESeq2 significant: {len(deseq_sig)}")

# ── TOP 10 UPREGULATED ────────────────────────────────────────
# Sort by fold change descending — largest positive FC first
edger_top_up  = edger_sig[edger_sig['log2FC'] > 0].nlargest(10, 'log2FC')['mirna'].tolist()
deseq_top_up  = deseq_sig[deseq_sig['log2FC'] > 0].nlargest(10, 'log2FC')['mirna'].tolist()
consensus_up  = [m for m in edger_top_up if m in deseq_top_up]

# ── TOP 10 DOWNREGULATED ──────────────────────────────────────
# Sort by fold change ascending — largest negative FC first
edger_top_dn  = edger_sig[edger_sig['log2FC'] < 0].nsmallest(10, 'log2FC')['mirna'].tolist()
deseq_top_dn  = deseq_sig[deseq_sig['log2FC'] < 0].nsmallest(10, 'log2FC')['mirna'].tolist()
consensus_dn  = [m for m in edger_top_dn if m in deseq_top_dn]

# ── PRINT RESULTS ─────────────────────────────────────────────
print("\n" + "="*55)
print("TOP 10 UPREGULATED — edgeR")
print("="*55)
print(edger_sig[edger_sig['log2FC']>0].nlargest(10,'log2FC')[['mirna','log2FC','fdr']].to_string(index=False))

print("\n" + "="*55)
print("TOP 10 UPREGULATED — PyDESeq2")
print("="*55)
print(deseq_sig[deseq_sig['log2FC']>0].nlargest(10,'log2FC')[['mirna','log2FC','fdr']].to_string(index=False))

print("\n" + "="*55)
print(f"CONSENSUS UPREGULATED (in BOTH top 10): {len(consensus_up)}")
print("="*55)
print(consensus_up)

print("\n" + "="*55)
print("TOP 10 DOWNREGULATED — edgeR")
print("="*55)
print(edger_sig[edger_sig['log2FC']<0].nsmallest(10,'log2FC')[['mirna','log2FC','fdr']].to_string(index=False))

print("\n" + "="*55)
print("TOP 10 DOWNREGULATED — PyDESeq2")
print("="*55)
print(deseq_sig[deseq_sig['log2FC']<0].nsmallest(10,'log2FC')[['mirna','log2FC','fdr']].to_string(index=False))

print("\n" + "="*55)
print(f"CONSENSUS DOWNREGULATED (in BOTH top 10): {len(consensus_dn)}")
print("="*55)
print(consensus_dn)

# ── BUILD FULL CONSENSUS TABLE ────────────────────────────────
all_consensus = consensus_up + consensus_dn

# Merge edgeR and PyDESeq2 stats for consensus miRNAs
edger_sub = edger[edger['mirna'].isin(all_consensus)][['mirna','log2FC','fdr']].rename(
    columns={'log2FC':'log2FC_edgeR','fdr':'fdr_edgeR'})
deseq_sub = deseq[deseq['mirna'].isin(all_consensus)][['mirna','log2FC','fdr']].rename(
    columns={'log2FC':'log2FC_deseq','fdr':'fdr_deseq'})

consensus_df = edger_sub.merge(deseq_sub, on='mirna')
consensus_df['direction'] = consensus_df['log2FC_edgeR'].apply(
    lambda x: 'UP in cancer' if x > 0 else 'DOWN in cancer'
)
consensus_df['mean_log2FC'] = (
    (consensus_df['log2FC_edgeR'] + consensus_df['log2FC_deseq']) / 2
).round(4)
consensus_df = consensus_df.sort_values('mean_log2FC', ascending=False)

print("\n" + "="*55)
print("FULL CONSENSUS TABLE")
print("="*55)
print(consensus_df[['mirna','log2FC_edgeR','log2FC_deseq',
                     'mean_log2FC','direction']].to_string(index=False))

# ── SAVE ──────────────────────────────────────────────────────
consensus_df.to_csv(f"{RESULTS_FOLDER}/consensus_top_mirnas.csv", index=False)
print(f"\nSaved: {RESULTS_FOLDER}/consensus_top_mirnas.csv")

# ── PLOT ──────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle("edgeR vs PyDESeq2 — Top miRNA Comparison", fontsize=13, fontweight='bold')

# Panel 1: Fold change comparison for all significant miRNAs
ax = axes[0]
merged = edger[['mirna','log2FC']].merge(
    deseq[['mirna','log2FC']], on='mirna', suffixes=('_edgeR','_deseq')
)
is_consensus = merged['mirna'].isin(all_consensus)
ax.scatter(merged.loc[~is_consensus,'log2FC_edgeR'],
           merged.loc[~is_consensus,'log2FC_deseq'],
           s=6, alpha=0.3, c='#CCCCCC', label='other')
ax.scatter(merged.loc[is_consensus,'log2FC_edgeR'],
           merged.loc[is_consensus,'log2FC_deseq'],
           s=40, alpha=0.9, c='#E24B4A', label=f'consensus top ({len(all_consensus)})')
corr = float(merged[['log2FC_edgeR','log2FC_deseq']].corr().iloc[0,1])
lim  = max(merged['log2FC_edgeR'].abs().max(),
           merged['log2FC_deseq'].abs().max()) * 1.1
ax.plot([-lim,lim],[-lim,lim],'r--',linewidth=1,label='x=y')
for _, row in merged[is_consensus].iterrows():
    ax.annotate(row['mirna'], (row['log2FC_edgeR'], row['log2FC_deseq']),
                fontsize=5, ha='left')
ax.set_xlabel("edgeR log2FC")
ax.set_ylabel("PyDESeq2 log2FC")
ax.set_title(f"Fold change agreement\nr = {corr:.3f}")
ax.legend(fontsize=7)

# Panel 2: Bar chart of consensus miRNAs by mean fold change
ax2 = axes[1]
colors = ['#E24B4A' if x > 0 else '#378ADD'
          for x in consensus_df['mean_log2FC']]
ax2.barh(consensus_df['mirna'], consensus_df['mean_log2FC'],
         color=colors, alpha=0.85)
ax2.axvline(0, color='black', linewidth=0.8)
ax2.axvline( 0.5, color='gray', linestyle=':', linewidth=0.8)
ax2.axvline(-0.5, color='gray', linestyle=':', linewidth=0.8)
ax2.set_xlabel("Mean log2 Fold Change (Cancer vs Healthy)")
ax2.set_title("Consensus miRNAs\n(red = up in cancer, blue = down)")

plt.tight_layout()
plt.savefig(f"{RESULTS_FOLDER}/consensus_comparison.png", dpi=150, bbox_inches='tight')
plt.close()
print(f"Plot saved: {RESULTS_FOLDER}/consensus_comparison.png")

print(f"""
{'='*55}
SUMMARY
{'='*55}
  Consensus upregulated in cancer:   {len(consensus_up)}
  Consensus downregulated in cancer: {len(consensus_dn)}
  Total consensus biomarkers:        {len(all_consensus)}

  These miRNAs are significant in BOTH edgeR and PyDESeq2
  and have the largest fold changes in both methods.
  These are the targets your device should detect.
""")
