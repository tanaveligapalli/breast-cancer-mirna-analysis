# breast-cancer-mirna-analysis
Analysis of GSE73002 RAW data by first building an expression matrix and filtering out spike probes etc. and then conducting differential gene extraction using pyDEseq2 and edgeR to indicate strength of signal despite differences in statistical methodology

Methods: edgeR, PyDESeq2, Welch t-test.
Dataset: 1278 breast cancer vs 2740 healthy controls.
Key output: 10 consensus upregulated, 10 downregulated miRNAs.

For subtype analysis,
Since classical PAM50 subtype marker miRNAs were not detected in serum at sufficient levels, we identified data-driven molecular subgroups within the cancer cohort using the 153 miRNAs that were significantly differentially expressed between cancer and healthy controls
