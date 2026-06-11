library(edgeR)

OUTPUT_FOLDER <- "/Users/tanaveligapalli/R/results"
MATRIX_PATH   <- "/Users/tanaveligapalli/Downloads/miRNA_results/final_dataset.csv"
# ── STEP 1: LOAD DATA ─────────────────────────────────────────
cat("Loading data...\n")

full_data <- read.csv(MATRIX_PATH, row.names = 1, check.names = FALSE)

cat("  Full dataset:", nrow(full_data), "patients x", ncol(full_data), "columns\n")

# Separate clinical columns from miRNA columns
mirna_cols <- grep("^hsa-", colnames(full_data), value = TRUE)
expr       <- full_data[, mirna_cols]
metadata   <- full_data[, c("diagnosis", "label")]

cat("  Expression:", nrow(expr), "patients x", ncol(expr), "miRNAs\n")
cat("  Diagnosis breakdown:\n")
print(table(metadata$diagnosis))

# ── STEP 2: REMOVE PROSTATE DISEASE ───────────────────────────
keep     <- metadata$diagnosis != "prostate disease"
expr     <- expr[keep, ]
metadata <- metadata[keep, ]

cat("  After removing prostate disease:", nrow(expr), "patients\n")

# ── STEP 3: BUILD GROUP LABELS ────────────────────────────────
group <- factor(
  ifelse(metadata$diagnosis == "breast cancer", "cancer", "healthy"),
  levels = c("healthy", "cancer")
)
cat("\n  Group counts:\n")
print(table(group))
 
# ── STEP 4: REVERSE LOG2 TRANSFORM → INTEGER COUNTS ──────────
# edgeR needs count data (non-negative integers)
# We reverse: counts = round(2^log2_value), floor at 1
# Same approach as PyDESeq2 — necessary approximation for microarray data
cat("\nPreparing count matrix...\n")
counts      <- round(2^as.matrix(expr))
counts      <- pmax(counts, 1)          # floor at 1, no zeros
counts      <- t(counts)                # edgeR wants miRNAs as rows, patients as cols
 
cat("  Count matrix:", nrow(counts), "miRNAs x", ncol(counts), "patients\n")
cat("  Value range:", min(counts), "to", max(counts), "\n")
 
# ── STEP 5: BUILD DGEList OBJECT ──────────────────────────────
# DGEList is edgeR's core data container
# It holds counts + group labels + normalisation factors
y <- DGEList(
  counts = counts,
  group  = group
)
 
# ── STEP 6: TMM NORMALISATION ─────────────────────────────────
# TMM = Trimmed Mean of M-values
# Scales each sample so that the majority of miRNAs appear unchanged
# More robust than DESeq2's median of ratios for miRNA data because
# it trims extreme fold changes before computing the scaling factor,
# preventing highly expressed outlier miRNAs from skewing normalisation
cat("Applying TMM normalisation...\n")
y <- calcNormFactors(y, method = "TMM")
 
cat("  Normalisation factors (should be close to 1.0):\n")
cat("  Min:", round(min(y$samples$norm.factors), 4),
    "Max:", round(max(y$samples$norm.factors), 4),
    "Mean:", round(mean(y$samples$norm.factors), 4), "\n")
 
# ── STEP 7: DISPERSION ESTIMATION ────────────────────────────
# edgeR estimates dispersion in three stages:
#
# Common dispersion: one shared value for all miRNAs
#   → gives a quick global picture of overdispersion
#
# Trended dispersion: smooth curve of dispersion vs mean expression
#   → accounts for the fact that lowly expressed miRNAs are more
#     variable than highly expressed ones
#
# Tagwise dispersion: per-miRNA dispersion shrunk toward the trend
#   → the final estimate used in testing
#   → robust=TRUE makes it resistant to outlier miRNAs distorting
#     the trend fit
 
cat("Estimating dispersions...\n")
y <- estimateDisp(y, robust = TRUE)
 
cat("  Common dispersion (BCV):", round(sqrt(y$common.dispersion), 4), "\n")
cat("  (BCV = Biological Coefficient of Variation)\n")
cat("  (0.1-0.4 is typical for human studies)\n")
 
# ── STEP 8: EXACT TEST ────────────────────────────────────────
# For a simple two-group comparison, the exact test is more powerful
# than the GLM Wald test that DESeq2 uses.
# It conditions on the total count for each miRNA and computes
# exact p-values from the negative binomial distribution.
# pair = c("healthy", "cancer") means: cancer relative to healthy
# so positive logFC = higher in cancer
 
cat("Running exact test...\n")
et <- exactTest(y, pair = c("healthy", "cancer"))
 
# Extract all results (n=Inf means return all miRNAs, not just top ones)
results        <- topTags(et, n = Inf, sort.by = "PValue")$table
results$mirna  <- rownames(results)
results        <- results[, c("mirna", "logFC", "logCPM", "PValue", "FDR")]
 
# Rename columns to match Python output convention
colnames(results) <- c("mirna", "log2FC", "logCPM", "p_value", "fdr")
 
# Add significance flag: FDR < 0.05 AND |log2FC| > 0.5
results$significant <- results$fdr < 0.05 & abs(results$log2FC) > 0.5
 
# ── STEP 9: PRINT SUMMARY ────────────────────────────────────
n_sig <- sum(results$significant)
n_up  <- sum(results$significant & results$log2FC > 0)
n_dn  <- sum(results$significant & results$log2FC < 0)
 
cat("\n==============================================\n")
cat("RESULTS\n")
cat("==============================================\n")
cat("Significant miRNAs (FDR<0.05, |logFC|>0.5):", n_sig, "\n")
cat("  Upregulated in cancer:  ", n_up, "\n")
cat("  Downregulated in cancer:", n_dn, "\n")
 
cat("\nTop 10 upregulated in cancer:\n")
top_up <- head(results[results$log2FC > 0 & results$significant, ], 10)
print(top_up[, c("mirna","log2FC","p_value","fdr")], row.names = FALSE)
 
cat("\nTop 10 downregulated in cancer:\n")
top_dn <- head(results[results$log2FC < 0 & results$significant, ], 10)
print(top_dn[, c("mirna","log2FC","p_value","fdr")], row.names = FALSE)
 
# ── STEP 10: SAVE RESULTS ────────────────────────────────────
out_path <- file.path(OUTPUT_FOLDER, "edgeR_results.csv")
write.csv(results, out_path, row.names = FALSE)
cat("\nSaved:", out_path, "\n")
 
# ── STEP 11: PLOTS ────────────────────────────────────────────
pdf(file.path(OUTPUT_FOLDER, "edgeR_plots.pdf"), width = 12, height = 5)
par(mfrow = c(1, 3))
 
# Plot 1: BCV plot
# Shows dispersion vs mean expression
# The tagwise (individual) dispersions should scatter around the trend line
plotBCV(y, main = "Biological Coefficient of Variation\n(tagwise should follow trend)")
 
# Plot 2: MA plot
# x-axis = average log CPM (expression level)
# y-axis = log fold change
# Red points = significant
# Flat cloud = good normalisation
# Tilted cloud = normalisation problem
plotMD(et,
       main = "MA plot\n(red = significant, flat cloud = good)",
       status = ifelse(results[rownames(et$table), "significant"],
                       "Significant", "NotSig"),
       col   = c("Significant" = "#E24B4A", "NotSig" = "#CCCCCC"),
       hl.pch = 16, hl.cex = 0.4)
abline(h = c(-0.5, 0.5), col = "gray", lty = 2)
 
# Plot 3: Top miRNA expression boxplots
# Visual sanity check: do the top miRNAs actually look different?
top_mirna <- results$mirna[1]   # most significant miRNA
top_counts <- log2(counts[top_mirna, ] + 1)
boxplot(top_counts ~ group,
        col    = c("#378ADD", "#E24B4A"),
        main   = paste("Top miRNA:", top_mirna,
                       "\nlog2FC =", round(results$log2FC[1], 2)),
        ylab   = "log2 expression",
        xlab   = "Group",
        notch  = TRUE)   # notch = visual CI for median comparison
 
dev.off()
cat("Plots saved: edgeR_plots.pdf\n")
 
# ── STEP 12: THREE-WAY COMPARISON TABLE ──────────────────────
# Read back the Python results and find miRNAs significant in ALL methods
# This is your most trustworthy biomarker list
 
welch_path  <- file.path(OUTPUT_FOLDER, "limma_results.csv")
deseq_path  <- file.path(OUTPUT_FOLDER, "pydeseq2_results.csv")
 
if (file.exists(welch_path) && file.exists(deseq_path)) {
  cat("\nBuilding three-way comparison...\n")
 
  welch  <- read.csv(welch_path)
  deseq  <- read.csv(deseq_path)
 
  welch_sig  <- welch$mirna[welch$significant  == TRUE]
  deseq_sig  <- deseq$mirna[deseq$significant  == TRUE]
  edger_sig  <- results$mirna[results$significant == TRUE]
 
  # Intersection of all three
  all_three  <- Reduce(intersect, list(welch_sig, deseq_sig, edger_sig))
 
  cat("  Significant in Welch only:          ", length(welch_sig), "\n")
  cat("  Significant in PyDESeq2 only:       ", length(deseq_sig), "\n")
  cat("  Significant in edgeR only:          ", length(edger_sig), "\n")
  cat("  Significant in ALL THREE methods:   ", length(all_three), "\n")
  cat("\n  These", length(all_three), "miRNAs are your most reliable biomarkers:\n")
  print(all_three)
 
  # Save the consensus list
  consensus <- results[results$mirna %in% all_three,
                       c("mirna", "log2FC", "p_value", "fdr")]
  consensus <- consensus[order(abs(consensus$log2FC), decreasing = TRUE), ]
  write.csv(consensus,
            file.path(OUTPUT_FOLDER, "consensus_biomarkers.csv"),
            row.names = FALSE)
  cat("\nConsensus biomarker list saved: consensus_biomarkers.csv\n")
}
 
cat("\n==============================================\n")
cat("DONE\n")
cat("==============================================\n")
cat("Files in", OUTPUT_FOLDER, ":\n")
cat("  edgeR_results.csv       — full edgeR results\n")
cat("  edgeR_plots.pdf         — BCV + MA + boxplot\n")
cat("  consensus_biomarkers.csv — significant in ALL 3 methods\n")

