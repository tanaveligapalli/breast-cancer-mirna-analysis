
# Why voom for this data:
#   - Data is microarray log2 intensities back-transformed to counts
#   - voom models mean-variance relationship and creates precision weights
#   - More statistically honest than raw DESeq2/edgeR on non-true-counts
#   - More powerful than plain limma which ignores count-data properties
#

library(limma)
library(edgeR)

# ── CONFIG ────────────────────────────────────────────────────
DATA_PATH     <- "/Users/tanaveligapalli/Downloads/miRNA_results/final_dataset.csv"
OUTPUT_FOLDER <- "/Users/tanaveligapalli/R/results"

# ── STEP 1: LOAD AND SPLIT DATA ───────────────────────────────
cat("Loading data...\n")
full_data  <- read.csv(DATA_PATH, row.names=1, check.names=FALSE)

# Separate miRNA columns from clinical
mirna_cols <- grep("^hsa-", colnames(full_data), value=TRUE)
expr       <- full_data[, mirna_cols]
metadata   <- full_data[, c("diagnosis", "label")]

# Remove prostate disease
keep       <- metadata$diagnosis != "prostate disease"
expr       <- expr[keep, ]
metadata   <- metadata[keep, ]

# Group factor
group <- factor(
  ifelse(metadata$diagnosis == "breast cancer", "cancer", "healthy"),
  levels = c("healthy", "cancer")
)

cat("  Cancer: ", sum(group=="cancer"),
    "| Healthy:", sum(group=="healthy"),
    "| miRNAs:", ncol(expr), "\n")

# ── STEP 2: REVERSE LOG2 → COUNTS ────────────────────────────
# voom works on counts, same as edgeR
# We reverse log2 to get approximate raw intensities
counts <- round(2^as.matrix(expr))
counts <- pmax(counts, 1)
counts <- t(counts)   # voom expects miRNAs as rows, patients as columns

cat("  Count matrix:", nrow(counts), "miRNAs x", ncol(counts), "patients\n")
cat("  Range:", min(counts), "to", max(counts), "\n")

# ── STEP 3: BUILD DGEList AND NORMALISE ───────────────────────
# Same start as edgeR — DGEList + TMM normalisation
# voom uses the normalisation factors from edgeR internally
cat("\nNormalising with TMM...\n")
y <- DGEList(counts=counts, group=group)
y <- normLibSizes(y, method = "TMM")

cat("  Norm factors — Min:", round(min(y$samples$norm.factors), 4),
    "Max:", round(max(y$samples$norm.factors), 4), "\n")

# ── STEP 4: DESIGN MATRIX ─────────────────────────────────────
# The design matrix tells limma what the experimental groups are
# model.matrix() converts your group factor into a numeric matrix
# that the linear model can use
#
# ~0+group gives you:
#   column 1 (grouphealthy): 1 if healthy, 0 if cancer
#   column 2 (groupcancer):  1 if cancer,  0 if healthy
# This is called a "means parameterisation" — each coefficient
# represents the mean expression in one group directly
#
# Alternative: ~group gives intercept + contrast
# We use ~0+group because it makes the contrast definition clearer

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)  # rename to "healthy" and "cancer"

cat("\nDesign matrix (first 3 rows):\n")
print(head(design, 3))



cat("\nRunning voom transformation...\n")

# Open PDF for all plots
pdf(file.path(OUTPUT_FOLDER, "voom_plots.pdf"), width=14, height=5)
par(mfrow=c(1,3))

# voom() with plot=TRUE draws the mean-variance plot automatically
v <- voom(y, design, plot=TRUE)
title("voom: Mean-Variance Trend\n(decreasing curve = correct)")

cat("  voom complete\n")
cat("  Weight range:", round(min(v$weights),4),
    "to", round(max(v$weights),4), "\n")
cat("  (higher weight = more precise measurement)\n")

# ── STEP 6: FIT LINEAR MODEL ──────────────────────────────────
# lmFit fits a linear model for each miRNA simultaneously
# It uses the precision weights from voom
# Each miRNA gets: expression ~ intercept + group_effect
# The coefficients are the mean log2-CPM for each group
cat("\nFitting linear models...\n")
fit <- lmFit(v, design)

# ── STEP 7: DEFINE CONTRAST ───────────────────────────────────
# A contrast specifies which comparison you want
# cancer - healthy means: fold change = cancer mean - healthy mean
# Positive contrast value = higher in cancer
#
# makeContrasts() takes your contrast as a text formula
# The names must match colnames(design) exactly
contrast_matrix <- makeContrasts(
  cancer_vs_healthy = cancer - healthy,
  levels = design
)

cat("  Contrast: cancer vs healthy\n")

# Apply contrast to the fitted model
fit2 <- contrasts.fit(fit, contrast_matrix)

# ── STEP 8: EMPIRICAL BAYES MODERATION ───────────────────────
# eBayes() applies limma's empirical Bayes moderation
# NOW this is statistically valid because voom's weights have
# already corrected the heteroscedasticity (unequal variances)
#
# What eBayes does:
#   - Estimates a global prior variance across all miRNAs
#   - Shrinks each miRNA's variance toward this prior
#   - Computes moderated t-statistics using shrunk variances
#   - Computes p-values from moderated t-distribution
#
# trend=TRUE: tells eBayes the prior variance may itself depend
# on expression level — appropriate after voom which models
# exactly this mean-variance relationship
# robust=TRUE: makes the prior estimation resistant to outlier miRNAs

cat("Applying empirical Bayes moderation...\n")
fit2 <- eBayes(fit2, trend=TRUE, robust=TRUE)

cat("  Prior df:", round(fit2$df.prior, 2), "\n")
cat("  (should be 3-30 for well-behaved data)\n")

# ── STEP 9: EXTRACT RESULTS ───────────────────────────────────
# topTable() extracts the statistical results
# number=Inf returns all miRNAs, not just top ones
# adjust="BH" applies Benjamini-Hochberg FDR correction
# sort.by="P" sorts by ascending p-value
results        <- topTable(fit2,
                           coef    = "cancer_vs_healthy",
                           number  = Inf,
                           adjust  = "BH",
                           sort.by = "P")
results$mirna  <- rownames(results)
results        <- results[, c("mirna","logFC","AveExpr",
                               "t","P.Value","adj.P.Val","B")]

# Rename to match our convention
colnames(results) <- c("mirna","log2FC","avg_expr",
                        "t_stat","p_value","fdr","B_stat")

# Significance flag
results$significant <- results$fdr < 0.05 & abs(results$log2FC) > 0.5

n_sig <- sum(results$significant)
n_up  <- sum(results$significant & results$log2FC > 0)
n_dn  <- sum(results$significant & results$log2FC < 0)

cat("\n==============================================\n")
cat("RESULTS\n")
cat("==============================================\n")
cat("Significant (FDR<0.05, |logFC|>0.5):", n_sig, "\n")
cat("  Upregulated in cancer:  ", n_up, "\n")
cat("  Downregulated in cancer:", n_dn, "\n")

cat("\nTop 10 upregulated:\n")
top_up <- head(results[results$log2FC > 0 & results$significant, ], 10)
print(top_up[, c("mirna","log2FC","p_value","fdr")], row.names=FALSE)

cat("\nTop 10 downregulated:\n")
top_dn <- head(results[results$log2FC < 0 & results$significant, ], 10)
print(top_dn[, c("mirna","log2FC","p_value","fdr")], row.names=FALSE)

# ── STEP 10: PLOTS ────────────────────────────────────────────
# Volcano plot
sig_col <- ifelse(results$significant,
                  ifelse(results$log2FC > 0, "#E24B4A", "#378ADD"),
                  "#CCCCCC")
plot(results$log2FC,
     -log10(results$p_value + 1e-10),
     col  = sig_col,
     pch  = 16,
     cex  = 0.4,
     xlab = "log2 Fold Change (Cancer vs Healthy)",
     ylab = "-log10(p-value)",
     main = paste("limma-voom Volcano\n(",
                  n_up, "up,", n_dn, "down)"))
abline(v=c(-0.5, 0.5), col="gray", lty=2)
abline(h=-log10(0.05), col="gray", lty=2)

# Label top 5
top5 <- head(results[results$significant, ], 5)
text(top5$log2FC,
     -log10(top5$p_value + 1e-10),
     labels = sub("hsa-","", top5$mirna),
     cex = 0.6, pos = 4)

# MA plot
plotMD(fit2,
       column = 1,
       status = ifelse(results[rownames(fit2$coefficients), "significant"],
                       "Sig", "NS"),
       col    = c("Sig"="#E24B4A", "NS"="#CCCCCC"),
       main   = "MA plot (limma-voom)\n(flat cloud = good normalisation)",
       hl.cex = 0.4)
abline(h=c(-0.5, 0.5), col="gray", lty=2)

dev.off()
cat("Plots saved: voom_plots.pdf\n")

# ── STEP 11: SAVE RESULTS ────────────────────────────────────
write.csv(results,
          file.path(OUTPUT_FOLDER, "voom_results.csv"),
          row.names=FALSE)
cat("Saved: voom_results.csv\n")

# ── STEP 12: THREE-WAY CONSENSUS ─────────────────────────────
# Find miRNAs significant in ALL THREE methods:
# Welch t-test, edgeR, AND limma-voom
welch_path <- file.path(OUTPUT_FOLDER, "limma_results.csv")
edger_path <- file.path(OUTPUT_FOLDER, "edgeR_results.csv")

if (file.exists(welch_path) && file.exists(edger_path)) {
    cat("\nBuilding three-way consensus...\n")

    welch <- read.csv(welch_path)
    edger <- read.csv(edger_path)

    welch_sig <- welch$mirna[welch$significant == TRUE]
    edger_sig <- edger$mirna[edger$significant == TRUE]
    voom_sig  <- results$mirna[results$significant == TRUE]

    all_three <- Reduce(intersect,
                        list(welch_sig, edger_sig, voom_sig))

    cat("  Welch significant:         ", length(welch_sig), "\n")
    cat("  edgeR significant:         ", length(edger_sig), "\n")
    cat("  limma-voom significant:    ", length(voom_sig),  "\n")
    cat("  ALL THREE consensus:       ", length(all_three), "\n")

    cat("\n  Consensus miRNAs:\n")
    consensus <- results[results$mirna %in% all_three,
                         c("mirna","log2FC","p_value","fdr")]
    consensus <- consensus[order(abs(consensus$log2FC),
                                 decreasing=TRUE), ]
    print(consensus, row.names=FALSE)

    write.csv(consensus,
              file.path(OUTPUT_FOLDER, "voom_consensus.csv"),
              row.names=FALSE)
    cat("\nSaved: voom_consensus.csv\n")
}

cat("\n==============================================\n")
cat("DONE\n")
cat("==============================================\n")
cat("voom_results.csv   — full limma-voom results\n")
cat("voom_plots.pdf     — voom + volcano + MA\n")
cat("voom_consensus.csv — significant in ALL 3 methods\n")
