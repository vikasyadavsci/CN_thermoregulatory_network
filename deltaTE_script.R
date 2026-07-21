ribo_counts <- read.delim("E:/Ribo-seq-analysis/FD_genome_analysis/Ribo-seq/gene_counts_matrix.tsv", row.names = 1)
rna_counts <- read.delim("E:/Ribo-seq-analysis/FD_genome_analysis/RNA-seq/gene_counts_matrix.tsv", row.names = 1)
sample_info <- read.delim("E:/Ribo-seq-analysis/FD_genome_analysis/Sample_metadata_info.txt", row.names = 1)

output_dir <- "TE_analysis_results"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Save as tab-delimited files for deltaTE
write.table(ribo_counts, "ribo_counts.txt", sep = "\t", quote = FALSE)
write.table(rna_counts,  "rna_counts.txt",  sep = "\t", quote = FALSE)
write.table(sample_info, "sample_info.txt", sep = "\t", quote = FALSE)

# Load required libraries
library(DESeq2)

# Load count matrices and sample information file
ribo_counts <- read.delim("ribo_counts.txt")
rna_counts <- read.delim("rna_counts.txt")
sample_info <- read.delim("sample_info.txt")

sample_info$Condition <- relevel(factor(sample_info$Condition), ref = "H25")

# Create DESeq2 object with interaction term
ddsMat <- DESeqDataSetFromMatrix(
  countData = cbind(ribo_counts, rna_counts),
  colData = sample_info,
  design = ~ Condition + SeqType + Condition:SeqType
)

# Run DESeq2
ddsMat <- DESeq(ddsMat)
resultsNames(ddsMat) # Check available contrasts

rna_info <- subset(sample_info, SeqType == "RNA")
ribo_info <- subset(sample_info, SeqType == "RIBO")

dds_rna <- DESeqDataSetFromMatrix(countData = rna_counts, colData = rna_info, design = ~ Condition)
dds_rna <- DESeq(dds_rna)
resultsNames(dds_rna)
# Extract results
res_C25_vs_H25 <- results(dds_rna, name = "Condition_C25_vs_H25")
res_C25_vs_H25$Comparison <- "C25_vs_H25"

res_C37_vs_H25 <- results(dds_rna, name = "Condition_C37_vs_H25")
res_C37_vs_H25$Comparison <- "C37_vs_H25"

res_H37_vs_H25 <- results(dds_rna, name = "Condition_H37_vs_H25")
res_H37_vs_H25$Comparison <- "H37_vs_H25"

# Save to CSV
write.csv(as.data.frame(res_C25_vs_H25), "C25_vs_H25-RNA.csv")
write.csv(as.data.frame(res_C37_vs_H25), "C37_vs_H25-RNA.csv")
write.csv(as.data.frame(res_H37_vs_H25), "H37_vs_H25-RNA.csv")
# Combine all into one data frame
merged_res <- rbind(res_C25_vs_H25, res_C37_vs_H25, res_H37_vs_H25)
# Save to CSV
write.csv(merged_res, "merged_comparisons-RNA.csv", row.names = TRUE)

dds_ribo <- DESeqDataSetFromMatrix(countData = ribo_counts, colData = ribo_info, design = ~ Condition)
dds_ribo <- DESeq(dds_ribo)
resultsNames(dds_ribo)
# Extract results
res_C25_vs_H25 <- results(dds_ribo, name = "Condition_C25_vs_H25")
res_C25_vs_H25$Comparison <- "C25_vs_H25"

res_C37_vs_H25 <- results(dds_ribo, name = "Condition_C37_vs_H25")
res_C37_vs_H25$Comparison <- "C37_vs_H25"

res_H37_vs_H25 <- results(dds_ribo, name = "Condition_H37_vs_H25")
res_H37_vs_H25$Comparison <- "H37_vs_H25"

# Save to CSV
write.csv(as.data.frame(res_C25_vs_H25), "C25_vs_H25-ribo.csv")
write.csv(as.data.frame(res_C37_vs_H25), "C37_vs_H25-ribo.csv")
write.csv(as.data.frame(res_H37_vs_H25), "H37_vs_H25-ribo.csv")
# Combine all into one data frame
merged_res <- rbind(res_C25_vs_H25, res_C37_vs_H25, res_H37_vs_H25)
# Save to CSV
write.csv(merged_res, "merged_comparisons-ribo.csv", row.names = TRUE)


analyze_condition <- function(cond, ddsMat, dds_rna, dds_ribo, padj_cutoff = 0.05) {
  # TE results (interaction term)
  te_name <- paste0("Condition", cond, ".SeqTypeRNA")
  res_te <- results(ddsMat, name = te_name)
  
  # RNA and Ribo results
  res_rna <- results(dds_rna, contrast = c("Condition", cond, "H25"))
  res_ribo <- results(dds_ribo, contrast = c("Condition", cond, "H25"))
  
  # Optional: shrink LFCs
  res_rna <- lfcShrink(dds_rna, coef = paste0("Condition_", cond, "_vs_H25"), res = res_rna)
  res_ribo <- lfcShrink(dds_ribo, coef = paste0("Condition_", cond, "_vs_H25"), res = res_ribo)
  
  # Classification
  forwarded <- rownames(res_te)[which(res_te$padj > padj_cutoff & res_ribo$padj < padj_cutoff & res_rna$padj < padj_cutoff)]
  exclusive <- rownames(res_te)[which(res_te$padj < padj_cutoff & res_ribo$padj < padj_cutoff & res_rna$padj > padj_cutoff)]
  both <- rownames(res_te)[which(res_te$padj < padj_cutoff & res_ribo$padj < padj_cutoff & res_rna$padj < padj_cutoff)]
  intensified <- rownames(res_te)[which(res_te$padj < padj_cutoff & res_ribo$padj < padj_cutoff & res_rna$padj < padj_cutoff &
                                          res_te$log2FoldChange * res_rna$log2FoldChange > 0)]
  buffered <- rownames(res_te)[which(res_te$padj < padj_cutoff & res_ribo$padj < padj_cutoff & res_rna$padj < padj_cutoff &
                                       res_te$log2FoldChange * res_rna$log2FoldChange < 0)]
  
  # Helper to build full file paths
  out <- function(filename) file.path(output_dir, filename)
  
  # Save results
  write.table(res_te[which(res_te$padj < padj_cutoff), ], out(paste0("DTEGs_", cond, ".txt")), quote = FALSE)
  write.table(res_rna[which(res_rna$padj < padj_cutoff), ], out(paste0("DTGs_", cond, ".txt")), quote = FALSE)
  
  write.table(forwarded, out(paste0("forwarded_genes_", cond, ".txt")), quote = FALSE, row.names = FALSE)
  write.table(exclusive, out(paste0("exclusive_genes_", cond, ".txt")), quote = FALSE, row.names = FALSE)
  write.table(both, out(paste0("both_regulation_genes_", cond, ".txt")), quote = FALSE, row.names = FALSE)
  write.table(intensified, out(paste0("intensified_genes_", cond, ".txt")), quote = FALSE, row.names = FALSE)
  write.table(buffered, out(paste0("buffered_genes_", cond, ".txt")), quote = FALSE, row.names = FALSE)
  
    return(list(
    DTEGs = res_te[which(res_te$padj < padj_cutoff), ],
    DTGs = res_rna[which(res_rna$padj < padj_cutoff), ],
    full_rna = res_rna,
    full_ribo = res_ribo,
    forwarded = forwarded,
    exclusive = exclusive,
    both = both,
    intensified = intensified,
    buffered = buffered
  ))
}


conditions <- c("C25", "C37", "H37")
results_list <- lapply(conditions, function(cond) {
  analyze_condition(cond, ddsMat, dds_rna, dds_ribo)
})
names(results_list) <- conditions



plot_regulation_scatter <- function(res_rna, res_ribo, forwarded, exclusive, intensified, buffered, condition, output_dir = ".") {
  # Extract log2 fold changes
  lfc_rna <- res_rna$log2FoldChange
  lfc_ribo <- res_ribo$log2FoldChange
  
  # Match gene names and remove NAs
  common_genes <- intersect(names(lfc_rna), names(lfc_ribo))
  lfc_rna <- lfc_rna[common_genes]
  lfc_ribo <- lfc_ribo[common_genes]
  
  valid_idx <- which(!is.na(lfc_rna) & !is.na(lfc_ribo))
  lfc_rna <- lfc_rna[valid_idx]
  lfc_ribo <- lfc_ribo[valid_idx]
  valid_genes <- names(lfc_rna)
  
  # Filter category lists to valid genes
  forwarded <- intersect(forwarded, valid_genes)
  exclusive <- intersect(exclusive, valid_genes)
  intensified <- intersect(intensified, valid_genes)
  buffered <- intersect(buffered, valid_genes)
  
  # Diagnostic
  cat("Saving", length(valid_genes), "genes for", condition, "to PDF\n")
  
  # Set up SVG output
  file_path <- file.path(output_dir, paste0("scatter_plot_", condition, "_vs_H25.pdf"))
  pdf(file = file_path, width = 7, height = 7)
  
  # Plot
  max_val <- max(abs(c(lfc_rna, lfc_ribo)), na.rm = TRUE)
  plot(x = lfc_rna, y = lfc_ribo,
       xlab = "RNA-seq log2 fold change",
       ylab = "Ribo-seq log2 fold change",
       main = paste("Translational Regulation:", condition, "vs H25"),
       asp = 1, pch = 16,
       col = rgb(0.5, 0.5, 0.5, 0.1),
       xlim = c(-max_val, max_val), ylim = c(-max_val, max_val), cex = 0.4)
  
  abline(h = 0, v = 0, col = "gray")
  abline(a = 0, b = 1, col = "gray", lty = 2)
  
  # Overlay categories
  points(x = lfc_rna[forwarded], y = lfc_ribo[forwarded], col = "blue", pch = 16)
  points(x = lfc_rna[exclusive], y = lfc_ribo[exclusive], col = "red", pch = 16)
  points(x = lfc_rna[intensified], y = lfc_ribo[intensified], col = "darkgreen", pch = 16)
  points(x = lfc_rna[buffered], y = lfc_ribo[buffered], col = "purple", pch = 16)
  
  legend("topleft", legend = c("Forwarded", "Exclusive", "Intensified", "Buffered"),
         col = c("blue", "red", "darkgreen", "purple"), pch = 16, cex = 0.8)
  
  dev.off()
}

for (cond in names(results_list)) {
  res <- results_list[[cond]]
  plot_regulation_scatter(
    res_rna = res$full_rna,
    res_ribo = res$full_ribo,
    forwarded = res$forwarded,
    exclusive = res$exclusive,
    intensified = res$intensified,
    buffered = res$buffered,
    condition = cond,
    output_dir = output_dir
  )
}

# DTEGs
dtegs_C25 <- rownames(results_list[["C25"]]$DTEGs)
dtegs_C37 <- rownames(results_list[["C37"]]$DTEGs)
dtegs_H37 <- rownames(results_list[["H37"]]$DTEGs)

# DTGs
dtgs_C25 <- rownames(results_list[["C25"]]$DTGs)
dtgs_C37 <- rownames(results_list[["C37"]]$DTGs)
dtgs_H37 <- rownames(results_list[["H37"]]$DTGs)

# Unique DTEGs
unique_dtegs_C37 <- setdiff(dtegs_C37, union(dtegs_C25, dtegs_H37))
write.table(unique_dtegs_C37, file.path(output_dir, "unique_DTEGs_C37.txt"), quote = FALSE, row.names = FALSE)

# Unique DTGs
unique_dtgs_C37 <- setdiff(dtgs_C37, union(dtgs_C25, dtgs_H37))
write.table(unique_dtgs_C37, file.path(output_dir, "unique_DTG_C37.txt"), quote = FALSE, row.names = FALSE)

# Overlap between unique DTEGs and DTGs for C37
overlap_unique_C37 <- intersect(unique_dtegs_C37, unique_dtgs_C37)

# Unique to DTEGs only
only_dtegs_C37 <- setdiff(unique_dtegs_C37, overlap_unique_C37)

# Unique to DTGs only
only_dtgs_C37 <- setdiff(unique_dtgs_C37, overlap_unique_C37)

# Save all three sets
write.table(overlap_unique_C37,
            file.path(output_dir, "shared_unique_DTEG_DTG_C37.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(only_dtegs_C37,
            file.path(output_dir, "unique_only_DTEGs_C37.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(only_dtgs_C37,
            file.path(output_dir, "unique_only_DTG_C37.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# Get full result tables for C37
res_te_C37 <- results_list[["C37"]]$DTEGs
res_rna_C37 <- results_list[["C37"]]$DTGs

# Subset to unique DTEGs and DTGs
res_te_C37_unique <- res_te_C37[unique_dtegs_C37, ]
res_rna_C37_unique <- res_rna_C37[unique_dtgs_C37, ]

# DTEGs
up_dtegs_C37 <- rownames(res_te_C37_unique[res_te_C37_unique$log2FoldChange > 0, ])
down_dtegs_C37 <- rownames(res_te_C37_unique[res_te_C37_unique$log2FoldChange < 0, ])

# DTGs
up_dtgs_C37 <- rownames(res_rna_C37_unique[res_rna_C37_unique$log2FoldChange > 0, ])
down_dtgs_C37 <- rownames(res_rna_C37_unique[res_rna_C37_unique$log2FoldChange < 0, ])

write.table(up_dtegs_C37, file.path(output_dir, "unique_up_DTEGs_C37.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(down_dtegs_C37, file.path(output_dir, "unique_down_DTEGs_C37.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(up_dtgs_C37, file.path(output_dir, "unique_up_DTG_C37.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(down_dtgs_C37, file.path(output_dir, "unique_down_DTG_C37.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

shared_up_C37 <- intersect(up_dtegs_C37, up_dtgs_C37)
shared_down_C37 <- intersect(down_dtegs_C37, down_dtgs_C37)

write.table(shared_up_C37, file.path(output_dir, "shared_up_DTEG_DTG_C37.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(shared_down_C37, file.path(output_dir, "shared_down_DTEG_DTG_C37.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# Exclusive Upregulated
exclusive_up_dtegs_C37 <- setdiff(up_dtegs_C37, up_dtgs_C37)
exclusive_up_dtgs_C37 <- setdiff(up_dtgs_C37, up_dtegs_C37)

# Exclusive Downregulated
exclusive_down_dtegs_C37 <- setdiff(down_dtegs_C37, down_dtgs_C37)
exclusive_down_dtgs_C37 <- setdiff(down_dtgs_C37, down_dtegs_C37)

# Save to files
write.table(exclusive_up_dtegs_C37,
            file.path(output_dir, "exclusive_up_DTEGs_C37.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(exclusive_up_dtgs_C37,
            file.path(output_dir, "exclusive_up_DTG_C37.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(exclusive_down_dtegs_C37,
            file.path(output_dir, "exclusive_down_DTEGs_C37.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(exclusive_down_dtgs_C37,
            file.path(output_dir, "exclusive_down_DTG_C37.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

library(ggplot2)

# Get full DESeq2 results for C37
res_rna <- results_list[["C37"]]$full_rna
res_ribo <- results_list[["C37"]]$full_ribo

# Get unique DTEGs and DTGs for C37
unique_dtegs_C37 <- setdiff(rownames(results_list[["C37"]]$DTEGs),
                            union(rownames(results_list[["C25"]]$DTEGs),
                                  rownames(results_list[["H37"]]$DTEGs)))

unique_dtgs_C37 <- setdiff(rownames(results_list[["C37"]]$DTGs),
                           union(rownames(results_list[["C25"]]$DTGs),
                                 rownames(results_list[["H37"]]$DTGs)))

# Combine unique genes
unique_genes <- union(unique_dtegs_C37, unique_dtgs_C37)

# Extract log2 fold changes
common_genes <- intersect(rownames(res_rna), rownames(res_ribo))
lfc_rna <- res_rna[common_genes, "log2FoldChange"]
lfc_ribo <- res_ribo[common_genes, "log2FoldChange"]

# Subset to unique genes
lfc_rna <- lfc_rna[unique_genes]
lfc_ribo <- lfc_ribo[unique_genes]

# Remove NAs
valid_idx <- which(!is.na(lfc_rna) & !is.na(lfc_ribo))
lfc_rna <- lfc_rna[valid_idx]
lfc_ribo <- lfc_ribo[valid_idx]
genes <- names(lfc_rna)

# Create data frame
df <- data.frame(
  Gene = genes,
  RNA_log2FC = lfc_rna,
  Ribo_log2FC = lfc_ribo,
  Quadrant = factor(
    ifelse(lfc_rna > 0 & lfc_ribo > 0, "Q1: Up-Up",
    ifelse(lfc_rna < 0 & lfc_ribo > 0, "Q2: Down-Up",
    ifelse(lfc_rna < 0 & lfc_ribo < 0, "Q3: Down-Down",
    ifelse(lfc_rna > 0 & lfc_ribo < 0, "Q4: Up-Down", "Other")))),
    levels = c("Q1: Up-Up", "Q2: Down-Up", "Q3: Down-Down", "Q4: Up-Down", "Other")
  )
)

# Plot
p <- ggplot(df, aes(x = RNA_log2FC, y = Ribo_log2FC, color = Quadrant)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c(
    "Q1: Up-Up" = "forestgreen",
    "Q2: Down-Up" = "purple",
    "Q3: Down-Down" = "firebrick",
    "Q4: Up-Down" = "orange",
    "Other" = "gray70"
  )) +
  labs(
    title = "Quadrant Scatter Plot: C37-Specific Genes",
    x = "RNA-seq log2 Fold Change",
    y = "Ribo-seq log2 Fold Change",
    color = "Regulation Pattern"
  ) +
  theme_minimal(base_size = 14)

# Save plot
ggsave(file.path(output_dir, "quadrant_scatter_C37.pdf"), plot = p, width = 7, height = 6)


library(VennDiagram)
library(grid)

# Create Venn for Upregulated genes
venn_up <- venn.diagram(
  x = list(
    DTEG_Up = up_dtegs_C37,
    DTG_Up = up_dtgs_C37
  ),
  filename = NULL,  # Don't save yet
  fill = c("skyblue", "salmon"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  main = "C37 Upregulated Genes: DTEG vs DTG"
)

pdf(file.path(output_dir, "venn_upregulated_C37.pdf"))
grid.draw(venn_up)
dev.off()

# Create Venn for Downregulated genes
venn_down <- venn.diagram(
  x = list(
    DTEG_Down = down_dtegs_C37,
    DTG_Down = down_dtgs_C37
  ),
  filename = NULL,
  fill = c("lightgreen", "plum"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  main = "C37 Downregulated Genes: DTEG vs DTG"
)

pdf(file.path(output_dir, "venn_downregulated_C37.pdf"))
grid.draw(venn_down)
dev.off()


library(VennDiagram)
library(grid)

# DTEG Venn
venn_dteg <- venn.diagram(
  x = list(C25 = dtegs_C25, C37 = dtegs_C37, H37 = dtegs_H37),
  filename = NULL,
  fill = c("skyblue", "orange", "seagreen3"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.pos = 0,
  cat.dist = 0.05,
  margin = 0.1
)

grid.newpage()
grid.draw(venn_dteg)

venn_dtg <- venn.diagram(
  x = list(C25 = dtgs_C25, C37 = dtgs_C37, H37 = dtgs_H37),
  filename = NULL,
  fill = c("skyblue", "orange", "seagreen3"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.pos = 0,
  cat.dist = 0.05,
  margin = 0.1
)

grid.newpage()
grid.draw(venn_dtg)

# C37 overlaps
shared_C37 <- intersect(dtegs_C37, dtgs_C37)
only_dtegs_C37 <- setdiff(dtegs_C37, dtgs_C37)
only_dtgs_C37 <- setdiff(dtgs_C37, dtegs_C37)

# Save these
write.table(shared_C37, file.path(output_dir, "C37_shared_DTEG_DTG.txt"), quote = FALSE, row.names = FALSE)
write.table(only_dtegs_C37, file.path(output_dir, "C37_only_DTEGs.txt"), quote = FALSE, row.names = FALSE)
write.table(only_dtgs_C37, file.path(output_dir, "C37_only_DTG.txt"), quote = FALSE, row.names = FALSE)

venn_c37 <- venn.diagram(
  x = list(DTEG = dtegs_C37, DTG = dtgs_C37),
  filename = NULL,  # Don't save to file
  fill = c("tomato", "steelblue"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.pos = 0,
  cat.dist = 0.05,
  margin = 0.1
)

grid.newpage()
grid.draw(venn_c37)

venn_unique_c37 <- venn.diagram(
  x = list(DTEG = unique_dtegs_C37, DTG = unique_dtgs_C37),
  filename = NULL,
  fill = c("tomato", "steelblue"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.pos = 0,
  cat.dist = 0.05,
  margin = 0.1
)

grid.newpage()
grid.draw(venn_unique_c37)












# Obtain fold changes for translation efficiency (TE)
res_C25 <- results(ddsMat, name = "Condition_C25_vs_H25")
res_C37 <- results(ddsMat, name = "Condition_C37_vs_H25")
res_H37 <- results(ddsMat, name = "Condition_H37_vs_H25")

# Store list of Differential Translation Efficiency Genes (DTEGs)
write.table(res[which(res$padj < 0.05), ], "DTEGs.txt", quote = FALSE)

# Run DESeq2 for RNA-seq to obtain Differentially Transcribed Genes (DTGs)
ddsMat_rna <- DESeqDataSetFromMatrix(
  countData = rna_counts,
  colData = sample_info[sample_info$SeqType == "RNA", ],
  design = ~ Condition
)
ddsMat_rna <- DESeq(ddsMat_rna)
resultsNames(ddsMat_rna)
res_rna <- results(ddsMat_rna, name = "Condition_C37_vs_C25")
res_rna <- lfcShrink(ddsMat_rna, coef = "Condition_C37_vs_C25", res = res_rna)

# Store list of DTGs
write.table(res_rna[which(res_rna$padj < 0.05), ], "DTGs.txt", quote = FALSE)

# Run DESeq2 for Ribo-seq to obtain RPF count changes
ddsMat_ribo <- DESeqDataSetFromMatrix(
  countData = ribo_counts,
  colData = sample_info[sample_info$SeqType == "RIBO", ],
  design = ~ Condition
)
ddsMat_ribo <- DESeq(ddsMat_ribo)
resultsNames(ddsMat_ribo)
res_ribo <- results(ddsMat_ribo, name = "Condition_C37_vs_C25")
res_ribo <- lfcShrink(ddsMat_ribo, coef = "Condition_C37_vs_C25", res = res_ribo)

# Categorize genes into regulation classes
forwarded <- rownames(res)[which(res$padj > 0.05 & res_ribo$padj < 0.05 & res_rna$padj < 0.05)]
exclusive <- rownames(res)[which(res$padj < 0.05 & res_ribo$padj < 0.05 & res_rna$padj > 0.05)]
both <- rownames(res)[which(res$padj < 0.05 & res_ribo$padj < 0.05 & res_rna$padj < 0.05)]

# Further categorize into intensified and buffered genes
intensified <- rownames(res)[which(res$padj < 0.05 & res_ribo$padj < 0.05 & res_rna$padj < 0.05 & 
                                     res[,2] * res_rna[,2] > 0)]
buffered <- rownames(res)[which(res$padj < 0.05 & res_ribo$padj < 0.05 & res_rna$padj < 0.05 & 
                                  res[,2] * res_rna[,2] < 0)]

# Save categorized gene lists
write.table(forwarded, "forwarded_genes.txt", quote = FALSE, row.names = FALSE)
write.table(exclusive, "exclusive_genes.txt", quote = FALSE, row.names = FALSE)
write.table(both, "both_regulation_genes.txt", quote = FALSE, row.names = FALSE)
write.table(intensified, "intensified_genes.txt", quote = FALSE, row.names = FALSE)
write.table(buffered, "buffered_genes.txt", quote = FALSE, row.names = FALSE)
print("Analysis completed. Gene lists saved.")

# Visualization of global translational and transcriptional regulation
max_val <- max(res_ribo[,2], res_rna[,2], na.rm = TRUE)
plot(y = res_ribo[,2], x = res_rna[,2],
     xlab = "RNA-seq log2 fold change",
     ylab = "Ribo-seq log2 fold change", asp = 1, pch = 16,
     col = rgb(128/255, 128/255, 128/255, 0.1),
     ylim = c(-max_val, max_val), xlim = c(-max_val, max_val), cex = 0.4)

abline(a = 0, b = 1, col = "gray")
abline(h = 0, v = 0, col = "gray")

points(y = res_ribo[forwarded, 2], x = res_rna[forwarded, 2],
       pch = 16, col = rgb(0, 0, 1, 1))

points(y = res_ribo[exclusive, 2], x = res_rna[exclusive, 2],
       pch = 16, col = rgb(1, 0, 0, 1))

points(y = res_ribo[intensified, 2], x = res_rna[intensified, 2],
       pch = 16, col = rgb(1, 0, 1, 1))

points(y = res_ribo[buffered, 2], x = res_rna[buffered, 2],
       pch = 16, col = rgb(1, 0, 1, 1))

print("Analysis completed. Gene lists saved. Visualization generated.")


