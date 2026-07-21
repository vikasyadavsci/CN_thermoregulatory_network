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

dds_ribo <- DESeqDataSetFromMatrix(countData = ribo_counts, colData = ribo_info, design = ~ Condition)
dds_ribo <- DESeq(dds_ribo)

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

library(ggplot2)

# Create output subdirectory
pca_dir <- file.path(output_dir, "PCA_plots")
if (!dir.exists(pca_dir)) dir.create(pca_dir)

# Variance-stabilizing transformation
vsd_rna <- vst(dds_rna, blind = TRUE)

# PCA plot
p_rna <- plotPCA(vsd_rna, intgroup = "Condition") +
  ggtitle("PCA: RNA-seq Samples")

# Save
ggsave(file.path(pca_dir, "PCA_RNAseq_C37.pdf"), plot = p_rna, width = 6, height = 5)


# Variance-stabilizing transformation
vsd_ribo <- vst(dds_ribo, blind = TRUE)

# PCA plot
p_ribo <- plotPCA(vsd_ribo, intgroup = "Condition") +
  ggtitle("PCA: Ribo-seq Samples")

# Save
ggsave(file.path(pca_dir, "PCA_Riboseq_C37.pdf"), plot = p_ribo, width = 6, height = 5)


library(pheatmap)

# Example: load gene list from file
genes_of_interest <- readLines("Mito-genes.txt")

# Extract vst-transformed matrices
vsd_rna_mat <- assay(vsd_rna)
vsd_ribo_mat <- assay(vsd_ribo)

# Subset to genes of interest
vsd_rna_subset <- vsd_rna_mat[rownames(vsd_rna_mat) %in% genes_of_interest, ]
vsd_ribo_subset <- vsd_ribo_mat[rownames(vsd_ribo_mat) %in% genes_of_interest, ]

# Optional: scale rows (z-score)
vsd_rna_scaled <- t(scale(t(vsd_rna_subset)))
vsd_ribo_scaled <- t(scale(t(vsd_ribo_subset)))

# Ensure sample order matches annotation
annotation_df_rna <- sample_info[colnames(vsd_rna_scaled), , drop = FALSE]
annotation_df_ribo <- sample_info[colnames(vsd_ribo_scaled), , drop = FALSE]

# Unified color scale
combined_range <- range(c(vsd_rna_scaled, vsd_ribo_scaled), na.rm = TRUE)
breaks_seq <- seq(combined_range[1], combined_range[2], length.out = 100)

# Create output directory
heatmap_dir <- file.path(output_dir, "heatmaps")
if (!dir.exists(heatmap_dir)) dir.create(heatmap_dir)

# RNA-seq heatmap
pdf(file.path(heatmap_dir, "heatmap_RNAseq_genes_of_interest.pdf"), width = 8, height = 10)
pheatmap(vsd_rna_scaled,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_col = annotation_df_rna,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 8,
         fontsize_col = 10,
         breaks = breaks_seq,
         main = "RNA-seq: Genes of Interest")
dev.off()

# Ribo-seq heatmap
pdf(file.path(heatmap_dir, "heatmap_Riboseq_genes_of_interest.pdf"), width = 8, height = 10)
pheatmap(vsd_ribo_scaled,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_col = annotation_df_ribo,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 8,
         fontsize_col = 10,
         breaks = breaks_seq,
         main = "Ribo-seq: Genes of Interest")
dev.off()

