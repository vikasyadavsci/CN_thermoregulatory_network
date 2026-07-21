library(DEP2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(pheatmap)

# Define the path to your CSV file
file_path <- "Path/to/normalized-WT-cna1-yak1-renamed-dup-split.csv"

# Read the CSV file into a data frame
data <- read.csv(file_path, header = TRUE)

# Format the modification information and generate modified-peptides identifier
unique_pho <- make_unique_ptm(data, gene_name = "Genes", 
                              protein_ID = "Accession", aa = "PTM_AA",
                              pos = "PTM_Loc")

DT::datatable(unique_pho[1:7, c("name", "ID", "Genes", "Accession", "PTM_AA", "PTM_Loc")],
              options = list(scrollX = TRUE, pageLength = 7))

# Take 'Intensity' columns
ecols <- grep("Intensity.", colnames(unique_pho))

# Construct a SE object
se_ptm <- make_se_parse(unique_pho, columns = ecols, 
                        mode = "delim", sep = "_", remove_prefix = TRUE, log2transform = TRUE)

# Perform differential analysis
diff_ptm <- test_diff(se_ptm, type = "control", control = "H25", fdr.type = "BH")
diff_ptmH <- test_diff(se_ptm, type = "manual", test = "H37_vs_H25", fdr.type = "BH")
diff_ptmA <- test_diff(se_ptm, type = "manual", test = "A37_vs_A25", fdr.type = "BH")
diff_ptmY <- test_diff(se_ptm, type = "manual", test = "Y37_vs_Y25", fdr.type = "BH")

# Add rejections
dep_ptm <- DEP2::add_rejections(diff_ptm, alpha = 0.05, lfc = 0.585)
dep_ptmH <- DEP2::add_rejections(diff_ptmH, alpha = 0.05, lfc = 0.585)
dep_ptmA <- DEP2::add_rejections(diff_ptmA, alpha = 0.05, lfc = 0.585)
dep_ptmY <- DEP2::add_rejections(diff_ptmY, alpha = 0.05, lfc = 0.585)

# Extract row data
row_data <- as.data.frame(rowData(dep_ptm))
row_dataH <- as.data.frame(rowData(dep_ptmH))
row_dataA <- as.data.frame(rowData(dep_ptmA))
row_dataY <- as.data.frame(rowData(dep_ptmY))

# Check the structure of row_data to find the correct columns
str(row_data)

library(VennDiagram)
library(grid)

# Extract upregulated significant results for each group
H37p_upregulated <- row_dataH[row_dataH$H37_vs_H25_significant == TRUE & row_dataH$H37_vs_H25_diff > 0.585, "name"]
A37p_upregulated <- row_dataA[row_dataA$A37_vs_A25_significant == TRUE & row_dataA$A37_vs_A25_diff > 0.585, "name"]
H37p_downregulated <- row_dataH[row_dataH$H37_vs_H25_significant == TRUE & row_dataH$H37_vs_H25_diff < -0.585, "name"]
Y37p_downregulated <- row_dataY[row_dataY$Y37_vs_Y25_significant == TRUE & row_dataY$Y37_vs_Y25_diff < -0.585, "name"]

# Create sets for the Venn diagram
H37up_set <- unique(H37p_upregulated)
A37up_set <- unique(A37p_upregulated)
H37down_set <- unique(H37p_downregulated)
Y37down_set <- unique(Y37p_downregulated)

# Check the length of each set to ensure they are populated
cat("Group 1 set length:", length(H37up_set), "\n")
cat("Group 2 set length:", length(A37up_set), "\n")
cat("Group 3 set length:", length(H37down_set), "\n")
cat("Group 4 set length:", length(Y37down_set), "\n")

# Peptides upregulated in A37 but not in H37
A37_specific_up <- setdiff(A37up_set, H37up_set)

# Define the output directory
output_dir <- "CNA1-specific analysis"

# Create the directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}


#Plot heatmap for only A37_up-peptides

A37_up <- row_data %>%
  filter(
    A37_vs_H25_significant & A37_vs_H25_diff > 0.585,
    !A25_vs_H25_significant,
    !H37_vs_H25_significant
  )

# Load gene list
gene_list <- readLines("Turbo-ID-analysis/CT37_CTF37_specific-genes-list.txt")

# Filter for genes of interest
filtered_A37_up <- A37_up %>%
  filter(Genes %in% gene_list)

# Get rownames of selected peptides
selected_names_A37_up <- rownames(filtered_A37_up)

# Extract expression matrix from SE object
expr_mat_A37_up <- assay(se_ptm)[selected_names_A37_up, ]

# Optional: scale rows
expr_scaled_A37_up <- t(scale(t(expr_mat_A37_up)))

# Set rownames for display
rownames(expr_scaled_A37_up) <- make.unique(filtered_A37_up$Genes)

# Group peptides by gene name for heatmap ordering
gene_names <- filtered_A37_up$Genes
order_index <- order(gene_names)  # Group by gene
expr_grouped <- expr_scaled_A37_up[order_index, ]
genes_ordered <- gene_names[order_index]

# Set row labels: display gene name only once
label_names <- genes_ordered
label_names[duplicated(label_names)] <- ""
gene_counts <- table(genes_ordered)
label_names_with_counts <- paste0(label_names, " (", gene_counts[label_names], ")")
label_names_with_counts[duplicated(label_names)] <- ""

# Create internal rownames for uniqueness
rownames(expr_grouped) <- make.unique(genes_ordered)

# Remove Y25 and Y37 columns from the expression matrix
samples_to_keep <- colnames(expr_grouped)[!grepl("Y25|Y37", colnames(expr_grouped))]
expr_grouped_filtered <- expr_grouped[, samples_to_keep]

# ─────────────────────────────────────────────
# 🔹 Motif Annotation
# ─────────────────────────────────────────────

# Combine motif lists with labels
motif1_genes <- motif1_genes %>% mutate(Motif = "PxIxIT")
motif2_genes <- motif2_genes %>% mutate(Motif = "LxVP")
motif_combined <- bind_rows(motif1_genes, motif2_genes)

# Handle duplicates (genes with both motifs)
# Define desired motif order
desired_order <- c("PxIxIT", "LxVP")
motif_status <- motif_combined %>%
  group_by(Gene) %>%
  summarise(Motif = paste(intersect(desired_order, unique(Motif)), collapse = " & "))

# Match to genes in heatmap
heatmap_genes <- data.frame(Gene = genes_ordered)
heatmap_genes$Motif <- motif_status$Motif[match(heatmap_genes$Gene, motif_status$Gene)]
heatmap_genes$Motif[is.na(heatmap_genes$Motif)] <- "None"
rownames(heatmap_genes) <- rownames(expr_grouped)

# Ensure motif column is a factor with desired order
heatmap_genes$Motif <- factor(heatmap_genes$Motif,
                              levels = c("PxIxIT & LxVP", "PxIxIT", "LxVP", "None"))

# ─────────────────────────────────────────────
# 🔹 Peptide Count Annotation
# ─────────────────────────────────────────────

peptide_counts <- gene_counts[genes_ordered]
heatmap_genes$Peptide_Count <- as.numeric(peptide_counts)

# ─────────────────────────────────────────────
# 🔹 Row Annotation Data Frame
# ─────────────────────────────────────────────

row_annotation <- heatmap_genes[, c("Motif", "Peptide_Count")]

# Define annotation colors
motif_colors <- c("PxIxIT & LxVP" = "#009E73",
                  "PxIxIT" = "#E69F00",
                  "LxVP" = "#56B4E9",
                  "None" = "gray80")

annotation_colors <- list(
  Motif = motif_colors,
  Peptide_Count = colorRampPalette(c("white", "darkred"))(100)
)

# ─────────────────────────────────────────────
# 🔹 Cluster by Motif Type
# ─────────────────────────────────────────────

motif_order <- order(row_annotation$Motif)
expr_clustered <- expr_grouped_filtered[motif_order, ]
label_clustered <- label_names_with_counts[motif_order]
row_annotation_clustered <- row_annotation[motif_order, , drop = FALSE]

# ─────────────────────────────────────────────
# 🔹 Plot Heatmap
# ─────────────────────────────────────────────

pheatmap(expr_clustered,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         labels_row = label_clustered,
         annotation_row = row_annotation_clustered,
         annotation_colors = annotation_colors,
         fontsize_col = 10,
         main = "A37-Upregulated Peptides Clustered by Motif Presence and Peptide Count")

# Save to PDF
pdf("Turbo-ID-analysis//A37-Upregulated Peptides TurboID shared Clustered by Motif Presence.pdf", width = 10, height = 8)

pheatmap(expr_clustered,
         cluster_rows = FALSE,  # Already ordered by motif
         cluster_cols = TRUE,
         show_rownames = TRUE,
         labels_row = label_clustered,
         annotation_row = row_annotation_clustered,
         annotation_colors = annotation_colors,
         fontsize_col = 10,
         main = "A37-Upregulated Peptides Clustered by Motif Presence")

dev.off()

write.csv(filtered_A37_up,
          file = "Turbo-ID-analysis/Proteins_shared_phosphoproteome_A37_up-only.csv",
          row.names = FALSE)


# Extract sample names (e.g., "A37", "H25") from column names
sample_names <- gsub("_\\d+$", "", colnames(expr_grouped_filtered))

# Create a mapping of columns to sample groups
sample_groups <- split(colnames(expr_grouped_filtered), sample_names)

# Function to average replicates for each sample
average_replicates <- function(mat, groups) {
  sapply(groups, function(cols) {
    rowMeans(mat[, cols, drop = FALSE])
  })
}

# Apply averaging
expr_averaged <- average_replicates(expr_grouped_filtered, sample_groups)

# Convert to matrix
expr_averaged <- as.matrix(expr_averaged)

expr_clustered_avg <- expr_averaged[motif_order, ]
label_clustered_avg <- label_names_with_counts[motif_order]
row_annotation_clustered <- row_annotation[motif_order, , drop = FALSE]

pheatmap(expr_clustered_avg,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         labels_row = label_clustered_avg,
         annotation_row = row_annotation_clustered,
         annotation_colors = annotation_colors,
         fontsize_col = 10,
         main = "Averaged A37-Upregulated Peptides by Sample")

# Save to PDF
pdf("Turbo-ID-analysis//A37-Upregulated Peptides TurboID shared Clustered by Motif Presence - averaged-heatmap.pdf", width = 10, height = 8)

pheatmap(expr_clustered_avg,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         labels_row = label_clustered_avg,
         annotation_row = row_annotation_clustered,
         annotation_colors = annotation_colors,
         fontsize_col = 10,
         main = "A37-Upregulated Peptides Clustered by Motif Presence averaged")

dev.off()
