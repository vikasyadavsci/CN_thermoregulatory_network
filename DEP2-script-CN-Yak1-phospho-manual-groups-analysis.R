library(DEP2)
library(dplyr)
library(ggplot2)

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


library(UpSetR)
# Create a combined list of unique elements
all_elements <- unique(c(H37up_set, A37up_set, H37down_set, Y37down_set))

# Create a data frame with binary indicators
data <- data.frame(
  Element = all_elements,
  H37up = ifelse(all_elements %in% H37up_set, 1, 0),
  A37up = ifelse(all_elements %in% A37up_set, 1, 0),
  H37down = ifelse(all_elements %in% H37down_set, 1, 0),
  Y37down = ifelse(all_elements %in% Y37down_set, 1, 0)
)

# Generate the UpSet plot
upset(data, sets = c("H37up", "A37up", "H37down", "Y37down"),
      keep.order = TRUE, 
      sets.bar.color = c("#ff7f00", "#1f78b4", "#ff7f90", "#33a02c"),
      matrix.color = "#1f78b4", 
      main.bar.color = "darkgrey",
      mainbar.y.label = "Intersection Size", 
      sets.x.label = "Set Size")

# Create a named list of sets
venn_list <- list(
  H37_Up = H37up_set,
  A37_Up = A37up_set,
  H37_Down = H37down_set,
  Y37_Down = Y37down_set
)

# Generate the Venn diagram
venn.plot <- venn.diagram(
  x = venn_list,
  filename = NULL,  # NULL returns a grid object instead of saving to file
  fill = c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00"),
  alpha = 0.5,
  cex = 1.2,
  cat.cex = 1.2,
  cat.pos = 0,
  cat.dist = 0.05,
  margin = 0.1,
  main = "Venn Diagram of Differential PTMs"
)

# Draw the diagram
grid.newpage()
grid.draw(venn.plot)


# Find peptides shared between Y37p and A37p exclusively
shared_Y37down_A37up <- intersect(Y37p_downregulated, A37p_upregulated)
exclusive_Y37down_A37up <- setdiff(shared_Y37down_A37up, c(H37p_upregulated, H37p_downregulated))

# Save to a CSV file
write.csv(exclusive_Y37down_A37up, file = "exclusive_peptides.csv", row.names = FALSE)

# Load the combined data file
exclusive_peptides <- read.csv("exclusive_peptides.csv")

# Merge with row_data to extract respective rows for each peptide
merged_data <- merge(exclusive_peptides, row_data, by.x = "x", by.y = "name")

# Save the merged data to a new CSV file
write.csv(merged_data, file = "merged_exclusive_peptides_data.csv", row.names = FALSE)

# Load PTM-exclusive peptides
ptm_exclusive <- read.csv("merged_exclusive_peptides_data.csv")

# Ensure TurboID data is in memory (already created in your script)
# upregulated_CT37_CTF37_exclusive

# Merge on gene symbol
merged_cross <- merge(ptm_exclusive, upregulated_CT37_CTF37_exclusive, 
                      by = "Genes", suffixes = c("_PTM", "_TurboID"))

# Or merge on protein accession
# merged_cross <- merge(ptm_exclusive, upregulated_CT37_CTF37_exclusive, 
#                       by = "Accession", suffixes = c("_PTM", "_TurboID"))

# Save to CSV
write.csv(merged_cross, "PTM_TurboID_overlap.csv", row.names = FALSE)

# Quick look
head(merged_cross)



#Heatmap generation

library(SummarizedExperiment)
library(dplyr)

# Extract assay matrix (log2 intensities)
mat <- assay(se_ptm)

# Extract row metadata
rd <- rowData(se_ptm)

# Keep only exclusive peptides
exclusive_idx <- rd$name %in% exclusive_Y37down_A37up
mat_exclusive <- mat[exclusive_idx, ]

# Add peptide names as rownames
rownames(mat_exclusive) <- rd$name[exclusive_idx]

mat_scaled <- t(scale(t(mat_exclusive)))

library(ComplexHeatmap)
library(circlize)

sample_info <- colData(se_ptm) %>% as.data.frame()

ha <- HeatmapAnnotation(
  Condition = sample_info$condition,
  Replicate = sample_info$replicate,
  col = list(
    Condition = c(
      H25="#1f78b4", H37="#a6cee3",
      A25="#33a02c", A37="#b2df8a",
      Y25="#e31a1c", Y37="#fb9a99"
    ),
    Replicate = c("1"="grey70", "2"="grey40", "3"="black")
  )
)

# Save to PDF
pdf("CN-Yak1/A37Up_Y37Down_peptides_heatmap.pdf", width = 10, height = 8)

Heatmap(
  mat_scaled,
  name = "Z-score",
  top_annotation = ha,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 8),
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  column_title = "Exclusive Y37-down ∩ A37-up Peptides",
  row_title = "Peptides"
)

dev.off()

# ─────────────────────────────────────────────
# 🔹 Motif Annotation
# ─────────────────────────────────────────────

# Map gene names back from unique_pho
gene_map <- unique_pho[, c("name", "Genes")]

genes_exclusive <- gene_map$Genes[match(rownames(mat_scaled), gene_map$name)]



library(readxl)
library(dplyr)

motif_file <- "Motif-analysis-FungiDB.xlsx"

motif1_genes <- read_excel(motif_file, sheet = "PxIxIT_neoformans")
motif2_genes <- read_excel(motif_file, sheet = "LxVP_neoformans")
# Keep only the Gene column
motif1_genes <- motif1_genes %>% select(Gene)
motif2_genes <- motif2_genes %>% select(Gene)

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
heatmap_genes <- data.frame(Gene = genes_exclusive)
heatmap_genes$Motif <- motif_status$Motif[match(heatmap_genes$Gene, motif_status$Gene)]
heatmap_genes$Motif[is.na(heatmap_genes$Motif)] <- "None"

# Ensure motif column is a factor with desired order
heatmap_genes$Motif <- factor(heatmap_genes$Motif,
                              levels = c("PxIxIT & LxVP", "PxIxIT", "LxVP", "None"))

# ─────────────────────────────────────────────
# 🔹 Peptide Count Annotation
# ─────────────────────────────────────────────

gene_counts <- table(genes_exclusive)
heatmap_genes$Peptide_Count <- as.numeric(gene_counts[genes_exclusive])

# ─────────────────────────────────────────────
# 🔹 Row Annotation Data Frame
# ─────────────────────────────────────────────

# Order rows by Motif (and optionally by gene within motif)
row_order <- order(heatmap_genes$Motif, heatmap_genes$Gene)

mat_ordered <- mat_scaled[row_order, ]
heatmap_genes_ordered <- heatmap_genes[row_order, ]


library(ComplexHeatmap)
library(circlize)

row_annot <- HeatmapAnnotation(
  df = heatmap_genes_ordered[, c("Motif", "Peptide_Count")],
  col = list(
    Motif = c(
      "PxIxIT & LxVP" = "#009E73",
      "PxIxIT" = "#E69F00",
      "LxVP" = "#56B4E9",
      "None" = "gray80"
    ),
    Peptide_Count = colorRamp2(
      c(min(heatmap_genes_ordered$Peptide_Count),
        max(heatmap_genes_ordered$Peptide_Count)),
      c("white", "darkred")
    )
  ),
  which = "row"
)


sample_info <- colData(se_ptm) %>% as.data.frame()

ha <- HeatmapAnnotation(
  Condition = sample_info$condition,
  Replicate = sample_info$replicate,
  col = list(
    Condition = c(
      H25="#1f78b4", H37="#a6cee3",
      A25="#33a02c", A37="#b2df8a",
      Y25="#e31a1c", Y37="#fb9a99"
    ),
    Replicate = c("1"="grey70", "2"="grey40", "3"="black")
  )
)

# ─────────────────────────────────────────────
# 🔹 Plot Heatmap
# ─────────────────────────────────────────────

pdf("CN-Yak1/A37Up_Y37Down_peptides_heatmap_clustered_by_motif.pdf", width = 10, height = 10)

Heatmap(
  mat_ordered,
  name = "Z-score",
  top_annotation = ha,
  left_annotation = row_annot,
  cluster_rows = FALSE,   # we already ordered by motif
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 8),
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  column_title = "Exclusive Y37-down ∩ A37-up Peptides",
  row_title = "Peptides"
)

dev.off()
