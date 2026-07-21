library(DEP2)
library(dplyr)
library(ggplot2)

# Define the path to your CSV file
file_path <- "Path/to/normalized-data-TurboID-renamed.csv"

# Read the CSV file into a data frame
data <- read.csv(file_path, header = TRUE)

# Format name(gene symbol) and id(protein ID). 
# Generate a unique names for each protein. names and ids are columns in table
unique_pg <- make_unique(data, names = "Genes", ids = "Accession", delim = ";")

# Take expression columns(LFQ intensity in this cases).
ecols <- grep("LFQ.Intensity.", colnames(unique_pg))

# Construct SE. The experiement design is exctracted from column.
se_pg <- DEP2::make_se_parse(unique_pg, columns = ecols, mode = "delim", 
                             sep = "_", remove_prefix = T, log2transform = T)

## Test on manual contrasts
diff_control <- test_diff(se_pg, type = "control", control = "TID25", fdr.type = "BH")

## Add significant rejections for features, based on
dep_control <- add_rejections(diff_control, alpha = 0.05, lfc = 1)

# Extract row data
row_data_control <- as.data.frame(rowData(dep_control))

# Check the structure of row_data to find the correct columns
str(row_data_control)

# Extract normalized intensity values from SummarizedExperiment object
pca_data <- t(assay(se_pg))

# Perform PCA
pca_result <- prcomp(pca_data, scale. = TRUE)

# Extract sample metadata from DEP2 object
metadata <- colData(se_pg)[, c("ID", "condition", "replicate")]

# Convert to a dataframe for merging
metadata_df <- as.data.frame(metadata)
colnames(metadata_df) <- c("ID", "Condition", "replicate")

# Create PCA dataframe and merge with metadata
pca_df <- data.frame(ID = colnames(se_pg),
                     PC1 = pca_result$x[, 1],
                     PC2 = pca_result$x[, 2])

pca_df <- merge(pca_df, metadata_df, by = "ID")

# Compute cluster boundary points using convex hull
hull_points <- pca_df %>%
  group_by(Condition) %>%
  slice(chull(PC1, PC2))

# Compute label positions by selecting the **furthest point** along PC1
label_positions <- hull_points %>%
  group_by(Condition) %>%
  filter(PC1 == max(PC1))  # Position labels at the rightmost point

# Generate PCA plot with convex hulls and shifted cluster labels
ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 3) +
  geom_polygon(data = hull_points, aes(fill = Condition), alpha = 0.2) +
  geom_text(data = label_positions, aes(label = Condition), hjust = -0.2, size = 4) +  # Shift labels outside
  labs(title = "PCA of PTM Intensities with Cluster Labels Next to Groups",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme_minimal()

# Generate heatmap of significant PTMs
plot_heatmap(dep_control)

# Load required libraries
library(pheatmap)
library(RColorBrewer)

# Extract significant proteins
sig_proteins <- get_signicant(dep_control)

# Extract row names from the S4 object
sig_names <- rownames(sig_proteins)

# Now subset the assay matrix
expr_matrix <- assay(dep_control)[sig_names, ]

# Optional: scale rows (z-score normalization)
scaled_matrix <- t(scale(t(expr_matrix)))

# Set annotation for columns (samples)
annotation_col <- data.frame(Condition = metadata_df$Condition)
rownames(annotation_col) <- metadata_df$ID

# Choose a color palette
heat_colors <- colorRampPalette(c("blue", "#F7F7F7", "red"))(100)

# Plot heatmap
pheatmap(scaled_matrix,
         color = heat_colors,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = annotation_col,
         show_rownames = FALSE,
         fontsize_col = 10,
         main = "Heatmap of Significant Proteins")

