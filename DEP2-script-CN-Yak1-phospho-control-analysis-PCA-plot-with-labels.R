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

# Add rejections
dep_ptm <- DEP2::add_rejections(diff_ptm, alpha = 0.05, lfc = 1.5)

# Extract normalized intensity values from SummarizedExperiment object
pca_data <- t(assay(se_ptm))

# Perform PCA
pca_result <- prcomp(pca_data, scale. = TRUE)

# Extract sample metadata from DEP2 object
metadata <- colData(se_ptm)[, c("ID", "condition", "replicate")]

# Convert to a dataframe for merging
metadata_df <- as.data.frame(metadata)
colnames(metadata_df) <- c("ID", "Condition", "replicate")

# Create PCA dataframe and merge with metadata
pca_df <- data.frame(ID = colnames(se_ptm),
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
plot_heatmap(dep_ptm)

