library(DEP2)
library(dplyr)
library(ggplot2)

# Define the path to your CSV file
file_path <- "Path/to/file.csv"

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

# Extract row data
row_data <- as.data.frame(rowData(dep_ptm))
# Remove peptides with "_BOVINE" or "_HUMAN" in their name
row_data <- row_data[!grepl("_BOVINE|_HUMAN", row_data$Accession), ]

# Check the structure of row_data to find the correct columns
str(row_data)

library(VennDiagram)
library(grid)

# Extract all relevant comparison columns
comparison_cols <- grep("_vs_H25_significant", colnames(row_data), value = TRUE)

# Initialize a list to store downregulated PTMs
downregulated_list <- list()

# Loop through each comparison
for (col in comparison_cols) {
  condition <- sub("_vs_H25_significant", "", col)
  diff_col <- paste0(condition, "_vs_H25_diff")
  
  downregulated <- row_data[row_data[[col]] == TRUE & row_data[[diff_col]] < -1.5, "name"]
  downregulated_list[[condition]] <- downregulated
}

# Extract only "37" and "25" groups
groups_37 <- grep("37$", names(downregulated_list), value = TRUE)
groups_25 <- grep("25$", names(downregulated_list), value = TRUE)

# Create separate lists
down_37 <- downregulated_list[groups_37]
down_25 <- downregulated_list[groups_25]

# Flatten all peptide lists
all_25 <- unique(unlist(down_25))
all_37 <- unique(unlist(down_37))

# Peptides shared across all "25" groups
shared_25 <- Reduce(intersect, down_25)

# Peptides shared between any two or more "37" groups
shared_37 <- names(which(table(unlist(down_37)) > 1))

# Peptides shared between "25" and "37"
shared_25_37 <- intersect(all_25, all_37)

# General shared pool (excluding core 25-shared)
shared_general <- setdiff(union(shared_25_37, shared_37), shared_25)

# Unique peptides for each "37" group
unique_37_cleaned <- lapply(down_37, function(peps) {
  setdiff(peps, union(shared_25, shared_general))
})

#Saving each file
# Create output directory
output_dir <- "unique_peptides_37_down_detailed"
# Create the directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Loop through each group
for (grp in names(unique_37_cleaned)) {
  # Get unique peptide names for this group
  peptide_names <- unique_37_cleaned[[grp]]
  
  # Subset full row_data for these peptides
  enriched_df <- row_data[row_data$name %in% peptide_names, ]
  
  # Save to CSV
  file_name <- paste0(output_dir, "/", grp, "_unique_peptides_detailed.csv")
  write.csv(enriched_df, file = file_name, row.names = FALSE)
}



library(ggforce)

# Count unique peptides per group
unique_counts <- sapply(unique_37_cleaned, length)

# Total shared peptides
shared_total <- length(shared_25) + length(shared_general)

# Create a data frame for plotting
unique_counts <- sapply(unique_37_cleaned, length)
groups <- names(unique_counts)
counts <- as.numeric(unique_counts)
angles <- seq(0, 2 * pi, length.out = length(groups) + 1)[- (length(groups) + 1)]

df <- data.frame(
  group = groups,
  count = counts,
  angle = angles,
  x = cos(angles) * 2,
  y = sin(angles) * 2
)

# Sunburst_Plot
ggplot() +
  # Central shared circle
  geom_circle(aes(x0 = 0, y0 = 0, r = 1), fill = "#cccccc", alpha = 0.5) +
  annotate("text", x = 0, y = 0, label = paste0("Shared\n", length(shared_25) + length(shared_general)), size = 5) +
  
  # Connecting lines from center to each group
  geom_segment(data = df, aes(x = 0, y = 0, xend = x, yend = y), linetype = "dotted", color = "gray50") +
  
  # Outer unique group circles with color gradient
  geom_circle(data = df, aes(x0 = x, y0 = y, r = 0.5, fill = count), alpha = 0.8) +
  scale_fill_gradient(low = "#a6cee3", high = "#1f78b4") +
  
  # Labels slightly offset from center
  geom_text(data = df, aes(x = x * 1.1, y = y * 1.1, label = paste0(group, "\n", count)), size = 4) +
  
  coord_fixed() +
  theme_void() +
  ggtitle("Starburst Venn: Unique Peptides per '37' Group") +
  theme(legend.position = "none")



# Initialize list to store genes per group
unique_genes_37 <- list()

for (grp in names(unique_37_cleaned)) {
  peptide_names <- unique_37_cleaned[[grp]]
  genes <- row_data[row_data$name %in% peptide_names, "Genes"]
  unique_genes_37[[grp]] <- unique(genes)
}

# Find the maximum number of genes in any group
max_len <- max(sapply(unique_genes_37, length))

# Pad each gene list with NA to match the longest one
padded_genes <- lapply(unique_genes_37, function(g) {
  length(g) <- max_len  # this pads with NAs
  return(g)
})

# Combine into a data frame
gene_df <- as.data.frame(padded_genes, stringsAsFactors = FALSE)

# Save to CSV
write.csv(gene_df, file = "unique_peptides_37_down_detailed/unique_genes_by_group.csv", row.names = FALSE, quote = FALSE)


# Step 1: Count peptides per gene for each group
gene_peptide_counts <- lapply(names(unique_37_cleaned), function(grp) {
  peptides <- unique_37_cleaned[[grp]]
  genes <- row_data[row_data$name %in% peptides, c("name", "Genes")]
  tbl <- table(genes$Genes)
  data.frame(Gene = names(tbl), Count = as.integer(tbl), stringsAsFactors = FALSE)
})
names(gene_peptide_counts) <- names(unique_37_cleaned)

# Step 2: Pad each group’s data frame to the same number of rows
max_len <- max(sapply(gene_peptide_counts, nrow))
padded_tables <- lapply(gene_peptide_counts, function(df) {
  len <- nrow(df)
  if (len < max_len) {
    df[(len + 1):max_len, ] <- NA
  }
  return(df)
})

# Step 3: Rename columns and combine into one data frame
named_tables <- mapply(function(df, grp) {
  colnames(df) <- c(paste0(grp, "_Gene"), paste0(grp, "_Count"))
  return(df)
}, padded_tables, names(padded_tables), SIMPLIFY = FALSE)

# Combine all into one wide data frame
final_df <- do.call(cbind, named_tables)

# Step 4: Save to CSV
write.csv(final_df, file = "unique_peptides_37_down_detailed/unique_genes_with_counts_separate_columns.csv", row.names = FALSE, quote = FALSE)

library(UpSetR)

# Create the binary matrix from your gene lists
gene_matrix <- fromList(unique_genes_37)

# Explicitly define the group order from the list names
group_order <- rev(names(unique_genes_37))

# Open PDF device
pdf("unique_peptides_37_down_detailed/upset_unique_genes.pdf", width = 12, height = 8)

# Plot the UpSet diagram
upset(gene_matrix,
      sets = group_order,
      keep.order = TRUE,
      order.by = "freq",
      text.scale = 1.3)

# Close the device
dev.off()



# Count how many groups each gene appears in
all_genes <- unlist(unique_genes_37)
gene_group_counts <- table(all_genes)

# For each group, keep only genes that appear in that group and nowhere else
unique_genes_only <- lapply(unique_genes_37, function(g) {
  g[gene_group_counts[g] == 1]
})

# For each group, extract peptides for its unique genes
for (grp in names(unique_genes_only)) {
  # Get unique genes for this group
  genes <- unique_genes_only[[grp]]
  
  # Get peptides in this group
  peptides <- unique_37_cleaned[[grp]]
  
  # Subset row_data for peptides that belong to the unique genes
  df <- row_data[row_data$name %in% peptides & row_data$Genes %in% genes, ]
  
  # Save to CSV
  file_name <- paste0(output_dir, "/", grp, "_unique_genes_with_peptides.csv")
  write.csv(df, file = file_name, row.names = FALSE)
}




# Extract peptides missing in X8537

# Define the groups of interest
group_include <- c("X8237", "X8337", "X8437")
group_exclude <- "X8537"

# Get peptides from each group
peptides_include <- lapply(up_37[group_include], unique)
peptides_exclude <- up_37[[group_exclude]]

# Find peptides common to all three included groups
shared_included <- Reduce(intersect, peptides_include)

# Remove peptides that are also in the excluded group
final_peptides <- setdiff(shared_included, peptides_exclude)

# View result
length(final_peptides)
head(final_peptides)

# Subset row_data for these peptides
subset_df <- row_data[row_data$name %in% final_peptides, ]

# Save to file
write.csv(subset_df, file = "unique_peptides_37_detailed/peptides_X8237_X8337_X8437_not_X8537.csv", row.names = FALSE)





# Create a mapping of peptide → protein
peptide_to_protein <- row_data[, c("name", "Accession")]

# Define groups
group_include <- c("X8237", "X8337", "X8437")
group_exclude <- "X8537"

# Peptides per group
peptides_by_group <- up_37[c(group_include, group_exclude)]

# All peptides across selected groups
all_peptides <- unique(unlist(peptides_by_group))

# Create binary presence matrix
presence_matrix <- sapply(peptides_by_group, function(peps) all_peptides %in% peps)
rownames(presence_matrix) <- all_peptides

# Peptides present in all 3 include groups AND absent in the exclude group
strict_peptides <- rownames(presence_matrix)[
  rowSums(presence_matrix[, group_include]) == length(group_include) &
    presence_matrix[, group_exclude] == FALSE
]

# Map peptides to proteins
strict_protein_map <- peptide_to_protein[peptide_to_protein$name %in% strict_peptides, ]

# Group peptides by protein
protein_peptides <- split(strict_protein_map$name, strict_protein_map$Accession)

# For each protein, get all its peptides from row_data
all_protein_peptides <- split(row_data$name, row_data$Accession)

# Keep only proteins where all their peptides are in the strict_peptides set
strict_proteins <- names(protein_peptides)[
  sapply(names(protein_peptides), function(prot) {
    all_peps <- all_protein_peptides[[prot]]
    all(all_peps %in% strict_peptides)
  })
]

# Subset row_data for these proteins
strict_protein_df <- row_data[row_data$Accession %in% strict_proteins, ]

# Save to CSV
write.csv(strict_protein_df, file = "unique_peptides_37_detailed/strict_unique_proteins_X8237_X8337_X8437_not_X8537.csv", row.names = FALSE)

# Summary
length(strict_proteins)
head(strict_proteins)

