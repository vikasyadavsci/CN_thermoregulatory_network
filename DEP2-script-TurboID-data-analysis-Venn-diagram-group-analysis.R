library(DEP2)
library(dplyr)
# Define the path to your CSV file
file_path <- "5675_SupplementalData_112921_normalized-data-TurboID-renamed.csv"

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
diff_CT25 <- test_diff(se_pg, type = "manual", test  = c("CT25_vs_TID25"), fdr.type = "BH")
diff_CTF25 <- test_diff(se_pg, type = "manual", test  = c("CTF25_vs_TID25"), fdr.type = "BH")
diff_CT37 <- test_diff(se_pg, type = "manual", test  = c("CT37_vs_TID37"), fdr.type = "BH")
diff_CTF37 <- test_diff(se_pg, type = "manual", test  = c("CTF37_vs_TID37"), fdr.type = "BH")

## Add significant rejections for features, based on 
dep_CT25 <- add_rejections(diff_CT25, alpha = 0.05, lfc = 1)
dep_CTF25 <- add_rejections(diff_CTF25, alpha = 0.05, lfc = 1)
dep_CT37 <- add_rejections(diff_CT37, alpha = 0.05, lfc = 1)
dep_CTF37 <- add_rejections(diff_CTF37, alpha = 0.05, lfc = 1)

# Extract row data
row_data_CT25 <- as.data.frame(rowData(dep_CT25))
row_data_CTF25 <- as.data.frame(rowData(dep_CTF25))
row_data_CT37 <- as.data.frame(rowData(dep_CT37))
row_data_CTF37 <- as.data.frame(rowData(dep_CTF37))

# Check the structure of row_data to find the correct columns
str(row_data_CT25)
str(row_data_CTF25)
str(row_data_CT37)
str(row_data_CTF37)


library(VennDiagram)
library(grid)

# Extract upregulated significant results for each group
CT25_upregulated <- row_data_CT25[row_data_CT25$CT25_vs_TID25_significant == TRUE & row_data_CT25$CT25_vs_TID25_diff > 1, "name"]
CTF25_upregulated <- row_data_CTF25[row_data_CTF25$CTF25_vs_TID25_significant == TRUE & row_data_CTF25$CTF25_vs_TID25_diff > 1, "name"]
CT37_upregulated <- row_data_CT37[row_data_CT37$CT37_vs_TID37_significant == TRUE & row_data_CT37$CT37_vs_TID37_diff > 1, "name"]
CTF37_upregulated <- row_data_CTF37[row_data_CTF37$CTF37_vs_TID37_significant == TRUE & row_data_CTF37$CTF37_vs_TID37_diff > 1, "name"]

# Create sets for the Venn diagram
CT25_set <- unique(CT25_upregulated)
CTF25_set <- unique(CTF25_upregulated)
CT37_set <- unique(CT37_upregulated)
CTF37_set <- unique(CTF37_upregulated)

# Check the length of each set to ensure they are populated
cat("Group 1 set length:", length(CT25_set), "\n")
cat("Group 2 set length:", length(CTF25_set), "\n")
cat("Group 3 set length:", length(CT37_set), "\n")
cat("Group 4 set length:", length(CTF37_set), "\n")

# Define the output directory
output_dir <- "Turbo-ID-analysis"
# Create the directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Generate the Venn diagram and save it to a file
venn.plot <- venn.diagram(
  x = list(CT25 = CT25_set, CTF25 = CTF25_set, CT37 = CT37_set, CTF37 = CTF37_set),
  category.names = c("CT25", "CTF25", "CT37", "CTF37"),
  fill = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c"),
  alpha = 0.5,
  cat.pos = 0,
  cat.dist = 0.05,
  lwd = 2,
  filename = NULL
)

# Save the Venn diagram to a file
pdf(file = "Turbo-ID-analysis/venn_diagram.pdf")
grid.draw(venn.plot)
dev.off()

# Display the Venn diagram
grid::grid.draw(venn.plot)

# Extract upregulated genes
upregulated_CT37 <- row_data_CT37 %>%
  filter(`CT37_vs_TID37_p.adj` < 0.05 & `CT37_vs_TID37_diff` > 1) %>%
  dplyr::select(name, Genes, Accession)

upregulated_CTF37 <- row_data_CTF37 %>%
  filter(`CTF37_vs_TID37_p.adj` < 0.05 & `CTF37_vs_TID37_diff` > 1) %>%
  dplyr::select(name, Genes, Accession)

upregulated_CT25 <- row_data_CT25 %>%
  filter(`CT25_vs_TID25_p.adj` < 0.05 & `CT25_vs_TID25_diff` > 1) %>%
  dplyr::select(name, Genes, Accession)

upregulated_CTF25 <- row_data_CTF25 %>%
  filter(`CTF25_vs_TID25_p.adj` < 0.05 & `CTF25_vs_TID25_diff` > 1) %>%
  dplyr::select(name, Genes, Accession)

# Combine genes
upregulated_CT37_CTF37 <- bind_rows(upregulated_CT37, upregulated_CTF37) %>%
  distinct()
upregulated_CT25_CTF25 <- bind_rows(upregulated_CT25, upregulated_CTF25) %>%
  distinct()

# Exclude genes upregulated 25 from the combined gene set
upregulated_CT37_CTF37_exclusive <- upregulated_CT37_CTF37 %>%
  filter(!(name %in% upregulated_CT25_CTF25$name))

# Define the full path for the output file
output_file <- file.path(output_dir, "CT37_CTF37_exclusive_proteins.csv")
# Write the data to CSV
write.csv(upregulated_CT37_CTF37_exclusive, output_file, row.names = FALSE)
cat("File saved to:", output_file, "\n")

# Retrieve sequences using UniProt.ws
library(UniProt.ws)
up <- UniProt.ws(taxId=235443)  # Cneo
ids <- upregulated_CT37_CTF37_exclusive$Accession
# Use mapUniProt to retrieve sequences
seqs <- mapUniProt(
  from = "UniProtKB_AC-ID",
  to = "UniProtKB",
  columns = c("accession", "id", "sequence"),
  query = ids
)

# Assuming your result is in a data frame called `seqs`
# Create FASTA-formatted strings
fasta_lines <- paste0(">sp|", seqs$Entry, "|", seqs$Entry.Name, "\n", seqs$Sequence)
# Define output path
fasta_file <- "Turbo-ID-analysis/CT37_CTF37_exclusive_proteins.fasta"
# Write to file
writeLines(fasta_lines, fasta_file)
cat("FASTA file saved to:", fasta_file, "\n")


# Perform inner join on the "Genes" column to keep shared genes and all values from both datasets
IP_Phos_shared <- upregulated_X37_X37F_exclusive %>%
  inner_join(exclusive_cnb37, by = "Genes")

genes_X37_X37F <- upregulated_X37_X37F_exclusive$Genes
genes_cnb37 <- exclusive_cnb37$Genes
venn.plot <- venn.diagram(
  x = list(X37_X37F = genes_X37_X37F, cnb37 = genes_cnb37),
  category.names = c("X37_X37F", "cnb37"),
  fill = c("red", "green"),
  alpha = 0.5,
  cat.pos = 0,
  cat.dist = 0.05,
  lwd = 2,
  filename = NULL
)
# Save the Venn diagram to a file
pdf(file = "venn_diagram1.pdf")
grid.draw(venn.plot)
dev.off()
# Display the Venn diagram
grid::grid.draw(venn.plot)


# Find the number of unique values in the Genes column of shared_genes
IP_Phos_shared_unique <- IP_Phos_shared %>%
  distinct(Genes)

IP_Phos_shared_unique_list <- unique(IP_Phos_shared_unique$Genes)
print(IP_Phos_shared_unique_list)



