library(DEP2)
library(dplyr)
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

## Test every sample versus a control
diff_pg <- test_diff(se_pg, type = "control", control = "TID25", fdr.type = "BH")

## Test on manual contrasts
diff_pg2 <- test_diff(se_pg, type = "manual", test  = c("X25_vs_X25F"), fdr.type = "BH")

## Add significant rejections for features, based on 
dep_pg <- add_rejections(diff_pg, alpha = 0.05, lfc = 1.5)

## get the significant subset
dep_pg_sig <- get_signicant(dep_pg)
nrow(dep_pg_sig)

## plot the cutoff line
plot_volcano(dep_pg, contrast = "TID37_vs_TID25", adjusted = F,
             add_threshold_line = "intersect", pCutoff  = 0.05, fcCutoff = 1.5)

plot_heatmap(dep_pg)
