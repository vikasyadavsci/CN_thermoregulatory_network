library(Rsamtools)
library(ggplot2)
library(dplyr)

bam_dir <- "Path/to/Ribo-seq-bam-files/"
bam_files <- list.files(bam_dir, pattern = "\\.bam$", full.names = TRUE)

# ---- FAST, MEMORY-SAFE READ LENGTH EXTRACTION ----
get_read_lengths <- function(bam_file, chunk = 1e6) {
  bf <- BamFile(bam_file, yieldSize = chunk)
  open(bf)
  
  lengths <- c()
  
  repeat {
    chunk_data <- scanBam(bf, param = ScanBamParam(what = "qwidth"))[[1]]
    if (length(chunk_data$qwidth) == 0) break
    lengths <- c(lengths, chunk_data$qwidth)
  }
  
  close(bf)
  
  tibble(length = lengths, file = basename(bam_file))
}

# Extract lengths for all BAMs
all_lengths <- bind_rows(lapply(bam_files, get_read_lengths))

# ---- COMBINED FACETED PLOT ----
p_all <- ggplot(all_lengths, aes(x = length)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
  facet_wrap(~ file, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Ribo-seq Read Length Distribution (≤ 50 nt)",
    x = "Read Length (nt)",
    y = "Count"
  ) +
  scale_x_continuous(limits = c(20, 50), breaks = seq(20, 50, 5))

print(p_all)

# ---- PER-FILE PDF OUTPUT ----
pdf_path <- file.path(bam_dir, "Ribo_read_length_histograms.pdf")
pdf(pdf_path, width = 7, height = 5)

for (f in unique(all_lengths$file)) {
  df <- filter(all_lengths, file == f)
  
  p <- ggplot(df, aes(x = length)) +
    geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
    scale_x_continuous(limits = c(20, 50), breaks = seq(20, 50, 5)) +
    theme_minimal() +
    labs(
      title = paste("Read Length Distribution:", f),
      x = "Read Length (nt)",
      y = "Count"
    )
  
  print(p)
}

dev.off()

# ---- SAVE ALL RPF HISTOGRAMS ON ONE PAGE ----
combined_pdf <- file.path(bam_dir, "Ribo_read_length_ALL_in_one_page.pdf")

p_all <- ggplot(all_lengths, aes(x = length)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
  facet_wrap(~ file, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Ribo-seq Read Length Distribution (All Samples)",
    x = "Read Length (nt)",
    y = "Count"
  ) +
  scale_x_continuous(limits = c(20, 50), breaks = seq(20, 50, 5))

pdf(combined_pdf, width = 12, height = 8)
print(p_all)
dev.off()

