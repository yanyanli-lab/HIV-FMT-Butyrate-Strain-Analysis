# Title: Figure 5 - Genomic Blueprints of the 15 Key Drivers
# Description: This script visualizes the functional gene profiles (Butyrate, 
# SOD, and Detoxification) of the top 15 active MAGs using a clustered heatmap.

library(dplyr)
library(tidyr)
library(pheatmap)
library(readxl)

# --- 1. Data Loading ---
# Ensure these files are present in your working directory
mag_abund <- read.csv("MASTER_MAGs_Abundance_Matrix.csv", row.names = 1, check.names = FALSE)
bin_ko_counts <- read.csv("GLOBAL_BIN_KO_COUNT_MATRIX.csv", row.names = 1, check.names = FALSE)
df_meta <- read_excel("sample_mapping.xlsx")

# --- 2. Targeted Selection: Identify Top 15 Burst Drivers (Consistent with Fig 4) ---
target_kos <- c("K01034", "K00929", "K01035")
all_colnames <- colnames(bin_ko_counts)
matched_cols <- all_colnames[unlist(lapply(target_kos, function(x) grep(x, all_colnames)))]
common_bins <- intersect(rownames(mag_abund), rownames(bin_ko_counts))

# Calculate Genomic Load and filter FMT group samples at Week 8
bin_gene_load <- rowSums(bin_ko_counts[common_bins, matched_cols, drop=FALSE])
w08_sample_ids <- df_meta %>% 
    filter(GROUP == "FMT" & (Timepoint == "W08" | Timepoint == "W07")) %>% 
    pull(1)

# Rank MAGs by realized power specifically at the W08 peak
w08_mean_abund <- rowMeans(mag_abund[common_bins, w08_sample_ids, drop=FALSE])
final_15_ids <- names(sort(w08_mean_abund * bin_gene_load, decreasing = TRUE)[1:15])

# --- 3. Gene Feature Extraction and Aggregation ---

# A. Butyrate Production Pathway (K01034, K00929, K01035)
b_cols <- all_colnames[grep("K01034|K00929|K01035", all_colnames)]
b_score <- rowSums(bin_ko_counts[, b_cols, drop=FALSE])

# B. Antioxidant Defense (Superoxide Dismutase, SOD)
# Targeted search for SOD-related KOs (e.g., K2157)
s_cols <- all_colnames[grep("K2157", all_colnames)]
s_score <- rowSums(bin_ko_counts[, s_cols, drop=FALSE])

# C. Detoxification Genes (e.g., K01795)
d_cols <- all_colnames[grep("K01795", all_colnames)]
d_score <- rowSums(bin_ko_counts[, d_cols, drop=FALSE])

# --- 4. Heatmap Matrix Construction ---
heatmap_data <- data.frame(
  `Butyrate Production` = b_score[final_15_ids],
  `Antioxidant (SOD)` = s_score[final_15_ids],
  `Detoxification` = d_score[final_15_ids],
  check.names = FALSE
)

# --- 5. Export and Plotting ---

# Exporting as high-resolution PDF for publication
pdf("Figure_5_Genomic_Heatmap.pdf", width = 7, height = 8)

pheatmap(
  heatmap_data,
  display_numbers = TRUE,       # Show gene copy numbers in cells
  number_format = "%.0f",       # Format as integers
  cluster_cols = FALSE,         # Keep functional columns in order
  cluster_rows = TRUE,          # Cluster strains by similarity
  # Red color ramp to highlight high gene loads
  color = colorRampPalette(c("white", "#FEE0D2", "#DE2D26"))(100),
  main = "Genomic Blueprints of the 15 Drivers (Validated)",
  fontsize_number = 11,
  cellwidth = 70,
  cellheight = 30,
  border_color = "grey80"
)

dev.off()

print("✅ Figure 5 Heatmap has been exported as a PDF!")
