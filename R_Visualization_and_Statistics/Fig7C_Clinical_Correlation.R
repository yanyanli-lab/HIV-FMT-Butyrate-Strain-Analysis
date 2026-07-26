# Title: Final Statistical Analysis and Clinical Correlation (Figure 7c)
# Description: This script implements the Functional Capacity Weighted Model (FCWM),
# calculates the Realized Power Score, and performs Spearman correlation with Nadir CD4.
# Corresponding to Methods Section 5.5 and 5.7.

# 1. Load Required Libraries
library(dplyr)
library(tidyr)
library(readxl)
library(stringr)
library(ggplot2)

# --- 2. Load Data Tables ---
# Note: Ensure input files are placed in the working directory
mag_abund <- read.csv("MASTER_MAGs_Abundance_Matrix.csv", row.names = 1, check.names = FALSE)
bin_ko_counts <- read.csv("GLOBAL_BIN_KO_COUNT_MATRIX.csv", row.names = 1, check.names = FALSE)
df_meta <- read_excel("sample_mapping.xlsx")
clinical_data <- read_excel("Info patients REFRESH.xlsx")

# --- 3. Implement FCWM Model (Step 5.5) ---
# Target enzymes for the butyrate-producing pathway
target_kos <- c("K01034", "K00929", "K01035")
ko_cols <- colnames(bin_ko_counts)[unlist(lapply(target_kos, function(x) grep(x, colnames(bin_ko_counts))))]
common_bins <- intersect(rownames(mag_abund), rownames(bin_ko_counts))

# Calculate "Genomic Load" for each MAG
bin_gene_load <- rowSums(bin_ko_counts[common_bins, ko_cols, drop=FALSE])

# Calculate "Realized Power Score" (Abundance * Genomic Load)
sample_total_powers <- colSums(mag_abund[common_bins, ] * bin_gene_load, na.rm = TRUE)

# Filter for FMT group at Week 8 (W08)
individual_w08_scores <- df_meta %>%
  mutate(Timepoint = toupper(Timepoint)) %>%
  filter(GROUP == "FMT" & (Timepoint == "W08" | Timepoint == "W07")) %>%
  mutate(Total_Power = sample_total_powers[run_accession]) %>%
  mutate(Patient_ID = str_extract(sample_alias, "P[0-9]+")) %>%
  select(Patient_ID, Total_Power)

# --- 4. Clinical Correlation Analysis (Step 5.7) ---
# Merging Butyrate Power Score with Nadir CD4 data
final_analysis_df <- individual_w08_scores %>%
  left_join(clinical_data, by = c("Patient_ID" = "ID")) %>%
  mutate(Nadir_CD4 = as.numeric(`Nadir CD4 T cell`)) %>%
  filter(!is.na(Nadir_CD4))

# Statistical Testing: Spearman Correlation
cor_res <- cor.test(final_analysis_df$Nadir_CD4, final_analysis_df$Total_Power, method = "spearman")
fit <- lm(Total_Power ~ Nadir_CD4, data = final_analysis_df)
r2_val <- round(summary(fit)$r.squared, 2)

# --- 5. Figure Generation (Figure 7c) ---
# Visualizing the relationship between Historical Nadir CD4 and W08 Butyrate Power
p_corr_nadir <- ggplot(final_analysis_df, aes(x = Nadir_CD4, y = Total_Power)) +
  # Linear regression line with confidence interval (Red style as in the manuscript)
  geom_smooth(method = "lm", color = "#D73027", fill = "#D73027", alpha = 0.15) +
  # Data points
  geom_point(size = 4, color = "#D73027", alpha = 0.7) +
  theme_bw(base_size = 14) +
  labs(
    title = "Clinical Impact: Nadir CD4 vs. Butyrate Power (W08)",
    subtitle = paste0("Spearman's Rho: ", round(cor_res$estimate, 2), " (P = ", round(cor_res$p.value, 4), ")"),
    x = "Historical Nadir CD4 Count (cells/uL)",
    y = "Realized Butyrate Power Score at W08"
  ) +
  # Annotate R-squared value
  annotate("text", x = min(final_analysis_df$Nadir_CD4), y = max(final_analysis_df$Total_Power), 
           label = paste0("R^2 = ", r2_val), hjust = 0, vjust = 1, fontface = "bold")

# --- 6. Save the Final Figure ---
ggsave("Figure_7c_NadirCD4_Correlation.pdf", p_corr_nadir, width = 7, height = 6)

print("✅ Statistical analysis and Figure 7c generation completed!")
