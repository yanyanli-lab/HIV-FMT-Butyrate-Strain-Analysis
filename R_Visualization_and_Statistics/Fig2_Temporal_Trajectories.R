# Title: Figure 2 - Strain-level Temporal Trajectories of Butyrate-producing MAGs
# Description: This script identifies the top 15 core strains by total realized potential 
# and visualizes their abundance trajectories across 84 longitudinal samples.

library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(stringr)

# --- 1. Data Loading ---
# Ensure these files are in the same directory as the script
mag_abund <- read.csv("MASTER_MAGs_Abundance_Matrix.csv", row.names = 1, check.names = FALSE)
bin_ko_counts <- read.csv("GLOBAL_BIN_KO_COUNT_MATRIX.csv", row.names = 1, check.names = FALSE)
df_meta <- read_excel("sample_mapping.xlsx")

# --- 2. Core Algorithm: Identify Top 15 Strains by Full Period Potential ---
target_kos <- c("K01034", "K00929", "K01035")
ko_cols <- colnames(bin_ko_counts)[unlist(lapply(target_kos, function(x) grep(x, colnames(bin_ko_counts))))]
common_bins <- intersect(rownames(mag_abund), rownames(bin_ko_counts))
bin_gene_load <- rowSums(bin_ko_counts[common_bins, ko_cols, drop=FALSE])

# Calculate average abundance across all 84 samples
avg_abundance_all <- rowMeans(mag_abund[common_bins, ], na.rm = TRUE)
full_period_power_list <- names(sort(avg_abundance_all * bin_gene_load, decreasing = TRUE)[1:15])

# --- 3. Data Cleaning and Normalization ---
mapping_clean <- df_meta %>%
    rename(Sample_ID = 1, Group = GROUP) %>%
    mutate(Timepoint = toupper(Timepoint)) %>%
    # Merging W07 into W08 for peak phase analysis
    mutate(Timepoint = ifelse(Timepoint == "W07", "W08", Timepoint)) %>%
    filter(Timepoint %in% c("W00", "W01", "W08", "W24")) %>%
    mutate(Patient_ID = str_extract(sample_alias, "R_P[0-9]+")) %>%
    mutate(Group = ifelse(tolower(Group) == "fmt", "FMT", "Placebo"))

plot_data <- mag_abund[full_period_power_list, ] %>%
    as.data.frame() %>% mutate(MAG_ID = rownames(.)) %>%
    pivot_longer(cols = -MAG_ID, names_to = "Sample_ID", values_to = "Abundance") %>%
    left_join(mapping_clean, by = "Sample_ID") %>%
    filter(!is.na(Group)) %>%
    mutate(Timepoint = factor(Timepoint, levels = c("W00", "W01", "W08", "W24"))) %>%
    mutate(MAG_ID = factor(MAG_ID, levels = full_period_power_list))

# --- 4. Calculate Log-mean Trend Lines ---
summary_data <- plot_data %>%
    group_by(MAG_ID, Group, Timepoint) %>%
    summarise(Trend_Val = exp(mean(log(Abundance + 1))) - 1, .groups = 'drop')

# --- 5. Figure Generation (ggplot2) ---
p_full <- ggplot() +
    # Background: Individual trajectories (thin lines)
    geom_line(data = plot_data, 
              aes(x = Timepoint, y = Abundance + 1, group = Patient_ID, color = Group), 
              alpha = 0.08, linewidth = 0.15) +
    # Foreground: Trend lines (thick lines)
    geom_line(data = summary_data, 
              aes(x = Timepoint, y = Trend_Val + 1, group = Group, color = Group), 
              linewidth = 0.9) +
    # Layout and Scales
    facet_wrap(~MAG_ID, scales = "fixed", ncol = 5) + 
    scale_y_log10(limits = c(1, 1000), breaks = c(1, 10, 100, 1000)) +
    # Standard color palette
    scale_color_manual(values = c("FMT" = "#D73027", "Placebo" = "#4575B4")) +
    theme_bw(base_size = 12) +
    theme(
        strip.background = element_rect(fill = "#F0F0F0", color = NA),
        strip.text = element_text(color = "black", size = 7, face = "bold"),
        panel.grid.major = element_line(color = "grey96", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_blank()
    ) +
    labs(
        title = "Longitudinal Dynamics of Top 15 Butyrate-Producing MAGs",
        y = "Abundance (log10 TPM + 1)", 
        x = "Weeks Post-FMT"
    )

# --- 6. Save and Preview ---
ggsave("Figure_2_Temporal_Trajectories.pdf", p_full, width = 10, height = 8)
print(p_full)
