# Title: Figure 4 - Strain-level Temporal Trajectories of the Butyrate-producing Burst
# Description: This script identifies the top 15 MAGs with the highest realized 
# power specifically at the Week 8 peak and tracks their abundance over time.

library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(stringr)

# --- 1. Data Loading ---
# Ensure these files are present in your working directory
mag_abund <- read.csv("MASTER_MAGs_Abundance_Matrix.csv", row.names = 1, check.names = FALSE)
bin_ko_counts <- read.csv("GLOBAL_BIN_KO_COUNT_MATRIX.csv", row.names = 1, check.names = FALSE)
df_meta <- read_excel("sample_mapping.xlsx")
# Note: gtdb_summary is loaded but can be used for further annotation if needed
gtdb_summary <- read.delim("gtdbtk.bac120.summary.tsv", sep="\t", header=TRUE, check.names=FALSE)

# --- 2. Targeted Selection: Identifying Top 15 Drivers of the W08 Burst ---
target_kos <- c("K01034", "K00929", "K01035")
ko_cols <- colnames(bin_ko_counts)[unlist(lapply(target_kos, function(x) grep(x, colnames(bin_ko_counts))))]
common_bins <- intersect(rownames(mag_abund), rownames(bin_ko_counts))
bin_gene_load <- rowSums(bin_ko_counts[common_bins, ko_cols, drop=FALSE])

# Define W08 sample IDs for the FMT group (including W07 as the peak window)
w08_samples <- df_meta %>% 
  filter(GROUP == "FMT" & (Timepoint == "W08" | Timepoint == "W07")) %>% 
  pull(1) 

# Identify the Top 15 MAGs based on mean realized power at Week 8
w08_mean_abund <- rowMeans(mag_abund[common_bins, w08_samples, drop=FALSE])
w08_power_list <- names(sort(w08_mean_abund * bin_gene_load, decreasing = TRUE)[1:15])

# --- 3. Data Cleaning and Preparation ---
mapping_clean <- df_meta %>%
  rename(Sample_ID = 1, Group = GROUP) %>%
  mutate(Timepoint = toupper(Timepoint)) %>%
  mutate(Timepoint = ifelse(Timepoint == "W07", "W08", Timepoint)) %>%
  filter(Timepoint %in% c("W00", "W01", "W08", "W24")) %>%
  mutate(Patient_ID = str_extract(sample_alias, "R_P[0-9]+")) %>%
  mutate(Group = ifelse(tolower(Group) == "fmt", "FMT", "Placebo"))

plot_data <- mag_abund[w08_power_list, ] %>%
  as.data.frame() %>% mutate(MAG_ID = rownames(.)) %>%
  pivot_longer(cols = -MAG_ID, names_to = "Sample_ID", values_to = "Abundance") %>%
  left_join(mapping_clean, by = "Sample_ID") %>%
  filter(!is.na(Group)) %>%
  mutate(Timepoint = factor(Timepoint, levels = c("W00", "W01", "W08", "W24"))) %>%
  mutate(MAG_ID = factor(MAG_ID, levels = w08_power_list))

# --- 4. Calculate Geometric Mean Trend Lines (Log-mean) ---
summary_data <- plot_data %>%
  group_by(MAG_ID, Group, Timepoint) %>%
  summarise(Trend_Val = exp(mean(log(Abundance + 1))) - 1, .groups = 'drop')

# --- 5. Figure Generation (Figure 4) ---
p4 <- ggplot() +
  # Background: Individual patient trajectories (thin lines)
  geom_line(data = plot_data, 
            aes(x = Timepoint, y = Abundance + 1, group = Patient_ID, color = Group), 
            alpha = 0.08, linewidth = 0.15) +
  # Foreground: Mean trend lines (thick lines)
  geom_line(data = summary_data, 
            aes(x = Timepoint, y = Trend_Val + 1, group = Group, color = Group), 
            linewidth = 0.9) +
  # Layout: Facet per MAG with a fixed Y-axis up to 1000
  facet_wrap(~MAG_ID, scales = "fixed", ncol = 5) + 
  scale_y_log10(limits = c(1, 1000), breaks = c(1, 10, 100, 1000)) +
  # Academic color scheme
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
    title = "Strain-level Temporal Trajectories of the Butyrate-producing Burst at Week 8",
    y = "Abundance (log10 TPM + 1)", 
    x = "Weeks Post-FMT"
  )

# --- 6. Final Export ---
ggsave("Figure_4_Burst_Trajectories.pdf", p4, width = 10, height = 8)

# Preview plot
print(p4)
