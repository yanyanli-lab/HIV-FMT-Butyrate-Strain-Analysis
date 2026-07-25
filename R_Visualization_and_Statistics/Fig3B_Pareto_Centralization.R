# Title: Figure 3B - Functional Centralization at Week 8 (Pareto Analysis)
# Description: This script calculates the cumulative contribution of MAGs 
# to the total butyrate productivity at the Week 8 functional burst.

library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)

# --- 1. Data Loading and Calculation of Realized Power at W08 ---
mag_abund <- read.csv("MASTER_MAGs_Abundance_Matrix.csv", row.names = 1, check.names = FALSE)
bin_ko_counts <- read.csv("GLOBAL_BIN_KO_COUNT_MATRIX.csv", row.names = 1, check.names = FALSE)
df_meta <- read_excel("sample_mapping.xlsx")

# Define target enzymes for butyrate pathway
target_kos <- c("K01034", "K00929", "K01035")
all_colnames <- colnames(bin_ko_counts)
matched_cols <- all_colnames[unlist(lapply(target_kos, function(x) grep(x, all_colnames)))]

# Calculate Genomic Load for each MAG
common_bins <- intersect(rownames(mag_abund), rownames(bin_ko_counts))
bin_gene_load <- rowSums(bin_ko_counts[common_bins, matched_cols, drop=FALSE])

# Filter for FMT group samples at Week 8 (including W07 as peak phase)
w08_sample_ids <- df_meta %>% 
  filter(GROUP == "FMT" & (Timepoint == "W08" | Timepoint == "W07")) %>% 
  pull(1)

# Calculate mean Realized Power across 279 MAGs at W08
w08_mean_abund_all <- rowMeans(mag_abund[common_bins, w08_sample_ids, drop=FALSE])
w08_all_realized_power <- w08_mean_abund_all * bin_gene_load

# --- 2. Pareto Accumulation Data Preparation ---
all_powers_sorted <- sort(w08_all_realized_power, decreasing = TRUE)
cumulative_share <- cumsum(all_powers_sorted) / sum(all_powers_sorted) * 100

plot_df <- data.frame(
  Rank = 1:length(all_powers_sorted), 
  Cumulative = cumulative_share
)

# --- 3. Figure Generation (Figure 3B) ---
p_pareto <- ggplot(plot_df, aes(x = Rank, y = Cumulative)) +
  # Shaded area
  geom_area(fill = "#BC3C29", alpha = 0.1) + 
  # Main curve
  geom_line(color = "#BC3C29", linewidth = 1.2) +
  # Annotation lines for the Top 15 threshold
  geom_vline(xintercept = 15, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = 62.71, linetype = "dashed", color = "grey40") +
  # Key intersection point
  annotate("point", x = 15, y = 62.71, color = "#BC3C29", size = 4) +
  # Core conclusion text
  annotate("text", x = 100, y = 50, 
           label = "Top 15 Drivers\nAccount for 62.71%", 
           color = "#BC3C29", fontface = "bold", size = 5) +
  # Layout styling
  scale_x_continuous(expand = c(0, 2)) +
  scale_y_continuous(expand = c(0, 2), limits = c(0, 102), breaks = seq(0, 100, 20)) +
  theme_bw(base_size = 14) +
  labs(
    title = "Functional Centralization at Week 8",
    x = "MAGs Ranked by Realized Power (n=279)",
    y = "Cumulative Contribution (%)"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# --- 4. Save Final Figure ---
ggsave("Figure_3B_Pareto_Curve.pdf", p_pareto, width = 7, height = 6)

# Preview
print(p_pareto)
