# Title: Figure 7A - In situ Replication Rates (iRep) at Week 8
# Description: This bar chart visualizes the replication rates of three sentinel 
# butyrate-producing MAGs, indicating active growth during the functional burst.

library(ggplot2)

# --- 1. Data Preparation ---
# Sentinel, Producer, and Helper strains were selected as representative drivers
irep_df <- data.frame(
  Label = c("Flavonifractor\n(Sentinel)", "Acidaminococcus\n(Producer)", "Alistipes\n(Helper)"),
  iRep = c(1.961, 1.898, 1.834),
  Color_Group = c("Sentinel", "Producer", "Helper")
)

# Maintain the specific order for logical storytelling
irep_df$Label <- factor(irep_df$Label, levels = irep_df$Label)

# --- 2. Figure Generation ---
p_irep <- ggplot(irep_df, aes(x = Label, y = iRep, fill = Color_Group)) +
  # Bar chart with defined borders for the Morandi palette
  geom_bar(stat = "identity", width = 0.6, color = "#4D4D4D", linewidth = 0.7) +
  
  # Horizontal baseline (1.0) indicating one replication per cell cycle
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "#636363", linewidth = 0.8) +
  
  # Value labels above each bar for precision
  geom_text(aes(label = iRep), vjust = -0.6, size = 5.5, fontface = "bold", color = "#333333") +
  
  # Color palette: Morandi-style professional scheme
  scale_fill_manual(values = c(
    "Sentinel" = "#8DD3C7", # Celadon Green
    "Producer" = "#FB8072", # Soft Coral
    "Helper"   = "#BEBADA"  # Lavender Purple
  )) +
  
  # Axis and theme settings
  scale_y_continuous(limits = c(0, 2.3), expand = c(0, 0)) +
  theme_classic(base_size = 14) +
  labs(
    title = "In situ Replication Rates (iRep) at Week 8",
    x = "", 
    y = "Index of Replication (iRep)"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(face = "bold", color = "black", size = 11),
    axis.title.y = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15, margin = margin(b = 15))
  )

# --- 3. Export to High-Resolution PDF ---
ggsave("Figure_7A_iRep_Growth.pdf", p_irep, width = 6, height = 5)

# Preview
print(p_irep)
