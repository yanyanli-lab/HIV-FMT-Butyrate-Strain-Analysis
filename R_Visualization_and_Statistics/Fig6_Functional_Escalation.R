# Title: Figure 6 - Functional Escalation of the Core 15 Consortium
# Description: This bar chart visualizes the non-linear increase in functional contribution 
# of the core strains from genomic potential to real-time productivity during the peak phase.

library(ggplot2)

# --- 1. Data Preparation ---
# Stage labels use concise English to ensure clarity on the X-axis
escalation_data <- data.frame(
  Stage = c("Genomic Potential", "Long-term Average", "Week 8 Peak"),
  Contribution = c(40.0, 47.55, 62.71)
)

# Force the order of the X-axis levels
escalation_data$Stage <- factor(escalation_data$Stage, levels = escalation_data$Stage)

# --- 2. Figure Generation ---
p_escalation <- ggplot(escalation_data, aes(x = Stage, y = Contribution, fill = Stage)) +
  # Main bars with black borders for definition
  geom_bar(stat = "identity", width = 0.5, color = "black", linewidth = 0.6) +
  
  # Guideline segments to illustrate the "escalation" steps
  geom_segment(aes(x = 1.25, y = 40.0, xend = 1.75, yend = 47.55), 
               linetype = "dashed", color = "grey30", linewidth = 0.5) +
  geom_segment(aes(x = 2.25, y = 47.55, xend = 2.75, yend = 62.71), 
               linetype = "dashed", color = "grey30", linewidth = 0.5) +
  
  # Value labels above each bar
  geom_text(aes(label = paste0(Contribution, "%")), vjust = -0.8, size = 5, fontface = "bold") +
  
  # Academic color scheme (Lancet/Journal style)
  scale_fill_manual(values = c("#BDD7EE", "#6BAED6", "#BC3C29")) +
  
  # Axis scaling and theme adjustments
  scale_y_continuous(limits = c(0, 80), expand = c(0, 0)) +
  theme_classic(base_size = 14) +
  labs(
    title = "Functional Escalation of Core 15 Consortium",
    y = "Contribution (%)",
    x = ""
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(face = "bold", color = "black", size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.margin = margin(t = 10, r = 10, b = 20, l = 10)
  )

# --- 3. Export to High-Resolution PDF ---
# Recommended dimensions: 7x5 inches for standard bar charts
ggsave("Figure_6_Functional_Escalation.pdf", p_escalation, width = 7, height = 5)

# Preview
print(p_escalation)
