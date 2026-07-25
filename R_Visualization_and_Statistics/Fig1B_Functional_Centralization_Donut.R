# Title: Figure 1B - Functional Centralization (Genomic Potential)
# Description: Generates a donut chart showing the contribution of Core 15 consortium
# to the total butyrate synthesis gene pool.

library(ggplot2)
library(dplyr)

# 1. Prepare Data
df_donut <- data.frame(
  Category = c("Core 15 Consortium", "Remaining 264 MAGs"),
  Value = c(40.0, 60.0)
)

# 2. Plotting (Advanced Blue Version)
p_donut_blue <- ggplot(df_donut, aes(x = 2, y = Value, fill = Category)) +
  # Draw the bars with a white border for a clean look
  geom_bar(stat = "identity", width = 1, color = "white", linewidth = 1.5) +
  coord_polar(theta = "y", start = 0) +
  # Set the inner hole (x-axis from 0.5 to 2.5 creates the donut effect)
  xlim(0.5, 2.5) +
  # Color scheme: Deep Blue (#0072B5) + Light Gray
  scale_fill_manual(values = c("Core 15 Consortium" = "#0072B5", "Remaining 264 MAGs" = "#EEEEEE")) +
  theme_void() + 
  # --- Center Text Annotations ---
  # Primary text: 40.0%
  annotate("text", x = 0.5, y = 0, label = "40.0%", 
           size = 11, fontface = "bold", color = "#0072B5") +
  # Secondary text: Genomic Potential (Adjusted x to avoid warnings)
  annotate("text", x = 0.8, y = 0, label = "Genomic Potential", 
           size = 4, color = "grey40", fontface = "italic", vjust = 4) +
  # -----------------------
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15, vjust = -1),
    plot.margin = margin(t = 20, r = 10, b = 20, l = 10)
  ) +
  labs(title = "Functional Centralization")

# Save the final output
ggsave("Figure_1B_Donut_Chart.pdf", p_donut_blue, width = 6, height = 7)

# Preview
print(p_donut_blue)
