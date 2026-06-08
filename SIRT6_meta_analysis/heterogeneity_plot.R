##################################################################################################
################################## meta_LFC vs heterogeneity (I²) plot ###########################
##################################################################################################

source("~/SIRT6_db/functions/custom_theme_ggplot2.R")

library(dplyr)
library(ggplot2)
library(cowplot)

meta_results <- read.csv("~/SIRT6_db/meta_analysis/meta_results_KO_human.csv")

# 1. Prepare data
meta_results <- meta_results %>%
  mutate(
    significance_group = case_when(
      FDR < 0.05 & meta_LFC > 0.58  ~ "UP",
      FDR < 0.05 & meta_LFC < -0.58 ~ "DOWN",
      FDR < 0.05                    ~ "Significant",
      TRUE                          ~ "Insignificant"
    ),
    
    # Size variable = -log10(FDR), safe + capped
    size_var = -log10(pmax(FDR, 1e-300)),
    size_var = pmin(size_var, 10)
  )


# 2. Labels (only significant genes)
label_data <- meta_results %>%
  filter(FDR < 0.05)

# 3. Plot
plot_i2 <- ggplot(meta_results, aes(x = meta_LFC, y = I2)) +
  
  
  # Insignificant genes 
  geom_point(
    data = meta_results %>% filter(significance_group == "Insignificant"),
    aes(size = size_var, 
        color = significance_group,
        fill = significance_group
    ),
    shape = 21,
    stroke = 0.8,
    alpha = 0.6
  ) +
  
  
  # Significant genes (filled)
  geom_point(
    data = meta_results %>% filter(significance_group != "Insignificant"),
    aes(
      size = size_var,
      color = significance_group,
      fill = significance_group
    ),
    shape = 21,
    stroke = 1,
    alpha = 0.9
  ) +
  
  
  # Threshold lines
  geom_hline(yintercept = 50, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "grey50") +
  
  
  # Labels
  geom_text_repel(
    data = label_data,
    # For multi-species meta-analysis
    # aes(label = human_gene_symbol),
    
    # For species-species meta-analysis
    aes(label = gene_symbol),
    size = 4,
    max.overlaps = 20
  ) +
  
  
  # Color scale (only for significant)
  scale_color_manual(
    values = c(
      "DOWN" = "steelblue",
      "UP" = "firebrick",
      "Significant" = "orange",
      "Insignificant" = "#DFDFDF"
    ),
    name = "DE genes"
  ) +
  
  
  # Fill scale
  scale_fill_manual(
    values = c(
      "DOWN" = "steelblue",
      "UP" = "firebrick",
      "Significant" = "orange",
      "Insignificant" = "#DFDFDF"
    ),
    guide = "none"
  ) +
  
  
  # Size scale
  scale_size_continuous(
    name = expression(-log[10](FDR)),
    range = c(1.5, 6),
    breaks = c(1, 2, 5, 10),
    labels = c("1", "2", "5", "10"),
    trans = "sqrt"
  ) +
  
  scale_x_continuous(
    breaks = c(-7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7),
    labels = c("-7", "-6", "-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4", "5", "6", "7")
  ) +
  coord_cartesian(xlim = c(-7, 7)) +
  
  # Labels and theme
  labs(
    x = "Meta log2 Fold Change",
    y = expression(Heterogeneity~(I^2))
  ) +
  
  theme_custom() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    legend.position = "right",
    axis.title    = element_text(size = 14),
    legend.title  = element_text(size = 14),
    legend.text   = element_text(size = 12)
  )

# Save plot without legend
plot_no_legend <- plot_i2 + theme(legend.position = "none")

ggsave(
  filename = "~/SIRT6_db/meta_analysis/meta_LFC_vs_I2_plot_KO_human.png",
  plot = plot_no_legend,
  width = 5,
  height = 5,
  dpi = 600
)

# Extract legend and save as separate figure
legend <- get_legend(plot_i2)

legend_plot <- ggdraw(legend) + 
  theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(
  filename = "~/SIRT6_db/meta_analysis/meta_LFC_vs_I2_legend_nested.png",
  plot = legend_plot,
  width = 2,
  height = 3,
  dpi = 600
)
