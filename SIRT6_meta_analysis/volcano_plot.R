##################################################################################################
################################## Volcano plot with heterogeneity (size) ########################
##################################################################################################

source("~/SIRT6_db/functions/custom_theme_ggplot2.R")

library(dplyr)
library(ggrepel)
library(ggplot2)

meta_results <- read.csv("~/SIRT6_db/meta_analysis/meta_results_KO_human.csv")

# 1. Prepare data
meta_results <- meta_results %>%
  mutate(
    # Size only for DEGs
    size_var = ifelse(
      diffexpressed == "Insignificant",
      NA,        # will handle separately
      I2
    )
  )


# 2. Volcano plot
volcano_plot <- ggplot(meta_results, aes(x = meta_LFC, y = -log10(FDR))) +
  
  
  # Threshold lines
  geom_vline(xintercept = c(-0.58, 0.58),
             linetype = "dashed", color = "gray") +
  
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", color = "gray") +
  
  
  # Insignificant genes (small grey dots)
  geom_point(
    data = meta_results %>% filter(diffexpressed == "Insignificant"),
    shape = 16,              # filled dot
    color = "grey70",
    size = 1.5,
    alpha = 0.6
  ) +
  
  
  # DEGs (circles with border only)
  geom_point(
    data = meta_results %>% filter(diffexpressed != "Insignificant"),
    aes(
      size = size_var,
      fill = diffexpressed,
      color = diffexpressed
    ),
    shape = 21,              # filled circle
    stroke = 1.2,           # thickness of border
    alpha = 0.9
  ) +
  
  
  # Labels
  geom_text_repel(
    data = meta_results %>% filter(diffexpressed != "Insignificant"),
    # For multi-species meta-analysis
    # aes(label = paste0(human_gene_symbol, " (", n_species, ")")),
    
    # For species-species meta-analysis
    aes(label = paste0(gene_symbol, " (", k_experiments, ")")),
    size = 4,
    max.overlaps = 10
  ) +
  
  
  # Color scale (border only)
  scale_fill_manual(
    values = c(
      "UP" = "firebrick",
      "DOWN" = "steelblue"
    ),
    name = "DE genes"
  ) +
  
  scale_color_manual(
    values = c("UP" = "firebrick", "DOWN" = "steelblue"),
    name = "DE genes"
  ) + 
  
  
  # Size scale (heterogeneity)
  scale_size_continuous(
    name = "Heterogeneity (I²)",
    range = c(2, 7),
    limits = c(0, 100),
    breaks = c(20, 50, 80),
    labels = c("<50%", "50–80%", ">80%")  # reversed interpretation
  ) +
  
  
  # Labels
  labs(
    x = "Meta log2 Fold Change",
    y = "-log10 FDR"
  ) +
  
  theme_custom() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    legend.position = "right",
    axis.title    = element_text(size = 14),
    legend.title  = element_text(size = 14),
    legend.text   = element_text(size = 12)
  )

# 3. Save
ggsave(
  filename = "~/SIRT6_db/meta_analysis/volcano_plot_SIRT6KOvsWT_human.png",
  plot = volcano_plot,
  width = 7,
  height = 5,
  dpi = 600
)
