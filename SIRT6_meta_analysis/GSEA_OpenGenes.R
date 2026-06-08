source("~/SIRT6_db/functions/custom_theme_ggplot2.R")

library(clusterProfiler)

############################################################
# 1. Load OpenGenes hallmark-gene mapping
############################################################

opengenes <- read.csv("~/SIRT6_db/meta_analysis/aging_databases/OpenGenesdb_gene2hallmarks.csv")

# Format for GSEA
opengenes_sets <- opengenes %>%
  dplyr::select(gs_name = hallmark, gene_symbol) %>%
  distinct()

cat("Hallmarks:", length(unique(opengenes_sets$gs_name)), "\n")
cat("Total gene-hallmark pairs:", nrow(opengenes_sets), "\n")
cat("Unique genes:", length(unique(opengenes_sets$gene_symbol)), "\n")

# Download the output table from meta-analysis
meta_results <- read.csv("~/SIRT6_db/meta_analysis/meta_results_KO.csv")

# Create ranked gene list (ranking rule: sign(LFC) * -log10(FDR))
gene_list <- meta_results %>%
  filter(!is.na(meta_LFC), !is.na(FDR)) %>%
  # Calculate the Pi-score
  mutate(stat_rank = sign(meta_LFC) * -log10(FDR)) %>%
  arrange(desc(stat_rank)) %>%
  distinct(human_gene_symbol, .keep_all = TRUE) %>%
  dplyr::select(human_gene_symbol, stat_rank)

gene_vector <- gene_list$stat_rank
names(gene_vector) <- gene_list$human_gene_symbol

sum(names(gene_vector) == "") # 1
sum(is.na(names(gene_vector)))
# remove this empty value
gene_vector <- gene_vector[names(gene_vector) != ""]
sum(is.na(names(gene_vector)))

# Check how many OpenGenes genes overlap with your ranked list
cat("Overlap with meta-analysis genes:", 
    sum(unique(opengenes_sets$gene_symbol) %in% names(gene_vector)), 
    "out of", length(unique(opengenes_sets$gene_symbol)), "\n")

############################################################
# 2. Run GSEA with OpenGenes hallmarks of aging
############################################################

set.seed(42)

gsea_opengenes <- GSEA(
  geneList = gene_vector,
  TERM2GENE = opengenes_sets,
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  minGSSize = 5,
  maxGSSize = 500,
  verbose = FALSE,
  seed = TRUE
)

############################################################
# 3. View results
############################################################

gsea_og_ko <- as.data.frame(gsea_opengenes)
print(gsea_og_ko %>% 
        dplyr::select(ID, NES, pvalue, p.adjust, setSize) %>%
        arrange(pvalue))

gsea_og_oe <- as.data.frame(gsea_opengenes)
print(gsea_og_oe %>% 
        dplyr::select(ID, NES, pvalue, p.adjust, setSize) %>%
        arrange(pvalue))




##################### Bubble plot ###############################################################
library(ggplot2)
library(dplyr)
library(tidyr)

############################################################
# 1. Combine KO and OE results
############################################################

ko_data <- gsea_og_ko %>%
  dplyr::select(ID, NES, p.adjust, setSize) %>%
  mutate(Condition = "SIRT6 KO")

oe_data <- gsea_og_oe %>%
  dplyr::select(ID, NES, p.adjust, setSize) %>%
  mutate(Condition = "SIRT6 OE")

plot_data <- bind_rows(ko_data, oe_data) %>%
  mutate(
    logFDR = -log10(p.adjust),
    significant = p.adjust < 0.05,
    Condition = factor(Condition, levels = c("SIRT6 KO", "SIRT6 OE"))
  )

# Order hallmarks by KO NES
hallmark_order <- ko_data %>%
  arrange(NES) %>%
  pull(ID)

plot_data$ID <- factor(plot_data$ID, levels = hallmark_order)

############################################################
# 2. Bubble plot
############################################################

gsea <- ggplot(plot_data, aes(x = Condition, y = ID)) +
  
  geom_point(
    aes(size = setSize, fill = NES),
    shape = 21,
    stroke = 0.3,
    color = "black"
  ) +
  
  # Mark significant with a border
  geom_point(
    data = plot_data %>% filter(significant),
    aes(size = setSize),
    shape = 21,
    fill = NA,
    color = "black",
    stroke = 1.5
  ) +
  
  scale_fill_gradient2(
    low = "steelblue",
    mid = "white",
    high = "firebrick",
    midpoint = 0,
    name = "NES"
  ) +
  
  scale_size_continuous(
    name = "Gene set size",
    range = c(3, 12),
    breaks = c(50, 100, 200, 400)
  ) +
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  theme_custom() +
  theme(
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.caption = element_text(size = 9, hjust = 0),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = "~/SIRT6_db/meta_analysis/GSEA_meta_analysis_OpenGenes.png",
  plot = gsea,
  width = 7,
  height = 4,
  dpi = 600
)













