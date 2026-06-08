source("~/SIRT6_db/functions/custom_theme_ggplot2.R")

library(dplyr)

#################################################################################################
######################## Meta-analysis on SIRT6 KO vs WT (only mouse) ########################
#################################################################################################

# Upload table for mouse meta-analysis consisting from only SIRT6 KO vs WT 
meta_ko <- read.csv("~/SIRT6_db/meta_analysis/table_for_meta_analysis_mouse_KO.csv")

# Keep only the genes present in ≥15 studies
gene_stats <- meta_ko %>%
  group_by(gene_id) %>%
  summarise(
    n_studies = n(),
    n_species = n_distinct(organism),
    .groups = "drop"
  )

genes_keep <- gene_stats %>%
  filter(n_studies >= 15) %>%
  pull(gene_id)

meta_ko <- meta_ko %>%
  filter(gene_id %in% genes_keep)

# Look through the final meta_ko table
cat("Total rows:", nrow(meta_ko), "\n")
cat("Unique human genes:", length(unique(meta_ko$gene_id)), "\n")
cat("Unique experiments:", length(unique(meta_ko$experiment_id)), "\n")
cat("Unique species:", length(unique(meta_ko$organism)), "\n")






############################################################
# Multilevel meta-analysis using rma.mv()
############################################################

library(data.table)
library(metafor)

############################################################
# 1. Meta-analysis function (multilevel model)
############################################################

run_meta_mv <- function(logFC, se, experiment_id) {
  
  tryCatch({
    
    # Remove invalid values
    valid_idx <- which(!is.na(logFC) & !is.na(se) & se > 0)
    
    logFC <- logFC[valid_idx]
    se <- se[valid_idx]
    experiment_id <- experiment_id[valid_idx]
    
    if (length(logFC) < 2) {
      stop("Not enough data")
    }
    
    # Multilevel model
    res <- rma.mv(
      yi = logFC,
      V = se^2,
      random = ~ 1 | experiment_id,
      method = "REML",
      control = list(iter.max = 1000, rel.tol = 1e-8)
    )
    
    # Total heterogeneity
    tau2_total <- sum(res$sigma2)
    
    # Mean sampling variance
    mean_vi <- mean(se^2)
    
    # I² (total)
    I2_total <- tau2_total / (tau2_total + mean_vi) * 100
    
    data.table(
      meta_LFC = as.numeric(res$b),
      meta_SE = as.numeric(res$se),
      CI_lower = as.numeric(res$ci.lb),
      CI_upper = as.numeric(res$ci.ub),
      pvalue = as.numeric(res$pval),
      tau2 = tau2_total,                 # total variance across levels
      I2 = I2_total,
      k_experiments = length(logFC)
    )
    
  }, error = function(e) {
    
    data.table(
      meta_LFC = NA_real_,
      meta_SE = NA_real_,
      CI_lower = NA_real_,
      CI_upper = NA_real_,
      pvalue = NA_real_,
      tau2 = NA_real_,
      I2 = NA_real_,
      k_experiments = length(logFC)
    )
  })
}

############################################################
# 2. Prepare input data
############################################################

# Convert to data.table (important for fast grouped computation)
setDT(meta_ko)

############################################################
# 3. Run meta-analysis per gene
############################################################

meta_results_mv <- meta_ko[
  ,
  run_meta_mv(
    log2FoldChange,
    lfcSE,
    experiment_id
  ),
  by = gene_id
]

############################################################
# 4. Diagnostics
############################################################

cat("Failed models (NA meta_LFC):\n")
print(sum(is.na(meta_results_mv$meta_LFC)))

cat("\nSummary of I² (heterogeneity):\n")
print(summary(meta_results_mv$I2))

############################################################
# 5. Multiple testing correction
############################################################

meta_results_mv[, FDR := p.adjust(pvalue, method = "BH")]

############################################################
# 6. Save results
############################################################

write.csv(
  meta_results_mv,
  "~/SIRT6_db/meta_analysis/meta_results_KO_mouse.csv",
  row.names = FALSE
)






######################################## Visualizations #######################################################
meta_results <- read.csv("~/SIRT6_db/meta_analysis/meta_results_KO_mouse.csv")
summary(meta_results$meta_LFC)

# Mark DEGs (up- and down-regulated)
meta_results$diffexpressed[meta_results$meta_LFC > 0.58 & meta_results$FDR < 0.05] <- "UP"
meta_results$diffexpressed[meta_results$meta_LFC < -0.58 & meta_results$FDR < 0.05] <- "DOWN"
# Fill NAs with "Insignificant"
meta_results$diffexpressed[is.na(meta_results$diffexpressed)] <- "Insignificant"

# Count DEGs
counts <- meta_results %>%
  dplyr::count(diffexpressed)

n_up <- counts$n[counts$diffexpressed == "UP"]
n_down <- counts$n[counts$diffexpressed == "DOWN"]
n_all <- n_up + n_down

# Extract gene symbols for significant genes
# 1. upload mouse table with gene symbol column
mouse <- read.csv("~/SIRT6_db/meta_analysis/table_for_meta_analysis_mouse_KO.csv")
# 2. Create mapping table
gene_map <- mouse %>%
  dplyr::select(gene_id, gene_symbol, human_gene_symbol) %>%
  distinct()
# 3. Join with meta results
meta_results <- meta_results %>%
  left_join(gene_map, by = "gene_id")

write.csv(
  meta_results,
  "~/SIRT6_db/meta_analysis/meta_results_KO_mouse.csv",
  row.names = FALSE
)

# save the table only with DEGs
diffgenes <- meta_results %>%
  filter(meta_results$diffexpressed != "Insignificant")

write.csv(diffgenes,
          "~/SIRT6_db/meta_analysis/meta_diffgenes_KO_mouse.csv",
          row.names = FALSE)

# Count the number of significant genes (FDR < 0.05)
siggenes <- meta_results %>%
  filter(meta_results$FDR < 0.05) 

write.csv(siggenes,
          "~/SIRT6_db/meta_analysis/meta_siggenes_KO_mouse.csv",
          row.names = FALSE)





################################################
## I. Volcano plot with heterogeneity (size) ###
################################################

library(ggplot2)
library(ggrepel)

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
    aes(label = paste0(gene_symbol, " (", k_experiments, ")")),
    size = 5,
    max.overlaps = 20
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
    y = "-log10 FDR",
    # title = paste0(
    # "Meta-analysis of SIRT6 KO vs WT in mouse (",
    # n_all, "DEGs: ", n_up, " UP)"
    # )
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
  filename = "~/SIRT6_db/meta_analysis/volcano_plot_SIRT6KOvsWT_mouse.png",
  plot = volcano_plot,
  width = 8,
  height = 5,
  dpi = 600
)






######################################################
###### II. meta_LFC vs heterogeneity (I²) plot ########
######################################################

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
    aes(label = gene_symbol),
    size = 5,
    max.overlaps = 6
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
    breaks = c(-3, -2, -1, 0, 1, 2, 3),
    labels = c("-3", "-2", "-1", "0", "1", "2", "3")
  ) +
  coord_cartesian(xlim = c(-3, 3)) +
  
  # Labels and theme
  labs(
    x = "Meta log2 Fold Change (SIRT6 KO vs WT) in mouse",
    y = expression(Heterogeneity~(I^2)),
    # title = "Conserved and context-specific SIRT6 transcriptional targets"
  ) +
  
  theme_custom() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    legend.position = "right",
    axis.title    = element_text(size = 14),
    legend.title  = element_text(size = 14),
    legend.text   = element_text(size = 12)
  )

library(cowplot)

# Save plot without legend
plot_no_legend <- plot_i2 + theme(legend.position = "none")

ggsave(
  filename = "~/SIRT6_db/meta_analysis/meta_LFC_vs_I2_plot_KO_mouse.png",
  plot = plot_no_legend,
  width = 6,
  height = 5,
  dpi = 600
)

# Extract legend and save as separate figure
legend <- get_legend(plot_i2)

legend_plot <- ggdraw(legend) + 
  theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(
  filename = "~/SIRT6_db/meta_analysis/meta_LFC_vs_I2_legend_KO_nested.png",
  plot = legend_plot,
  width = 2,
  height = 3,
  dpi = 600
)










#############################################
###### 4. Forest plot for key genes #########
#############################################

library(metafor)
library(dplyr)

# Add human_gene_symbol column to meta_ko table
meta_ko <- meta_ko %>%
  left_join(gene_map, by = "human_gene_id")

# Add biological system column to meta_ko table
df <- read.csv("~/SIRT6_db/DE_results/Summary_table_and_circular_plot/summary_DE_table.csv")
# mapping table (exclude GSE102830 experiment id to prevent many-to-many mapping)
system_map <- df %>%
  dplyr::select(experiment_id, biological_system) %>%
  distinct() %>%
  filter(experiment_id != "GSE102830")
# Join the tables
meta_ko <- meta_ko %>%
  left_join(system_map, by = "experiment_id")
# Fill GSE102830 biological system values using stratum values from meta_ko
meta_ko <- meta_ko %>%
  mutate(
    biological_system = ifelse(
      experiment_id == "GSE102830",
      stratum,
      biological_system
    )
  )
# Clean renaming for GSE102830
meta_ko <- meta_ko %>%
  mutate(
    biological_system = case_when(
      experiment_id == "GSE102830" & biological_system == "Brain"  ~ "Nervous",
      experiment_id == "GSE102830" & biological_system == "Heart"  ~ "Cardiovascular",
      experiment_id == "GSE102830" & biological_system == "Kidney" ~ "Metabolic",
      experiment_id == "GSE102830" & biological_system == "Liver"  ~ "Metabolic",
      experiment_id == "GSE102830" & biological_system == "Lung"   ~ "Respiratory",
      experiment_id == "GSE102830" & biological_system == "Muscle" ~ "Musculoskeletal",
      experiment_id == "GSE102830" & biological_system == "Thymus" ~ "Immune",
      TRUE ~ biological_system
    )
  )
unique(meta_ko$biological_system)


# Add a new column to meta_ko with short names for organisms
organism_map <- c("Mus musculus" = "Mus m.")

meta_ko <- meta_ko %>%
  mutate(organism_short = organism_map[organism])

unique(meta_ko$organism_short)

# Define genes for the forest plot
genes <- c(
  "ZNF518A",
  "ACTG2",
  "APLN",
  "MCF2",
  "GAL",
  "CFAP141"
)

# Create a function for constructing a forest plot
plot_forest_gene <- function(gene_name, data) {
  
  ############################################################
  # 1. Prepare data
  ############################################################
  
  gene_data <- data %>%
    filter(human_gene_symbol == gene_name) %>%
    arrange(organism, log2FoldChange)
  
  ############################################################
  # 2. Multilevel meta-analysis (IMPORTANT UPDATE)
  ############################################################
  
  res <- rma.mv(
    yi = gene_data$log2FoldChange,
    V  = gene_data$lfcSE^2,
    random = ~ 1 | experiment_id,
    data = gene_data,
    method = "REML"
  )
  
  ############################################################
  # 3. Labels (WITH biological system)
  ############################################################
  
  labels <- paste0(
    gene_data$organism_short,
    " | ",
    gene_data$experiment_id,
    " | ",
    gene_data$biological_system
  )
  
  ############################################################
  # 4. Margins (explanation above)
  ############################################################
  
  par(mar = c(4, 4, 5, 2))  # more top space for title
  
  ############################################################
  # 5. Forest plot
  ############################################################
  
  forest(
    res,
    slab = labels,
    
    # renamed axis
    xlab = "log2 Fold Change (SIRT6 KO vs WT)",
    
    main = gene_name,
    refline = 0,
    cex = 0.8,
    
    header = "Species | Experiment | System",
    
    mlab = "",
    shade = TRUE
  )
  
  ############################################################
  # 6. Add meta-analysis diamond
  ############################################################
  
  addpoly(
    res,
    row = -1,
    mlab = "",
    col = "firebrick",
    border = "firebrick"
  )
  
  ############################################################
  # 7. Add heterogeneity stats (adapted for rma.mv)
  ############################################################
  
  # Extract variance components
  tau2 <- sum(res$sigma2)
  
  # Approximate I2
  I2 <- 100 * tau2 / (tau2 + mean(gene_data$lfcSE^2))
  
  usr <- par("usr")
  
  text(
    x = usr[1],   # ⬅ move further LEFT
    y = usr[3] + 1,
    pos = 4,
    cex = 0.8,
    bquote(paste(
      "Multilevel model: ",
      I^2, " = ", .(formatC(I2, digits = 1, format = "f")), "%, ",
      tau^2, " = ", .(formatC(tau2, digits = 3, format = "f"))
    ))
  )
}

# run the function for each of the genes from the list
for (g in genes) {
  plot_forest_gene(g, meta_ko)
}







# Gene set enrichment analysis (GSEA)
library(clusterProfiler)

# Create ranked gene list
gene_list <- meta_results %>%
  filter(!is.na(meta_LFC)) %>%
  arrange(desc(meta_LFC)) %>%
  distinct(human_gene_symbol, .keep_all = TRUE) %>%
  dplyr::select(human_gene_symbol, meta_LFC)

gene_vector <- gene_list$meta_LFC
names(gene_vector) <- gene_list$human_gene_symbol

sum(names(gene_vector) == "") # 1
sum(is.na(names(gene_vector)))
# remove this empty value
gene_vector <- gene_vector[names(gene_vector) != ""]
sum(is.na(names(gene_vector)))

# Load MSigDB Hallmark pathways
msig_hallmark <- msigdbr(             # download gene sets
  species = "Homo sapiens", 
  category = "H"                      # select hallmark gene sets
)

hallmark_sets <- msig_hallmark %>%
  dplyr::select(gs_name, gene_symbol)

# Run GSEA
gsea_res <- GSEA(
  geneList = gene_vector,
  TERM2GENE = hallmark_sets,
  pAdjustMethod = "BH",
  verbose = FALSE
)

df <- as.data.frame(gsea_res)



# Visualizations

############################################################
# GSEA visualization (publication-quality)
############################################################

library(dplyr)
library(ggplot2)
library(stringr)

############################################################
# 1. Prepare GSEA dataframe
############################################################

gsea_df <- as.data.frame(gsea_res)

gsea_df <- gsea_df %>%
  mutate(
    ########################################################
    # Core enrichment size → GeneRatio
    ########################################################
    core_size = str_count(core_enrichment, "/") + 1,
    GeneRatio = core_size / setSize,
    
    ########################################################
    # Direction (Activated vs Suppressed)
    ########################################################
    direction = ifelse(
      NES > 0,
      "Activated (NES > 0)",
      "Suppressed (NES < 0)"
    ),
    
    ########################################################
    # Transform FDR → -log10(FDR)
    ########################################################
    logFDR = -log10(p.adjust)
  )

############################################################
# 2. Clean pathway names (IMPORTANT for readability)
############################################################

gsea_df$Description <- gsea_df$Description %>%
  gsub("HALLMARK_", "", .) %>%
  str_replace_all("_", " ")

############################################################
# 3. Select TOP pathways per direction
############################################################

gsea_df <- gsea_df %>%
  group_by(direction) %>%
  slice_max(order_by = abs(NES), n = 8) %>%
  ungroup()

############################################################
# 4. Fix ordering (shared across facets)
############################################################

gsea_df$Description <- factor(
  gsea_df$Description,
  levels = gsea_df %>%
    arrange(GeneRatio) %>%
    pull(Description)
)

############################################################
# 5. Plot
############################################################

ggplot(gsea_df, aes(x = GeneRatio, y = Description)) +
  
  ########################################################
# Points
########################################################
geom_point(
  aes(size = setSize, color = logFDR),
  alpha = 0.9
) +
  
  ########################################################
# Facets (Activated vs Suppressed)
########################################################
facet_wrap(~ direction, scales = "free_x") +
  
  ########################################################
# Color scale (significance)
########################################################
scale_color_gradient(
  low = "lightblue",
  high = "firebrick",
  name = "-log10(FDR)"
) +
  
  ########################################################
# Size scale (gene set size)
########################################################
scale_size_continuous(
  name = "Gene set size",
  range = c(3, 8)
) +
  
  ########################################################
# Labels
########################################################
labs(
  x = "Gene Ratio",
  y = "Pathway",
  title = "GSEA of SIRT6 targets: activated and suppressed pathways"
) +
  
  ########################################################
# Theme
########################################################
theme_minimal() +
  
  theme(
    strip.text = element_text(face = "bold", size = 12),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )



# Enrichment plot
# Visualizing the activation of chronic inflammation
p_inflammation <- gseaplot2(
  gsea_res,
  geneSetID = c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", 
                "HALLMARK_INFLAMMATORY_RESPONSE", 
                "HALLMARK_INTERFERON_GAMMA_RESPONSE"),
  title = "SIRT6 loss triggers chronic inflammation in mouse",
  color = c("#B2182B", "#D6604D", "#F4A582"),
  pvalue_table = TRUE,
  ES_geom = "line" # Cleaner look for multiple pathways
)

p_inflammation




# Visualizing the suppression of cell cycle/proliferation
p_cell_cycle <- gseaplot2(
  gsea_res,
  geneSetID = c("HALLMARK_E2F_TARGETS", 
                "HALLMARK_G2M_CHECKPOINT"),
  title = "SIRT6 deficiency leads to mouse cell cycle arrest",
  color = c("#08306B", "#2171B5"), # Deep Blue palette
  pvalue_table = TRUE
)
p_cell_cycle




# Generate the grouping for Structural/Fibrotic changes
p_structural <- gseaplot2(
  gsea_res,
  geneSetID = c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", 
                "HALLMARK_APICAL_JUNCTION"),
  title = "SIRT6 loss induces tissue remodeling",
  color = c("#B2182B", "#D6604D"), # Emerald and Navy contrast
  pvalue_table = TRUE,
  base_size = 11,      # Adjusts text size for better fit
  rel_heights = c(1.5, 0.5, 1) # Balances the ES curve, the rug, and the metric
)

print(p_structural)




# Enrichment of hallmarks of aging
library(dplyr)
library(ggplot2)

# 1. Define the Mapping Dictionary (as established)
mouse_aging_mapper <- c(
  # --- INFLAMMAGING ---
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"          = "Chronic inflammation",
  "HALLMARK_INFLAMMATORY_RESPONSE"            = "Chronic inflammation",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE"        = "Chronic inflammation",
  "HALLMARK_ALLOGRAFT_REJECTION"              = "Chronic inflammation",
  
  # --- CELLULAR SENESCENCE ---
  "HALLMARK_E2F_TARGETS"                      = "Cellular senescence",
  "HALLMARK_G2M_CHECKPOINT"                   = "Cellular senescence",
  
  # --- MITOCHONDRIAL DYSFUNCTION ---
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION"        = "Mitochondrial dysfunction",
  
  # --- DEREGULATED NUTRIENT SENSING ---
  "HALLMARK_HYPOXIA"                          = "Deregulated nutrient sensing",
  
  # --- GENOMIC INSTABILITY & TISSUE INTEGRITY ---
  "HALLMARK_APOPTOSIS"                        = "Genomic instability",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" = "Loss of tissue integrity",
  "HALLMARK_APICAL_JUNCTION"                  = "Loss of tissue integrity"
)


# 2. Apply the mapping to ALL results and clean names
all_aging_results <- gsea_res@result %>%
  mutate(Aging_Hallmark = mouse_aging_mapper[ID]) %>%
  filter(!is.na(Aging_Hallmark)) %>%
  mutate(Description = gsub("HALLMARK_", "", Description) %>% gsub("_", " ", .))

# 3. Create the Faceted Plot colored by P-Adjusted
p_aging <- ggplot(all_aging_results, aes(x = NES, y = reorder(Description, NES), fill = p.adjust)) +
  geom_bar(stat = "identity") +
  # Facet by the Hallmark of Aging category
  facet_grid(Aging_Hallmark ~ ., scales = "free_y", space = "free_y") +
  # Significance scale: Red (Significant) to Blue (Less Significant)
  scale_fill_gradient(low = "#E64B35FF", high = "#4DBBD5FF", name = "p.adj") +
  theme_minimal() +
  labs(
    title = "Hallmarks of aging profile in mouse",
    x = "Normalized enrichment score (NES)",
    y = NULL
  ) +
  theme(
    strip.text.y = element_text(angle = 0, face = "bold", size = 9),
    axis.text.y = element_text(size = 8),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.8, "lines"),
    plot.title = element_text(face = "bold", size = 14)
  )

print(p_aging)
