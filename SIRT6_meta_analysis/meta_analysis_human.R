source("~/SIRT6_db/functions/custom_theme_ggplot2.R")

library(dplyr)

meta_ko <- read.csv("~/SIRT6_db/meta_analysis/table_for_meta_analysis_human_KO.csv")

unique(meta$contrast)

#################################################################################################
######################## Meta-analysis on SIRT6 KO vs WT (only human) ########################
#################################################################################################

# Keep only the genes present in ≥4 studies
gene_stats <- meta_ko %>%
  group_by(gene_id) %>%
  summarise(
    n_studies = n(),
    n_species = n_distinct(organism),
    .groups = "drop"
  )

genes_keep <- gene_stats %>%
  filter(n_studies == 7) %>%
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
library(dplyr)

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
      control = list(iter.max = 5000, rel.tol = 1e-8)
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
  "~/SIRT6_db/meta_analysis/meta_results_KO_human.csv",
  row.names = FALSE
)







######################################## Visualizations #######################################################
meta_results <- read.csv("~/SIRT6_db/meta_analysis/meta_results_KO_human.csv")
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
# 1. upload human table with gene symbol column
human <- read.csv("~/SIRT6_db/meta_analysis/table_for_meta_analysis_human_KO.csv")
# 2. Create mapping table
gene_map <- human %>%
  dplyr::select(gene_id, gene_symbol) %>%
  distinct()
# 3. Join with meta results
meta_results <- meta_results %>%
  left_join(gene_map, by = "gene_id")

write.csv(
  meta_results,
  "~/SIRT6_db/meta_analysis/meta_results_KO_human.csv",
  row.names = FALSE
)

# save the table only with DEGs
diffgenes <- meta_results %>%
  filter(meta_results$diffexpressed != "Insignificant")

table(diffgenes$diffexpressed)

write.csv(diffgenes,
          "~/SIRT6_db/meta_analysis/meta_diffgenes_KO_human.csv",
          row.names = FALSE)

# Count the number of significant genes (FDR < 0.05)
siggenes <- meta_results %>%
  filter(meta_results$FDR < 0.05) 

write.csv(siggenes,
          "~/SIRT6_db/meta_analysis/meta_siggenes_KO_human.csv",
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
    aes(label = human_gene_symbol),
    size = 3,
    max.overlaps = 15
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
    x = "meta log2 Fold Change",
    y = "-log10 FDR",
    # title = paste0(
    # "Meta-analysis of SIRT6 KO vs WT in human (",
    # n_all, " DEGs: ", n_up, " UP, ", n_down, " DOWN)"
    # )
  ) +
  
  theme_custom()

# 3. Save
ggsave(
  filename = "~/SIRT6_db/meta_analysis/volcano_plot_SIRT6KOvsWT_human.png",
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
    aes(label = human_gene_symbol),
    size = 3,
    max.overlaps = 10
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
    x = "Meta log2 Fold Change (SIRT6 KO vs WT) in human",
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
  filename = "~/SIRT6_db/meta_analysis/meta_LFC_vs_I2_plot_KO_human.png",
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




num <- diffgenes[diffgenes$I2 == 0, ] # 0 genes
cons <- diffgenes[diffgenes$I2 < 50, ] # 5 genes
spec <- diffgenes[diffgenes$I2 > 50, ] # 27 genes





#####################################################
### 5. Functional enrichment analysis - MSigDB ######
#####################################################

# Load packages
library(dplyr)
library(clusterProfiler)

# install.packages("msigdbr")
library(msigdbr)

library(ggplot2)
library(enrichplot)

# Prepare gene sets
deg_genes_up <- meta_results %>% # up DEGs
  filter(FDR < 0.05 & meta_LFC > 0.58)

deg_genes_down <- meta_results %>% # down DEGs
  filter(FDR < 0.05 & meta_LFC < - 0.58)

sig_genes_up <- meta_results %>% # up significant genes
  filter(FDR < 0.05 & meta_LFC > 0)

sig_genes_down <- meta_results %>% # down significant genes
  filter(FDR < 0.05 & meta_LFC < 0)

background_genes <- meta_results$human_gene_symbol # background genes - all genes tested in meta-analysis

# Load MSigDB Hallmark pathways
msig_hallmark <- msigdbr(             # download gene sets
  species = "Homo sapiens", 
  category = "H"                      # select hallmark gene sets
)

hallmark_sets <- msig_hallmark %>%
  dplyr::select(gs_name, gene_symbol)

# Functional enrichment for DEGs
# Up
sum(unique(deg_genes_up$human_gene_symbol) %in% hallmark_sets$gene_symbol) # only 18, not 45
deg_up_enrich <- enricher(             # enricher uses hypergeometric test
  gene = deg_genes_up$human_gene_symbol,
  universe = background_genes,
  TERM2GENE = hallmark_sets 
)
dotplot(deg_up_enrich, showCategory = 15) +
  ggtitle("Hallmark enrichment of SIRT6 up DEGs")

# Down
sum(unique(deg_genes_down$human_gene_symbol) %in% hallmark_sets$gene_symbol) # only 14, not 33
deg_down_enrich <- enricher(             # enricher uses hypergeometric test
  gene = deg_genes_down$human_gene_symbol,
  universe = background_genes,
  TERM2GENE = hallmark_sets 
)
df_deg_down_enrich <- as.data.frame(deg_down_enrich) # no enriched gene sets 
dotplot(deg_down_enrich, showCategory = 15) +
  ggtitle("Hallmark enrichment of SIRT6 down DEGs")


# Functional enrichment for significant genes
# Up
sig_up_enrich <- enricher(             # enricher uses hypergeometric test
  gene = sig_genes_up$human_gene_symbol,
  universe = background_genes,
  TERM2GENE = hallmark_sets 
)
dotplot(sig_up_enrich, showCategory = 15) +
  ggtitle("Hallmark enrichment of SIRT6 up significant genes")

# Down
sig_down_enrich <- enricher(             # enricher uses hypergeometric test
  gene = sig_genes_down$human_gene_symbol,
  universe = background_genes,
  TERM2GENE = hallmark_sets 
)
dotplot(sig_down_enrich, showCategory = 15) +
  ggtitle("Hallmark enrichment of SIRT6 down significant genes")





##############################################################
### 6. Functional enrichment analysis - GO and KEGG db #######
##############################################################

library(clusterProfiler)
library(org.Hs.eg.db)

# Convert gene symbols -> Entrez IDs (required for GO)
convert_ids <- function(genes) {
  bitr(
    genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )$ENTREZID
}

# Prepare gene lists
deg_up_ids <- convert_ids(deg_genes_up$human_gene_symbol)
deg_down_ids <- convert_ids(deg_genes_down$human_gene_symbol)

sig_up_ids <- convert_ids(sig_genes_up$human_gene_symbol)
sig_down_ids <- convert_ids(sig_genes_down$human_gene_symbol)

bg_ids <- convert_ids(background_genes)

# GO enrichment 
go_deg_up <- enrichGO(
  gene = deg_up_ids,
  universe = bg_ids, # background genes 
  OrgDb = org.Hs.eg.db,
  ont = "BP", # ontology = biological processes 
  pAdjustMethod = "BH",
  readable = TRUE
)

go_deg_down <- enrichGO(
  gene = deg_down_ids,
  universe = bg_ids, # background genes 
  OrgDb = org.Hs.eg.db,
  ont = "BP", # ontology = biological processes 
  pAdjustMethod = "BH",
  readable = TRUE
)

go_sig_up <- enrichGO(
  gene = sig_up_ids,
  universe = bg_ids, # background genes 
  OrgDb = org.Hs.eg.db,
  ont = "BP", # ontology = biological processes 
  pAdjustMethod = "BH",
  readable = TRUE
)

go_sig_down <- enrichGO(
  gene = sig_down_ids,
  universe = bg_ids, # background genes 
  OrgDb = org.Hs.eg.db,
  ont = "BP", # ontology = biological processes 
  pAdjustMethod = "BH",
  readable = TRUE
)


# Visualization
dotplot(go_deg_up, showCategory = 15) +
  ggtitle("GO enrichment: SIRT6 up DEGs")

dotplot(go_deg_down, showCategory = 15) +
  ggtitle("GO enrichment: SIRT6 down DEGs")




# KEGG enrichment
kegg_deg_up <- enrichKEGG(
  gene = deg_up_ids,
  universe = bg_ids,
  organism = "hsa",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2
)

barplot(kegg_deg_up, showCategory = 10) +
  ggtitle("KEGG pathways (Up DEGs)")

kegg_deg_down <- enrichKEGG(
  gene = deg_down_ids,
  universe = bg_ids,
  organism = "hsa",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2
)

barplot(kegg_deg_down, showCategory = 10) +
  ggtitle("KEGG pathways (Down DEGs)")

kegg_sig_up <- enrichKEGG(
  gene = sig_up_ids,
  universe = bg_ids,
  organism = "hsa",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2
)

barplot(kegg_sig_up, showCategory = 10) +
  ggtitle("KEGG pathways (Up significant genes)")

kegg_sig_down <- enrichKEGG(
  gene = sig_down_ids,
  universe = bg_ids,
  organism = "hsa",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2
)

barplot(kegg_deg_down, showCategory = 10) +
  ggtitle("KEGG pathways (Down significant genes)")






###########################################################
### 7. Functional enrichment analysis - Reactome db #######
###########################################################

#BiocManager::install("ReactomePA")
library(ReactomePA)

# Reactome enrichment
react_deg_up <- enrichPathway(
  gene = deg_up_ids,
  universe = bg_ids,
  organism = "human",
  pAdjustMethod = "BH",
  readable = TRUE
)

react_deg_down <- enrichPathway(
  gene = deg_down_ids,
  universe = bg_ids,
  organism = "human",
  pAdjustMethod = "BH",
  readable = TRUE
)


# Visualization
dotplot(react_deg_up, showCategory = 15) +
  ggtitle("Reactome: SIRT6 up DEGs")

dotplot(react_deg_down, showCategory = 15) +
  ggtitle("Reactome: SIRT6 down DEGs")


react_sig_up <- enrichPathway(
  gene = sig_up_ids,
  universe = bg_ids,
  organism = "human",
  pAdjustMethod = "BH",
  readable = TRUE
)

react_sig_down <- enrichPathway(
  gene = sig_down_ids,
  universe = bg_ids,
  organism = "human",
  pAdjustMethod = "BH",
  readable = TRUE
)


# Visualization
dotplot(react_sig_up, showCategory = 15) +
  ggtitle("Reactome: SIRT6 up significant genes")

dotplot(react_sig_down, showCategory = 15) +
  ggtitle("Reactome: SIRT6 down significant genes")





#######################################################
### 8. Merge enrichment results from 3 databases ######
#######################################################

# Prepare MSigDB enrichment results (convert each result into a unified format)
# Use only that results for whom pathways were found
msig_up <- as.data.frame(deg_up_enrich) %>%
  mutate(
    gene_set = "DEGs UP",
    database = "MSigDB"
  )


# Prepare GO enrichment results (convert each result into a unified format)
go_up_df <- as.data.frame(go_deg_up) %>%
  mutate(
  gene_set = "DEGs UP",
  database = "GO"
  )

go_down_df <- as.data.frame(go_deg_down) %>%
  mutate(
    gene_set = "DEGs DOWN",
    database = "GO"
  )


# Prepare KEGG enrichment results (convert each result into a unified format)
kegg_deg_down_df <- as.data.frame(kegg_deg_down) %>%
  mutate(
    gene_set = "DEGs DOWN",
    database = "KEGG"
  )

# Combine everything
combined <- bind_rows(
  msig_up,
  go_up_df,
  go_down_df,
  kegg_deg_down_df
)

# Clean columns
combined <- combined %>%
  mutate(
    logFDR = -log10(p.adjust),
    pathway = Description
  )

# Keep top pathways per group
combined <- combined %>%
  group_by(gene_set, database) %>%
  slice_min(p.adjust, n = 5) %>%
  ungroup()

# Order gene sets
combined$gene_set <- factor(
  combined$gene_set,
  levels = c("DEGs DOWN", "DEGs UP")
)

# Make a 3-dimensional bubble plot
ggplot(combined, aes(x = gene_set, y = pathway)) +
  
  geom_point(
    aes(
      size = logFDR,
      color = database
    )
  ) +
  
  scale_color_manual(
    values = c(
      "MSigDB" = "firebrick",
      "GO" = "darkgreen",
      "KEGG" = "royalblue3"
    )
  ) +
  
  labs(
    x = "Gene set",
    y = "Pathway",
    size = "-log10(FDR)",
    color = "Database",
    title = "Integrated functional enrichment of SIRT6 targets"
  ) +
  
  theme_minimal() +
  
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# Save the table
write.csv(combined,
          "~/SIRT6_db/meta_analysis/combined.csv",
          row.names = FALSE)






########################################################################
#### 9. Hallmarks of aging enrichment analysis using Open Genes db ####
########################################################################

library(dplyr)
library(readr)
library(tidyr)
library(stringr)

# Load the table with genes and hallmarks from Open Genes db
opengenes <- read.csv("~/SIRT6_db/meta_analysis/aging_databases/OpenGenesdb_gene2hallmarks.csv")

# Construct gene sets (hallmark → genes)
aging_sets <- opengenes %>%
  group_by(hallmark) %>%
  summarise(genes = list(unique(gene_symbol))) %>%
  deframe()

############################################################
# 3 Enrichment function (hypergeometric test)
############################################################

run_enrichment <- function(gene_list, set_name, gene_set) {
  
  gene_list <- unique(gene_list)
  
  overlap_genes <- intersect(gene_list, gene_set)
  
  k <- length(overlap_genes)
  M <- length(gene_set)
  N <- length(background_genes)
  n <- length(gene_list)
  
  pval <- phyper(
    q = k - 1,
    m = M,
    n = N - M,
    k = n,
    lower.tail = FALSE
  )
  
  data.frame(
    hallmark = set_name,
    overlap = k,
    input_genes = n,
    hallmark_size = M,
    pvalue = pval,
    overlap_genes = paste(overlap_genes, collapse = "/")
  )
}

############################################################
# 4 Run enrichment across all hallmarks
############################################################

run_all_sets <- function(gene_vector) {
  bind_rows(
    lapply(names(aging_sets), function(name) {
      run_enrichment(
        gene_vector,
        name,
        aging_sets[[name]]
      )
    })
  )
}

# Run for your gene sets
deg_up_res   <- run_all_sets(deg_genes_up$human_gene_symbol)
deg_down_res <- run_all_sets(deg_genes_down$human_gene_symbol)
sig_up_res   <- run_all_sets(sig_genes_up$human_gene_symbol)
sig_down_res <- run_all_sets(sig_genes_down$human_gene_symbol)

# Add labels
deg_up_res$gene_set   <- "DEG UP"
deg_down_res$gene_set <- "DEG DOWN"
sig_up_res$gene_set   <- "SIG UP"
sig_down_res$gene_set <- "SIG DOWN"

# Combine all
results <- bind_rows(
  deg_up_res,
  deg_down_res,
  sig_up_res,
  sig_down_res
)

############################################################
# 5 Multiple testing correction
############################################################

results$FDR <- p.adjust(results$pvalue, method = "BH")

############################################################
# 6 Save results
############################################################

write_csv(
  results,
  "~/SIRT6_db/meta_analysis/opengenes_enrichment_results.csv"
)

############################################################
# 7 Visualization (heatmap)
############################################################

plot_data <- results %>%
  mutate(
    # Calculate -log10(p) for the heatmap fill
    logP = -log10(pvalue),
    logP = ifelse(is.infinite(logP), max(logP[is.finite(logP)]), logP),
    
    # Create the label: Overlap count + asterisk if p < 0.05
    sig_label = ifelse(pvalue < 0.05, 
                       paste0(overlap, "*"), 
                       as.character(overlap))
  )


############################################################
# Order hallmarks by strongest signal
############################################################

hallmark_order <- plot_data %>%
  group_by(hallmark) %>%
  summarise(max_signal = max(logP)) %>%
  arrange(desc(max_signal)) %>%
  pull(hallmark)

plot_data$hallmark <- factor(plot_data$hallmark, levels = hallmark_order)

############################################################
# Order gene sets (important for interpretation)
############################################################

plot_data$gene_set <- factor(
  plot_data$gene_set,
  levels = c("DEG DOWN", "SIG DOWN", "DEG UP", "SIG UP")
)

############################################################
# Plot heatmap
############################################################

ggplot(plot_data, aes(x = gene_set, y = hallmark)) +
  
  # Heatmap tiles
  geom_tile(aes(fill = logP), color = "white") +
  
  # Text labels with asterisks
  geom_text(aes(label = sig_label), size = 3.5, fontface = "bold") +
  
  # Visual separator between UP and DOWN sets
  geom_vline(xintercept = 2.5, color = "grey") +
  
  # Color scale
  scale_fill_gradient(
    low = "white",
    high = "firebrick",
    name = "-log10(pvalue)"
  ) +
  
  labs(
    x = "Gene set",
    y = "Open Genes hallmark",
    title = "Open Genes enrichment of SIRT6 targets",
    caption = "* p < 0.05 (nominal p-value)"
  ) +
  
  theme_minimal() +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.caption = element_text(hjust = 0, face = "italic", size = 9)
  )








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
# remove this NA and empty values
gene_vector <- gene_vector[!is.na(names(gene_vector))]
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
p <- gseaplot2(
  gsea_res,
  geneSetID = c("HALLMARK_E2F_TARGETS", 
                "HALLMARK_G2M_CHECKPOINT", 
                "HALLMARK_APOPTOSIS"),
  title = "SIRT6 deficiency triggers cell cycle arrest and programmed death in human",
  color = c("#D6604D", "#2171B5", "#08306B"),
  pvalue_table = TRUE,
  ES_geom = "line" # Cleaner look for multiple pathways
)

p







# Enrichment of hallmarks of aging
library(dplyr)
library(ggplot2)

# 1. Define the Mapping Dictionary (as established)
human_aging_mapper <- c(
  # --- CELLULAR SENESCENCE (Suppressed in KO) ---
  "HALLMARK_E2F_TARGETS"                       = "Cellular senescence",
  "HALLMARK_G2M_CHECKPOINT"                    = "Cellular senescence",
  
  # --- GENOMIC INSTABILITY (Activated in KO) ---
  "HALLMARK_APOPTOSIS"                         = "Genomic instability",
  
  # --- STEM CELL EXHAUSTION & REGENERATION (Suppressed in KO) ---
  "HALLMARK_MYOGENESIS"                        = "Stem cell exhaustion",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"  = "Loss of tissue integrity",
  
  # --- ALTERED COMMUNICATION / SASP (Activated in KO) ---
  "HALLMARK_PROTEIN_SECRETION"                 = "Altered intercellular communication"
)


# 2. Apply the mapping to ALL results and clean names
all_aging_results <- gsea_res@result %>%
  mutate(Aging_Hallmark = human_aging_mapper[ID]) %>%
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
    title = "Hallmarks of aging profile in human",
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




