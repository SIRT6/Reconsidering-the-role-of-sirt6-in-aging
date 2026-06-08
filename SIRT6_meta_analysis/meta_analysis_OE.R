source("~/SIRT6_db/functions/custom_theme_ggplot2.R")

library(dplyr)

meta <- read.csv("~/SIRT6_db/meta_analysis/table_for_meta_analysis.csv")

unique(meta$contrast)

#################################################################################################
############### Meta-analysis on SIRT6 OE/SIRT6 K3R OE vs WT (rma.mv() function) ################
#################################################################################################

# Keep only SIRT6 OE/SIRT6 K3R OE
meta_oe <- meta %>%
  filter(contrast == c("SIRT6.OE_vs_WT", "SIRT6.OE.K3R_vs_WT"))

unique(meta_oe$contrast)
unique(meta_oe$organism) # Drosophila melanogaster, Mus musculus, Homo sapiens

# Check the number of genes present in ≥2 and ≥3 organisms
gene_stats <- meta_oe %>%
  group_by(human_gene_id) %>%
  summarise(
    n_studies = n(),
    n_species = n_distinct(organism),
    .groups = "drop"
  )

cat("\nGenes present in ≥2 species:\n")
print(sum(gene_stats$n_species >= 2))

cat("\nGenes present in ≥3 species:\n")
print(sum(gene_stats$n_species >= 3))

# Keep only genes present in ≥2 species and ≥2 studies
genes_keep <- gene_stats %>%
  filter(n_species >= 2 & n_studies >= 2) %>%
  pull(human_gene_id)

meta_oe <- meta_oe %>%
  filter(human_gene_id %in% genes_keep)

# Look through the final meta_ko table
cat("Total rows:", nrow(meta_oe), "\n")
cat("Unique human genes:", length(unique(meta_oe$human_gene_id)), "\n")
cat("Unique experiments:", length(unique(meta_oe$experiment_id)), "\n")
cat("Unique species:", length(unique(meta_oe$organism)), "\n")





##########################################################################################
# Construct a covariance matrix that shows phylogenetic relationship between the organisms
##########################################################################################
library(ape)

# Upload Newick tree with divergence times
phylo_tree_3 <- read.tree("~/SIRT6_db/meta_analysis/three_species_tree.tre")

# Clean up tip labels for display
phylo_tree_3$tip.label <- gsub("_", " ", phylo_tree_3$tip.label)

# Plot
plot(phylo_tree_3, 
     type = "phylogram",
     edge.width = 2,
     label.offset = 5,
     font = 3,               
     cex = 1.2)

nodelabels(c("686", "87"), 
           cex = 0.8, bg = "lightyellow", frame = "rect")

axisPhylo()
title(xlab = "Divergence time (Mya)")

# Compute the phylogenetic correlation matrix
phylo_cor <- vcv(phylo_tree_3, corr = TRUE)

all(unique(meta_oe$organism) %in% rownames(phylo_cor))







############################################################
# Meta-analysis function (nested multilevel model) #########
############################################################

# Load libraries
library(dplyr)
library(data.table)
# install.packages("metafor")
library(metafor)

run_meta_mv <- function(logFC, se, experiment_id, species, phylo_cor) {
  
  tryCatch({
    
    # Remove invalid values
    valid_idx <- which(!is.na(logFC) & !is.na(se) & se > 0)
    
    logFC <- logFC[valid_idx]
    se <- se[valid_idx]
    experiment_id <- experiment_id[valid_idx]
    species <- species[valid_idx]
    
    if (length(logFC) < 2) {
      stop("Not enough data")
    }
    
    # Nested multilevel meta-analysis model with phylogenetic correction
    res <- rma.mv(
      yi = logFC,
      V = se^2,
      random = ~ 1 | species/experiment_id, # random factors - species and experiment id (the same - list(~ 1 | species, ~ 1 | species:experiment_id))
      R = list(species = phylo_cor),
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
      tau2_phylo = res$sigma2[1],
      tau2_experiment = res$sigma2[2],
      I2 = I2_total,
      k_experiments = length(logFC),
      n_species = uniqueN(species)
    )
    
  }, error = function(e) {
    
    data.table(
      meta_LFC = NA_real_,
      meta_SE = NA_real_,
      CI_lower = NA_real_,
      CI_upper = NA_real_,
      pvalue = NA_real_,
      tau2 = NA_real_,
      tau2_phylo = NA_real_,
      tau2_experiment = NA_real_,
      I2 = NA_real_,
      k_experiments = length(logFC),
      n_species = uniqueN(species)
    )
  })
}


# Convert to data.table (important for fast grouped computation)
setDT(meta_oe)

# Run meta-analysis per gene
meta_results_oe <- meta_oe[
  ,
  run_meta_mv(
    log2FoldChange,
    lfcSE,
    experiment_id,
    organism,
    phylo_cor
  ),
  by = human_gene_id
]


# Diagnostics
cat("Failed models (NA meta_LFC):\n")
print(sum(is.na(meta_results_oe$meta_LFC)))

cat("\nSummary of I² (heterogeneity):\n")
print(summary(meta_results_oe$I2))

cat("\nSpecies coverage:\n")
print(table(meta_results_oe$n_species))


# Multiple testing correction
meta_results_oe[, FDR := p.adjust(pvalue, method = "BH")]


# Save results
write.csv(
  meta_results_oe,
  "~/SIRT6_db/meta_analysis/meta_results_OE.csv",
  row.names = FALSE
)









##################################################################################################
########################################## Visualizations ########################################
##################################################################################################

meta_results_oe <- read.csv("~/SIRT6_db/meta_analysis/meta_results_OE.csv")

######################################
######## Prepare the dataframe #######
######################################

# Mark DEGs (up- and down-regulated)
meta_results_oe$diffexpressed[meta_results_oe$meta_LFC > 0.58 & meta_results_oe$FDR < 0.05] <- "UP"
meta_results_oe$diffexpressed[meta_results_oe$meta_LFC < -0.58 & meta_results_oe$FDR < 0.05] <- "DOWN"
# Fill NAs with "No change"
meta_results_oe$diffexpressed[is.na(meta_results_oe$diffexpressed)] <- "Insignificant"

# Extract gene symbols for significant genes
# 1. upload orthologs table with gene symbol column
orthologs <- read.csv("~/SIRT6_db/meta_analysis/ortholog_map_1to1.csv")
# 2. Create mapping table
gene_map <- orthologs %>%
  dplyr::select(human_gene_id, human_gene_symbol) %>%
  distinct()
# 3. Join with meta results
meta_results_oe <- meta_results_oe %>%
  left_join(gene_map, by = "human_gene_id")

# save the output table from meta-analysis
write.csv(meta_results_oe,
          "~/SIRT6_db/meta_analysis/meta_results_OE.csv",
          row.names = FALSE)

# save the table only with DEGs
diffgenes <- meta_results_oe %>%
  filter(meta_results_oe$diffexpressed != "Insignificant")
write.csv(diffgenes,
          "~/SIRT6_db/meta_analysis/meta_diffgenes_OE_rma_mv_nested.csv",
          row.names = FALSE)

# Count DEGs
counts <- meta_results_oe %>%
  dplyr::count(diffexpressed)

n_up <- counts$n[counts$diffexpressed == "UP"]
n_down <- counts$n[counts$diffexpressed == "DOWN"]
n_all <- n_up + n_down

# Count the number of significant genes (FDR < 0.05)
siggenes <- meta_results_oe %>%
  filter(meta_results_oe$FDR < 0.05)

write.csv(siggenes,
          "~/SIRT6_db/meta_analysis/meta_siggenes_OE_rma_mv_nested.csv",
          row.names = FALSE)






################################################
## I. Volcano plot with heterogeneity (size) ###
################################################

library(ggplot2)
library(ggrepel)

# 1. Prepare data
meta_results_oe <- meta_results_oe %>%
  mutate(
    # Size only for DEGs
    size_var = ifelse(
      diffexpressed == "Insignificant",
      NA,        # will handle separately
      I2
    )
  )


# 2. Volcano plot
volcano_plot <- ggplot(meta_results_oe, aes(x = meta_LFC, y = -log10(FDR))) +
  
  # Threshold lines
  geom_vline(xintercept = c(-0.58, 0.58),
             linetype = "dashed", color = "gray") +
  
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", color = "gray") +
  
  
  # Insignificant genes (small grey dots)
  geom_point(
    data = meta_results_oe %>% filter(diffexpressed == "Insignificant"),
    shape = 16,              # filled dot
    color = "grey70",
    size = 1.5,
    alpha = 0.6
  ) +
  
  
  # DEGs (circles with border only)
  geom_point(
    data = meta_results_oe %>% filter(diffexpressed != "Insignificant"),
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
    data = meta_results_oe %>% filter(diffexpressed != "Insignificant"),
    aes(label = paste0(human_gene_symbol, " (", n_species, ")")),
    size = 5,
    max.overlaps = 7
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
    # "Meta-analysis of SIRT6 OE vs WT (",
    # n_all, " DEGs: ", n_up, " UP, ", n_down, " DOWN)"
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
  filename = "~/SIRT6_db/meta_analysis/volcano_plot_SIRT6OEvsWT_nested.png",
  plot = volcano_plot,
  width = 8,
  height = 5,
  dpi = 1200
)







######################################################
###### II. meta_LFC vs heterogeneity (I²) plot ########
######################################################

# 1. Prepare data
meta_results_oe <- meta_results_oe %>%
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
label_data <- meta_results_oe %>%
  filter(FDR < 0.05)

# 3. Plot
plot_i2 <- ggplot(meta_results_oe, aes(x = meta_LFC, y = I2)) +
  
  
  # Insignificant genes 
  geom_point(
    data = meta_results_oe %>% filter(significance_group == "Insignificant"),
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
    data = meta_results_oe %>% filter(significance_group != "Insignificant"),
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
    size = 5,
    max.overlaps = 5
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
  
  # scale_x_continuous(
    # breaks = c(-3, -2, -1, 0, 1, 2, 3),
    # labels = c("-3", "-2", "-1", "0", "1", "2", "3")
  # ) +
  # coord_cartesian(xlim = c(-3, 3)) +
  
  
  # Labels and theme
  labs(
    x = "Meta log2 Fold Change (SIRT6 OE vs WT)",
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
  filename = "~/SIRT6_db/meta_analysis/meta_LFC_vs_I2_plot_OE_nested.png",
  plot = plot_no_legend,
  width = 7,
  height = 5,
  dpi = 1200
)

# Extract legend and save as separate figure
legend <- get_legend(plot_i2)

legend_plot <- ggdraw(legend) + 
  theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(
  filename = "~/SIRT6_db/meta_analysis/meta_LFC_vs_I2_legend_OE_nested.png",
  plot = legend_plot,
  width = 2,
  height = 3,
  dpi = 600
)










##################################################################################################
############################### Gene set enrichment analysis (GSEA) ##############################
##################################################################################################

library(clusterProfiler)
library(msigdbr)

# Create ranked gene list (ranking rule: sign(LFC) * -log10(FDR))
gene_list <- meta_results_oe %>%
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

# Load MSigDB Hallmark pathways
msig_hallmark <- msigdbr(             # download gene sets
  species = "Homo sapiens", 
  category = "H"                      # select hallmark gene sets
)

hallmark_sets <- msig_hallmark %>%
  dplyr::select(gs_name, gene_symbol)

set.seed(42)

# Run GSEA
gsea_res_oe <- GSEA(
  geneList = gene_vector,
  TERM2GENE = hallmark_sets,
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  verbose = FALSE,
  seed = TRUE
)





#####################################
# GSEA visualizations - dot plot ####
#####################################

# 1. Prepare GSEA dataframe
gsea_oe <- as.data.frame(gsea_res_oe)

gsea_oe <- gsea_oe %>%
  mutate(

    # Core enrichment size → GeneRatio
    core_size = str_count(core_enrichment, "/") + 1,
    GeneRatio = core_size / setSize,
    
    # Direction (Activated vs Suppressed)
    direction = ifelse(
      NES > 0,
      "Activated (NES > 0)",
      "Suppressed (NES < 0)"
    ),
    
    # Transform FDR → -log10(FDR)
    logFDR = -log10(p.adjust)
  )

# 2. Clean pathway names for readability
gsea_oe$Description <- gsea_oe$Description %>%
  gsub("HALLMARK_", "", .) %>%
  str_replace_all("_", " ")

# 3. Fix ordering (shared across facets)
gsea_oe$Description <- factor(
  gsea_oe$Description,
  levels = gsea_oe %>%
    arrange(GeneRatio) %>%
    pull(Description)
)

# Filter only significant pathways
gsea_oe_sig <- gsea_oe %>%
  filter(p.adjust < 0.05)

# Fix ordering by absolute NES
gsea_oe_sig$Description <- factor(
  gsea_oe_sig$Description,
  levels = gsea_oe_sig %>%
    arrange(abs(NES)) %>%
    pull(Description)
)

# 4. Plot
barplot <- ggplot(gsea_oe_sig, aes(x = NES, y = Description, fill = logFDR)) +
  
  geom_bar(stat = "identity", width = 0.7) +
  
  geom_vline(xintercept = 0, linetype = "solid", color = "black", alpha = 0.5) +
  
  annotate("text", x = max(gsea_oe_sig$NES) * 0.6, y = Inf, 
           label = "Activated", fontface = "bold", 
           size = 4, vjust = -1, color = "black") +
  
  annotate("text", x = min(gsea_oe_sig$NES) * 0.6, y = Inf, 
           label = "Suppressed", fontface = "bold", 
           size = 4, vjust = -1, color = "black") +
  
  scale_x_continuous(
    limits = c(-2, 2),
    breaks = c(-2, -1, 0, 1, 2)
  ) +
  
  coord_cartesian(clip = "off", xlim = c(-2, 2)) +
  
  scale_fill_gradient(
    low = "lightblue",
    high = "firebrick",
    name = expression(-log[10](FDR))
  ) +
  
  labs(
    x = "Normalized enrichment score (NES)",
    y = NULL,
    # title = "GSEA of SIRT6 targets: activated and suppressed pathways"
  ) +
  
  theme_custom() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.margin = margin(t = 20, r = 10, b = 10, l = 10)
  )

# Save
ggsave(
  filename = "~/SIRT6_db/meta_analysis/gsea_barplot_OE.png",
  plot = barplot,
  width = 7,
  height = 5,
  dpi = 600
)





############################################
# GSEA visualizations - enrichment plot ####
############################################

library(enrichplot)

# 1. Visualizing suppression of inflammation, stress and hormonal signaling
p_inflammaging_oe <- gseaplot2(
  gsea_res_oe,
  geneSetID = c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", 
                "HALLMARK_INFLAMMATORY_RESPONSE", 
                "HALLMARK_UV_RESPONSE_DN",
                "HALLMARK_ESTROGEN_RESPONSE_EARLY"),
  title = "SIRT6 suppresses inflammatory, stress, and hormonal signaling programs",
  color = c("#00468BFF", "#2171B5", "#0099B4FF", "#00A087FF"),
  pvalue_table = TRUE
)
p_inflammaging_oe


# 2. Visualizing activation of mitochondrial bioenergetics and DNA damage response
p_metabolism_oe <- gseaplot2(
  gsea_res_oe,
  geneSetID = c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", 
                "HALLMARK_UV_RESPONSE_UP"
                ),
  title = "SIRT6 promotes mitochondrial bioenergetics and DNA damage response",
  color = c("#E64B35FF", "#DC0000FF"), # Warm palette for Activation
  pvalue_table = TRUE
)
p_metabolism_oe



# 3. Visualizing UV response pathways
p_uv_oe <- gseaplot2(
  gsea_res_oe,
  geneSetID = c("HALLMARK_UV_RESPONSE_DN", 
                "HALLMARK_UV_RESPONSE_UP"
  ),
  title = "UV response pathways are bidirectionally regulated by SIRT6",
  color = c("#00468BFF", "#DC0000FF"), # Warm palette for Activation
  pvalue_table = TRUE
)
p_uv_oe






#############################################################################
# GSEA visualizations - comparison enrichment plot of hallmarks of aging ####
#############################################################################

# 1. Mapping dictionary for the pathways from KO and OE GSEA
aging_mapper <- c(
  # --- CHRONIC INFLAMMATION ---
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"          = "Chronic inflammation",
  "HALLMARK_INFLAMMATORY_RESPONSE"             = "Chronic inflammation",
  "HALLMARK_COMPLEMENT"                        = "Chronic inflammation",
  "HALLMARK_IL2_STAT5_SIGNALING"               = "Chronic inflammation",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING"           = "Chronic inflammation",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE"         = "Chronic inflammation",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE"         = "Chronic inflammation",
  
  # --- ALTERED INTERCELLULAR COMMUNICATION ---
  "HALLMARK_ALLOGRAFT_REJECTION"               = "Altered intercellular communication",
  "HALLMARK_TGF_BETA_SIGNALING"                = "Altered intercellular communication",
  "HALLMARK_COAGULATION"                       = "Altered intercellular communication",
  
  # --- CELLULAR SENESCENCE ---
  "HALLMARK_E2F_TARGETS"                       = "Cellular senescence",
  "HALLMARK_G2M_CHECKPOINT"                    = "Cellular senescence",
  "HALLMARK_MYC_TARGETS_V1"                    = "Cellular senescence",
  "HALLMARK_MYC_TARGETS_V2"                    = "Cellular senescence",
  "HALLMARK_MITOTIC_SPINDLE"                   = "Cellular senescence",
  "HALLMARK_P53_PATHWAY"                       = "Cellular senescence",
  
  # --- MITOCHONDRIAL DYSFUNCTION ---
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION"         = "Mitochondrial dysfunction",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"   = "Mitochondrial dysfunction",
  
  # --- METABOLIC ALTERATIONS ---
  "HALLMARK_FATTY_ACID_METABOLISM"             = "Metabolic alterations",
  "HALLMARK_ADIPOGENESIS"                      = "Metabolic alterations",
  "HALLMARK_PEROXISOME"                        = "Metabolic alterations",
  "HALLMARK_XENOBIOTIC_METABOLISM"             = "Metabolic alterations",
  "HALLMARK_CHOLESTEROL_HOMEOSTASIS"           = "Metabolic alterations",
  "HALLMARK_GLYCOLYSIS"                        = "Metabolic alterations",
  "HALLMARK_BILE_ACID_METABOLISM"              = "Metabolic alterations",
  "HALLMARK_HEME_METABOLISM"                   = "Metabolic alterations",
  "HALLMARK_PANCREAS_BETA_CELLS"               = "Metabolic alterations",
  
  # --- DEREGULATED NUTRIENT SENSING ---
  "HALLMARK_HYPOXIA"                           = "Deregulated nutrient sensing",
  "HALLMARK_MTORC1_SIGNALING"                  = "Deregulated nutrient sensing",
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING"          = "Deregulated nutrient sensing",
  
  # --- GENOMIC INSTABILITY ---
  "HALLMARK_APOPTOSIS"                         = "Genomic instability",
  "HALLMARK_UV_RESPONSE_DN"                    = "Genomic instability",
  "HALLMARK_UV_RESPONSE_UP"                    = "Genomic instability",
  "HALLMARK_DNA_REPAIR"                        = "Genomic instability",
  
  # --- LOSS OF PROTEOSTASIS ---
  "HALLMARK_PROTEIN_SECRETION"                 = "Loss of proteostasis",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"         = "Loss of proteostasis",
  
  # --- CELLULAR SIGNALING ---
  "HALLMARK_KRAS_SIGNALING_UP"                 = "Cellular signaling",
  "HALLMARK_KRAS_SIGNALING_DN"                 = "Cellular signaling",
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING"        = "Cellular signaling",
  "HALLMARK_ESTROGEN_RESPONSE_EARLY"           = "Cellular signaling",
  "HALLMARK_ESTROGEN_RESPONSE_LATE"            = "Cellular signaling",
  "HALLMARK_HEDGEHOG_SIGNALING"                = "Cellular signaling",
  "HALLMARK_NOTCH_SIGNALING"                   = "Cellular signaling",
  "HALLMARK_ANDROGEN_RESPONSE"                 = "Cellular signaling",
  
  # --- LOSS OF TISSUE INTEGRITY ---
  "HALLMARK_APICAL_JUNCTION"                   = "Loss of tissue integrity",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" = "Loss of tissue integrity",
  "HALLMARK_APICAL_SURFACE"                    = "Loss of tissue integrity",
  "HALLMARK_ANGIOGENESIS"                      = "Loss of tissue integrity",
  
  # --- STEM CELL EXHAUSTION ---
  "HALLMARK_MYOGENESIS"                        = "Stem cell exhaustion",
  "HALLMARK_SPERMATOGENESIS"                   = "Stem cell exhaustion"
)

# 2. Extract and Combine Data
# We extract ID and NES, then join them to keep all pathways (full_join)
ko_data <- gsea_df %>% dplyr::select(ID, NES_KO = NES, padj_KO = p.adjust) # gsea_df - result of gsea of SIRT6 KO
oe_data <- gsea_oe %>% dplyr::select(ID, NES_OE = NES, padj_OE = p.adjust)

plot_data <- full_join(ko_data, oe_data, by = "ID") %>%
  filter(ID %in% names(aging_mapper)) %>%
  mutate(Hallmark = aging_mapper[ID]) %>%
  pivot_longer(
    cols = c(NES_KO, NES_OE),
    names_to = "Condition",
    values_to = "NES"
  ) %>%
  # Pivot p.adjust the same way
  mutate(
    padj = ifelse(Condition == "NES_KO", padj_KO, padj_OE)
  ) %>%
  mutate(
    Description = gsub("HALLMARK_", "", ID) %>% gsub("_", " ", .),
    Condition = ifelse(Condition == "NES_KO", "SIRT6 KO", "SIRT6 OE"),
    sig_label = ifelse(!is.na(padj) & padj < 0.05, "*", "")
  )

# Change the description column
plot_data <- plot_data %>%
  mutate(
  # 2. Convert to sentence case (e.g., "Tnfa signaling via nfkb")
  Description = stringr::str_to_sentence(Description),
  
  # 3. Restore all specific MSigDB uppercase terms
  Description = stringr::str_replace_all(
    Description, 
    c(
      # Gene acronyms & Pathways
      "Tnfa"     = "TNFA",
      "nfkb"     = "NFKB",
      "Il2"      = "IL2",
      "stat5"    = "STAT5",
      "Il6"      = "IL6",
      "jak"      = "JAK",
      "stat3"    = "STAT3",
      "Tgf"      = "TGF",
      "E2f"      = "E2F",
      "G2m"      = "G2M",
      "Myc"      = "MYC",
      "P53"      = "p53",       # Kept as standard lowercase 'p'
      "Mtorc1"   = "mTORC1",    # Standard capitalization
      "Pi3k"     = "PI3K",
      "mtor"   = "mTOR",
      "akt"      = "AKT",
      "Mtor"     = "mTOR",      # Standard capitalization
      "Uv"       = "UV",
      "Dna"      = "DNA",
      "Kras"     = "KRAS",
      "Wnt"      = "WNT",
      
      # Directional and Version Suffixes
      "dn$"      = "DN",        # Matches 'Dn' only at the very end of text
      "up$"      = "UP",        # Matches 'Up' only at the very end of text
      "v1$"      = "V1",        # Matches 'V1' at the end
      "v2$"      = "V2"         # Matches 'V2' at the end
    )
  )
)



# 3. The 1st variant
plot <- ggplot(plot_data, aes(x = NES, y = reorder(Description, NES), fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8) +
  # Faceting by Hallmark groups the pathways near each other
  facet_grid(Hallmark ~ ., scales = "free_y", space = "free_y") +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", alpha = 0.5) +
  scale_fill_manual(values = c("SIRT6 KO" = "#FF8C00", "SIRT6 OE" = "#556B2F")) +
  theme_custom() +
  labs(
    title = "SIRT6: a conserved bidirectional regulator of the hallmarks of aging",
    subtitle = "SIRT6 KO and SIRT6 OE meta-analysis comparison",
    x = "Normalized Enrichment Score (NES)",
    y = NULL
  ) +
  theme(
    strip.text.y = element_text(angle = 0, face = "bold", size = 8, hjust = 0),
    axis.text.y = element_text(size = 8),
    panel.spacing = unit(0.5, "lines"),
    legend.position = "bottom"
  )

# 4. The 2nd variant
plot <- ggplot(plot_data, aes(x = NES, y = reorder(Description, NES), fill = Condition)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(
    aes(
      label = sig_label,
      hjust = ifelse(NES >= 0, -0.1, 1.1)  # place outside the bar tip
    ),
    size = 3,
    fontface = "bold",
    color = "black",
    vjust = 0.7
  ) +
  facet_grid(Hallmark ~ Condition, scales = "free_y", space = "free_y") +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", alpha = 0.5) +
  scale_fill_manual(values = c("SIRT6 KO" = "#FF8C00", "SIRT6 OE" = "#556B2F")) +
  theme_custom() +
  labs(
    # title = "SIRT6 as a bidirectional regulator of the hallmarks of aging",
    # subtitle = "SIRT6 KO and SIRT6 OE meta-analysis comparison",
    x = "Normalized Enrichment Score (NES)",
    y = NULL
  ) +
  theme(
      # Explicitly set color = "black" for all text elements
      plot.title   = element_text(hjust = 0.5, color = "black"),
      axis.title.x = element_text(color = "black"),
      axis.text.x  = element_text(color = "black"),
      axis.text.y  = element_text(size = 8, color = "black"),
      strip.text.x = element_text(face = "bold", size = 10, color = "black"),
      strip.text.y = element_text(angle = 0, face = "bold", size = 10, hjust = 0, color = "black"),
      
      panel.spacing = unit(0.5, "lines"), 
      legend.position = "none"
  )

ggsave(
  filename = "~/SIRT6_db/meta_analysis/GSEA_KO_and_OE.png",
  plot = plot,
  width = 9,
  height = 8,
  dpi = 600
)






#############################################################################
######## Prepare meta table with GSEA results for the article ###############
#############################################################################







