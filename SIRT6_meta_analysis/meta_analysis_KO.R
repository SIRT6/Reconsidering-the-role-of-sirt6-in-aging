source("~/SIRT6_db/functions/custom_theme_ggplot2.R")

meta <- read.csv("~/SIRT6_db/meta_analysis/table_for_meta_analysis.csv")

unique(meta$contrast)

# Number of unique genes that are not present in human
length(setdiff(unique(meta$human_gene_id), 
               meta$human_gene_id[meta$organism == "Homo sapiens"]))

##################################################################################################
###################### Meta-analysis on SIRT6 KO vs WT only (rma.mv() function) ##################
##################################################################################################

# Keep only SIRT6 KO vs WT experiments
meta_ko <- meta %>%
  filter(contrast == "SIRT6.KO_vs_WT")

unique(meta_ko$contrast)

# Check the number of genes present in ≥2 and ≥3 organisms
gene_stats <- meta_ko %>%
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

# Keep only genes present in ≥3 species and ≥2 studies
genes_keep <- gene_stats %>%
  filter(n_species >= 3 & n_studies >= 2) %>%
  pull(human_gene_id)

meta_ko <- meta_ko %>%
  filter(human_gene_id %in% genes_keep)

# Number of unique genes that are not present in human
length(setdiff(unique(meta_ko$human_gene_id), 
               meta_ko$human_gene_id[meta_ko$organism == "Homo sapiens"]))

# Look through the final meta_ko table
cat("Total rows:", nrow(meta_ko), "\n")
cat("Unique human genes:", length(unique(meta_ko$human_gene_id)), "\n")
cat("Unique experiments:", length(unique(meta_ko$experiment_id)), "\n")
cat("Unique species:", length(unique(meta_ko$organism)), "\n")





##########################################################################################
# Construct a covariance matrix that shows phylogenetic relationship between the organisms
##########################################################################################

library(ape)

# Newick tree with divergence times (in Mya from TimeTree.org)
# Topology: ((((Homo, Macaca), (Mus, Rattus)), Sus), Drosophila)
tree_text <- "((((Homo_sapiens:28.8,Macaca_fascicularis:28.8):58.2,(Mus_musculus:13.1,Rattus_norvegicus:13.1):73.9):7,Sus_scrofa:94):592,Drosophila_melanogaster:686);"
phylo_tree <- read.tree(text = tree_text)

# Clean up tip labels for display
phylo_tree$tip.label <- gsub("_", " ", phylo_tree$tip.label)

# Plot
plot(phylo_tree, 
     type = "phylogram",
     edge.width = 2,
     label.offset = 5,
     font = 3,               
     cex = 1.2)

nodelabels(c("686", "94", "87", "28.8", "13.1"), 
           cex = 0.8, bg = "lightyellow", frame = "rect")

# Add a time scale bar
axisPhylo()
title(xlab = "Divergence time (Mya)")

# Compute the phylogenetic correlation matrix
phylo_cor <- vcv(phylo_tree, corr = TRUE)

all(unique(meta_ko$organism) %in% rownames(phylo_cor))





############################################################
# Meta-analysis function (nested multilevel model) #########
############################################################

# Load libraries
library(dplyr)
library(data.table)
# install.packages("metafor")
library(metafor)

run_meta_mv <- function(logFC, se, experiment_id, species, phylo_cor) {
  
  tryCatch({ # protect the script from crashing
    
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
setDT(meta_ko)


# Run meta-analysis per gene 
meta_results_mv <- meta_ko[
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
print(sum(is.na(meta_results_mv$meta_LFC)))

cat("\nSummary of I² (heterogeneity):\n")
print(summary(meta_results_mv$I2))

cat("\nSpecies coverage:\n")
print(table(meta_results_mv$n_species))


# Multiple testing correction 
meta_results_mv[, FDR := p.adjust(pvalue, method = "BH")]


# Save results
write.csv(
  meta_results_mv,
  "~/SIRT6_db/meta_analysis/meta_results_KO_rma_mv_nested.csv",
  row.names = FALSE
)






##################################################################################################
########################################## Visualizations ########################################
##################################################################################################

meta_results_ko <- read.csv("~/SIRT6_db/meta_analysis/meta_results_KO_rma_mv_nested.csv")

##############################
## Prepare the dataframe #####
##############################

# Mark DEGs (up- and down-regulated)
meta_results_ko$diffexpressed[meta_results_ko$meta_LFC > 0.58 & meta_results_ko$FDR < 0.05] <- "UP"
meta_results_ko$diffexpressed[meta_results_ko$meta_LFC < -0.58 & meta_results_ko$FDR < 0.05] <- "DOWN"
# Fill NAs with "Insignificant"
meta_results_ko$diffexpressed[is.na(meta_results_ko$diffexpressed)] <- "Insignificant"

# save the table only with DEGs
diffgenes <- meta_results_ko %>%
  filter(meta_results_ko$diffexpressed != "Insignificant")
write.csv(diffgenes,
          "~/SIRT6_db/meta_analysis/meta_diffgenes_KO_rma_mv.csv",
          row.names = FALSE)

# Count DEGs
counts <- meta_results_ko %>%
  dplyr::count(diffexpressed)

n_up <- counts$n[counts$diffexpressed == "UP"]
n_down <- counts$n[counts$diffexpressed == "DOWN"]
n_all <- n_up + n_down

# Extract gene symbols for significant genes
# 1. upload orthologs table with gene symbol column
orthologs <- read.csv("~/SIRT6_db/meta_analysis/ortholog_map_1to1.csv")
# 2. Create mapping table
gene_map <- orthologs %>%
  dplyr::select(human_gene_id, human_gene_symbol) %>%
  distinct()
# 3. Join with meta results
meta_results_ko <- meta_results_ko %>%
  left_join(gene_map, by = "human_gene_id")

# Select significant genes for labeling
library(ggrepel)
sig_genes <- meta_results_ko %>%
  filter(diffexpressed != "Insignificant")

# Count the number of significant genes (FDR < 0.05)
df <- meta_results_ko %>%
  filter(meta_results_ko$FDR < 0.05) # the number of rows is 605 => 605 significant genes



################################################
## I. Volcano plot with heterogeneity (size) ###
################################################

library(ggplot2)
library(ggrepel)

# 1. Prepare data
meta_results_ko <- meta_results_ko %>%
  mutate(
    # Size only for DEGs
    size_var = ifelse(
      diffexpressed == "Insignificant",
      NA,        # will handle separately
      I2
    )
  )


# 2. Volcano plot
volcano_plot <- ggplot(meta_results_ko, aes(x = meta_LFC, y = -log10(FDR))) +
  
  
# Threshold lines
geom_vline(xintercept = c(-0.58, 0.58),
           linetype = "dashed", color = "gray") +
  
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", color = "gray") +
  

# Insignificant genes (small grey dots)
geom_point(
  data = meta_results_ko %>% filter(diffexpressed == "Insignificant"),
  shape = 16,              # filled dot
  color = "grey70",
  size = 1.5,
  alpha = 0.6
) +
  
  
# DEGs (circles with border only)
geom_point(
  data = meta_results_ko %>% filter(diffexpressed != "Insignificant"),
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
  data = meta_results_ko %>% filter(diffexpressed != "Insignificant"),
  aes(label = paste0(human_gene_symbol, " (", n_species, ")")),
  size = 3,
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
  x = "meta log2 Fold Change",
  y = "-log10 FDR",
  # title = paste0(
    # "Meta-analysis of SIRT6 KO vs WT (",
    # n_all, " DEGs: ", n_up, " UP, ", n_down, " DOWN)"
  # )
) +
  
  theme_custom()

# 3. Save
ggsave(
  filename = "~/SIRT6_db/meta_analysis/volcano_plot_SIRT6KOvsWT_nested.png",
  plot = volcano_plot,
  width = 8,
  height = 5,
  dpi = 600
)





######################################################
###### II. meta_LFC vs heterogeneity (I²) plot ########
######################################################

# 1. Prepare data
meta_results_ko <- meta_results_ko %>%
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
label_data <- meta_results_ko %>%
  filter(FDR < 0.05)

# 3. Plot
plot_i2 <- ggplot(meta_results_ko, aes(x = meta_LFC, y = I2)) +
  
  
# Insignificant genes 
geom_point(
  data = meta_results_ko %>% filter(significance_group == "Insignificant"),
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
  data = meta_results_ko %>% filter(significance_group != "Insignificant"),
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
  
  
# Labels and theme
labs(
  x = "Meta log2 Fold Change (SIRT6 KO vs WT)",
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
  filename = "~/SIRT6_db/meta_analysis/meta_LFC_vs_I2_plot_KO_nested.png",
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
  filename = "~/SIRT6_db/meta_analysis/meta_LFC_vs_I2_legend_KO_nested.png",
  plot = legend_plot,
  width = 2,
  height = 3,
  dpi = 600
)






#############################################
###### III. Forest plot for key genes #######
#############################################

# 1. Add human_gene_symbol column to meta_ko table
meta_ko <- meta_ko %>%
  left_join(gene_map, by = "human_gene_id")

# 2. Change double experiment ids names 
meta_ko <- meta_ko %>%
  mutate(
    experiment_id = if_else(
      str_detect(experiment_id, ","),
      {
        parts <- str_split(experiment_id, ",")
        first <- map_chr(parts, 1)
        second_last_digit <- map_chr(parts, ~ str_sub(.x[2], -1))
        paste0(first, ",", second_last_digit)
      },
      experiment_id
    )
  )

# 3. Add biological system column to meta_ko table
df <- read.csv("~/SIRT6_db/DE_results/Summary_table_and_circular_plot/summary_DE_table.csv")

sum(is.na(df))

# 4. keep only SIRT6 KO vs WT
df <- df %>%
  filter(contrast == "SIRT6 KO vs WT")

# 5. Create a mapping table 
system_map <- df %>%
  dplyr::select(experiment_id, stratum, biological_system) %>%
  distinct()

# 6. Join the tables by two columns 
meta_ko <- meta_ko %>%
  left_join(system_map, by = c("experiment_id", "stratum"))

# 7. Add a new column to meta_ko with short names for organisms
organism_map <- c(
  "Homo sapiens" = "Homo s.",
  "Mus musculus" = "Mus m.",
  "Rattus norvegicus" = "Rattus n.",
  "Sus scrofa" = "Sus s.",
  "Macaca fascicularis" = "Macaca f.",
  "Drosophila melanogaster" = "Drosophila m."
)
meta_ko <- meta_ko %>%
  mutate(organism_short = organism_map[organism])

# 8. Define genes for the forest plot
genes <- c(
  "APLN",
  "ZNF518A"
)

library(stringr)

# 9. Create a function for constructing a forest plot
plot_forest_gene <- function(gene_name, data, phylo_cor) {
  
  # Prepare data
  gene_data <- data %>%
    filter(human_gene_symbol == gene_name) %>%
    arrange(organism, log2FoldChange)
  
  # Nested multilevel meta-analysis with phylogenetic correction
  res <- rma.mv(
    yi = log2FoldChange,
    V  = lfcSE^2,
    random = ~ 1 | organism/experiment_id,
    R = list(organism = phylo_cor),
    data = gene_data,
    method = "REML",
    control = list(iter.max = 1000, rel.tol = 1e-8)
  )
  
  # Build aligned labels 
  species_width <- max(nchar(gene_data$organism_short), na.rm = TRUE)
  exp_width     <- max(nchar(gene_data$experiment_id), na.rm = TRUE)
  
  # Create aligned labels
  labels <- sprintf(
    paste0("%-", species_width, "s | %-", exp_width, "s | %s"),
    gene_data$organism_short,
    gene_data$experiment_id,
    gene_data$biological_system
  )
  
  # Header aligned the same way
  header <- sprintf(
    paste0("%-", species_width, "s | %-", exp_width, "s | %s"),
    "Species", "Experiment", "System"
  )
  
  # 4. Margins + monospace font 
  par(mar = c(4, 4, 5, 2))
  par(family = "mono")
  
  # Forest plot
  forest(
    res,
    slab = labels,
    xlab = "log2 Fold Change (SIRT6 KO vs WT)",
    main = gene_name,
    refline = 0,
    cex = 1.2,
    header = header,
    mlab = "",
    shade = TRUE,
    alim = c(-4, 4),           # fixed axis range for effect sizes
    xlim = c(-10, 8)           # total plot width (left edge to right edge, includes labels)
  )
  
  # Add meta-analysis diamond
  addpoly(
    res,
    row = -1,
    mlab = "",
    col = "firebrick",
    border = "firebrick"
  )
  
  # Heterogeneity stats (multilevel)
  tau2_total <- sum(res$sigma2)
  tau2_phylo <- res$sigma2[1]
  tau2_experiment <- res$sigma2[2]
  I2 <- 100 * tau2_total / (tau2_total + mean(gene_data$lfcSE^2))
  
  usr <- par("usr")
  line_height = 1.1
  base_y <- usr[3] + 0.3
  
  text(x = usr[1], y = base_y + line_height, pos = 4, cex = 1,
       "Multilevel model:")
  
  text(x = usr[1], y = base_y, pos = 4, cex = 1,
       bquote(I^2 ~ "=" ~ .(formatC(I2, digits = 1, format = "f")) * "%," ~
                tau[phylo]^2 ~ "=" ~ .(formatC(tau2_phylo, digits = 3, format = "e")) * "," ~
                tau[exp]^2 ~ "=" ~ .(formatC(tau2_experiment, digits = 3, format = "e"))))
}


# 10. Save 
out_dir <- "~/SIRT6_db/meta_analysis/"

# APLN
png(
  filename = paste0(out_dir, "forest_APLN.png"),
  width = 9300,
  height = 4000,
  res = 600
)

par(family = "mono")

plot_forest_gene("APLN", meta_ko, phylo_cor)

dev.off()

# ZNF518A
png(
  filename = paste0(out_dir, "forest_ZNF518A.png"),
  width = 9500,
  height = 6000,
  res = 600
)

par(family = "mono")

plot_forest_gene("ZNF518A", meta_ko, phylo_cor)

dev.off()






##################################################################################################
############################### Gene set enrichment analysis (GSEA) ##############################
##################################################################################################

library(clusterProfiler)

# Create ranked gene list (ranking rule: sign(LFC) * -log10(FDR))
gene_list <- meta_results_ko %>%
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

# Run GSEA
gsea_res_ko <- GSEA(
  geneList = gene_vector,
  TERM2GENE = hallmark_sets,
  pAdjustMethod = "BH",
  verbose = FALSE
)




#####################################
# GSEA visualizations - dot plot ####
#####################################

# 1. Prepare GSEA dataframe
gsea_df <- as.data.frame(gsea_res_ko)

gsea_df <- gsea_df %>%
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
gsea_df$Description <- gsea_df$Description %>%
  gsub("HALLMARK_", "", .) %>%
  str_replace_all("_", " ")

# 3. Fix ordering (shared across facets)
gsea_df$Description <- factor(
  gsea_df$Description,
  levels = gsea_df %>%
    arrange(GeneRatio) %>%
    pull(Description)
)

# 4. Plot
ggplot(gsea_df, aes(x = GeneRatio, y = Description)) +
  
# Points
geom_point(
  aes(size = setSize, color = logFDR),
  alpha = 0.9
) +
  
# Facets (Activated vs Suppressed)
facet_wrap(~ direction, scales = "free_x") +
  
# Color scale (significance)
scale_color_gradient(
  low = "lightblue",
  high = "firebrick",
  name = "-log10(FDR)"
) +
  
# Size scale (gene set size)
scale_size_continuous(
  name = "Gene set size",
  range = c(3, 8)
) +
  
# Labels
labs(
  x = "Gene Ratio",
  y = "Pathway",
  title = "GSEA of SIRT6 targets: activated and suppressed pathways"
) +
  
# Theme
theme_minimal() +
  
  theme(
    strip.text = element_text(face = "bold", size = 12),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )





############################################
# GSEA visualizations - enrichment plot ####
############################################

library(enrichplot)

# Top activated pathways
p1 <- gseaplot2(
  gsea_res_k,
  geneSetID = "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  title = "TNFα signaling via NF-κB",
  color = "#E64B35",               # A professional, publication-friendly red
  base_size = 14,                  # Slightly larger font for readability
  pvalue_table = TRUE,             # Overlays the NES, p-value, and p.adjust
  rel_heights = c(2, 0.5, 0.8),    # Gives more vertical space to the enrichment curve
  ES_geom = "line"                 # Can be "line" or "dot"
)

print(p1)


p2 <- gseaplot2(
  gsea_res_k,
  geneSetID = "HALLMARK_HYPOXIA",
  title = "Hypoxia response",
  color = "#E64B35",               # A professional, publication-friendly red
  base_size = 14,                  # Slightly larger font for readability
  pvalue_table = TRUE,             # Overlays the NES, p-value, and p.adjust
  rel_heights = c(2, 0.5, 0.8),    # Gives more vertical space to the enrichment curve
  ES_geom = "line"                 # Can be "line" or "dot"
)

print(p2)



# Top suppressed pathways
p4 <- gseaplot2(
  gsea_res_k,
  geneSetID = "HALLMARK_E2F_TARGETS",
  title = "E2F targets",
  color = "#E64B35",               # A professional, publication-friendly red
  base_size = 14,                  # Slightly larger font for readability
  pvalue_table = TRUE,             # Overlays the NES, p-value, and p.adjust
  rel_heights = c(2, 0.5, 0.8),    # Gives more vertical space to the enrichment curve
  ES_geom = "line"                 # Can be "line" or "dot"
)

print(p4)

p5 <- gseaplot2(
  gsea_res_k,
  geneSetID = "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  title = "Oxidative phosphorylation",
  color = "#E64B35",               # A professional, publication-friendly red
  base_size = 14,                  # Slightly larger font for readability
  pvalue_table = TRUE,             # Overlays the NES, p-value, and p.adjust
  rel_heights = c(2, 0.5, 0.8),    # Gives more vertical space to the enrichment curve
  ES_geom = "line"                 # Can be "line" or "dot"
)

print(p5)



# Overlaying related pathways to show a metabolic shift
p_combined <- gseaplot2(
  gsea_res_ko,
  geneSetID = c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_HYPOXIA"),
  title = "SIRT6 metabolic reprogramming",
  color = c("#E64B35", "#4DBBD5"), # Contrast Blue (OxPhos) vs Red (Hypoxia)
  pvalue_table = TRUE
)

p_combined



# Visualizing the suppression of cell cycle/proliferation
p_cell_cycle <- gseaplot2(
  gsea_res_ko,
  geneSetID = c("HALLMARK_E2F_TARGETS", 
                "HALLMARK_G2M_CHECKPOINT", 
                "HALLMARK_MYC_TARGETS_V1",
                "HALLMARK_MYC_TARGETS_V2"),
  title = "SIRT6 deficiency leads to cell cycle arrest",
  color = c("#08306B", "#2171B5", "#6BAED6", "#6BA"), # Using a distinct palette
  pvalue_table = TRUE
)

p_cell_cycle


# Visualizing the activation of inflammatory signaling (Inflammaging)
p_inflammaging <- gseaplot2(
  gsea_res_ko,
  geneSetID = c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", 
                "HALLMARK_INFLAMMATORY_RESPONSE", 
                "HALLMARK_COMPLEMENT"),
  title = "SIRT6 deficiency leads to enhanced inflammatory signaling",
  color = c("#E64B35FF", "#F39B7FFF", "#DC0000FF"), # Warm palette for activation
  pvalue_table = TRUE
)

print(p_inflammaging)




##################################################################
# GSEA visualizations - enrichment plot of hallmarks of aging ####
##################################################################

# 1. Define the mapping dictionary
aging_mapper <- c(
  # --- INFLAMMAGING & ALTERED COMMUNICATION ---
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"          = "Chronic inflammation",
  "HALLMARK_INFLAMMATORY_RESPONSE"            = "Chronic inflammation",
  "HALLMARK_COMPLEMENT"                       = "Chronic inflammation",
  "HALLMARK_IL2_STAT5_SIGNALING"              = "Chronic inflammation",
  "HALLMARK_ALLOGRAFT_REJECTION"              = "Altered intercellular communication",
  
  # --- CELLULAR SENESCENCE & PROLIFERATION ---
  "HALLMARK_E2F_TARGETS"                      = "Cellular senescence",
  "HALLMARK_G2M_CHECKPOINT"                   = "Cellular senescence",
  "HALLMARK_MYC_TARGETS_V1"                   = "Cellular senescence",
  "HALLMARK_MYC_TARGETS_V2"                   = "Cellular senescence",
  
  # --- METABOLIC ALTERATIONS & MITOCHONDRIAL DYSFUNCTION ---
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION"        = "Mitochondrial dysfunction",
  "HALLMARK_FATTY_ACID_METABOLISM"            = "Metabolic alterations",
  "HALLMARK_ADIPOGENESIS"                     = "Metabolic alterations",
  "HALLMARK_HYPOXIA"                          = "Deregulated nutrient sensing",
  
  # --- GENOMIC INSTABILITY & CELL FATE ---
  "HALLMARK_APOPTOSIS"                        = "Genomic instability",
  "HALLMARK_UV_RESPONSE_DN"                   = "Genomic instability",
  
  # --- SIGNALING & TISSUE INTEGRITY ---
  "HALLMARK_KRAS_SIGNALING_UP"                = "Cellular signaling",
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING"       = "Cellular signaling",
  "HALLMARK_APICAL_JUNCTION"                  = "Loss of tissue integrity",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" = "Loss of tissue integrity"
)


# 2. Apply the mapping to all results and clean names
all_aging_results <- gsea_res_ko@result %>%
  mutate(Aging_Hallmark = aging_mapper[ID]) %>%
  filter(!is.na(Aging_Hallmark)) %>%
  mutate(Description = gsub("HALLMARK_", "", Description) %>% gsub("_", " ", .))

# 3. Create the faceted plot colored by p-adjusted
p_aging <- ggplot(all_aging_results, aes(x = NES, y = reorder(Description, NES), fill = p.adjust)) +
  geom_bar(stat = "identity") +
  # Facet by the Hallmark of Aging category
  facet_grid(Aging_Hallmark ~ ., scales = "free_y", space = "free_y") +
  # Significance scale: Red (Significant) to Blue (Less Significant)
  scale_fill_gradient(low = "#E64B35FF", high = "#4DBBD5FF", name = "p.adj") +
  theme_minimal() +
  labs(
    title = "Hallmarks of Aging profile",
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
