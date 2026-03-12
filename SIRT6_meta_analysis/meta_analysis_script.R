meta <- read.csv("~/SIRT6_db/meta_analysis/table_for_meta_analysis.csv")

unique(meta$contrast)

##########################################################################################################################
###################################### Meta-analysis on SIRT6 KO vs WT only ###############################################
##########################################################################################################################

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

# Keep only genes present in ≥2 species and ≥2 studies
genes_keep <- gene_stats %>%
  filter(n_species >= 2 & n_studies >= 2) %>%
  pull(human_gene_id)

meta_ko <- meta_ko %>%
  filter(human_gene_id %in% genes_keep)

# Look through the final meta_ko table
cat("Total rows:", nrow(meta_ko), "\n")
cat("Unique human genes:", length(unique(meta_ko$human_gene_id)), "\n")
cat("Unique experiments:", length(unique(meta_ko$experiment_id)), "\n")
cat("Unique species:", length(unique(meta_ko$organism)), "\n")

# Load libraries
library(dplyr)
# install.packages("metafor")
library(metafor)
library(data.table)


# Meta-analysis function
run_meta <- function(logFC, se, species, gene) {
  
  tryCatch({           # protect the script from crashing
    
    res <- rma(       # run meta-analysis (this fits a random-effects meta-analysis model)
      yi = logFC,     # effect sizes = LFC
      sei = se,       # standard error
      method = "REML" # tells how to estimate the between-study variance tau^2 that affects the weights used in the meta-analysis. Large τ² → studies differ a lot, Small τ² → studies are consistent. So if heterogeneity is large, no single study dominates.
    )
    
    data.table(       # extracting results 
      meta_LFC = as.numeric(res$b), # meta_LFC = Σ(w_i * logFC_i) / Σ(w_i), where w_i = 1 / (SE_i² + τ²)
      meta_SE = as.numeric(res$se),
      CI_lower = as.numeric(res$ci.lb), # 95% confidence interval: meta_LFC ± 1.96 * SE
      CI_upper = as.numeric(res$ci.ub), # If CI does not include 0 → significant effect
      pvalue = as.numeric(res$pval),
      tau2 = as.numeric(res$tau2), # between-study variance
      I2 = as.numeric(res$I2), # percentage of variability due to heterogeneity: low I² → conserved effect across species, high I² → species-specific effects
      k_experiments = as.numeric(res$k),
      n_species = as.numeric(uniqueN(species))
    )
    
  }, error = function(e) {
    
    data.table(
      meta_LFC = as.numeric(NA),
      meta_SE = as.numeric(NA),
      CI_lower = as.numeric(NA),
      CI_upper = as.numeric(NA),
      pvalue = as.numeric(NA),
      tau2 = as.numeric(NA),
      I2 = as.numeric(NA),
      k_experiments = as.numeric(length(logFC)),
      n_species = as.numeric(uniqueN(species))
    )
  })
}

# Convert meta_ko from dataframe to data table
setDT(meta_ko)

# Run meta-analysis per genes
meta_results <- meta_ko[
  ,
  run_meta(
    log2FoldChange,
    lfcSE,
    organism,
    unique(human_gene_id)
  ),
  by = human_gene_id
]

# Count the number of genes with NA meta_LFC (failed genes) -> 9 (small -> good)
sum(is.na(meta_results$meta_LFC))

# Show distribution of heterogeneity across genes
summary(meta_results$I2)

# Inspect Species Coverage
table(meta_results$n_species)

# Multiple testing correction because I tested ~15k genes
meta_results$FDR <- p.adjust(meta_results$pvalue, method = "BH") # Benjamini–Hochberg correction

# Write the table to csv file
write.csv(meta_results,
          "~/SIRT6_db/meta_analysis/meta_results_KO.csv",
          row.names = FALSE)







################################################# Visualizations #######################################################
meta_results <- read.csv("~/SIRT6_db/meta_analysis/meta_results_KO.csv")
summary(meta_results$meta_LFC)

######################################
## 1. Meta-analysis volcano plot #####
######################################

# Mark DEGs (up- and down-regulated)
meta_results$diffexpressed[meta_results$meta_LFC > 0.58 & meta_results$FDR < 0.05] <- "UP"
meta_results$diffexpressed[meta_results$meta_LFC < -0.58 & meta_results$FDR < 0.05] <- "DOWN"
# Fill NAs with "Insignificant"
meta_results$diffexpressed[is.na(meta_results$diffexpressed)] <- "Insignificant"

# save the table only with DEGs
diffgenes <- meta_results %>%
  filter(meta_results$diffexpressed != "Insignificant")
write.csv(diffgenes,
          "~/SIRT6_db/meta_analysis/meta_diffgenes_KO.csv",
          row.names = FALSE)

# Count DEGs
counts <- meta_results %>%
  count(diffexpressed)

n_up <- counts$n[counts$diffexpressed == "UP"]
n_down <- counts$n[counts$diffexpressed == "DOWN"]
n_all <- n_up + n_down

# Extract gene symbols for significant genes
# 1. upload orthologs table with gene symbol column
orthologs <- read.csv("~/SIRT6_db/meta_analysis/ortholog_map_1to1.csv")
# 2. Create mapping table
gene_map <- orthologs %>%
  select(human_gene_id, human_gene_symbol) %>%
  distinct()
# 3. Join with meta results
meta_results <- meta_results %>%
  left_join(gene_map, by = "human_gene_id")

# Select significant genes for labeling
library(ggrepel)
sig_genes <- meta_results %>%
  filter(diffexpressed != "Insignificant")

# Volcano plot
ggplot(data=meta_results, aes(x = meta_LFC, y = -log10(FDR))) +
  geom_vline(xintercept = c(-0.58, 0.58), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  
  geom_point(aes(color = diffexpressed)) +
  
  # add labels
  geom_text_repel(
    data = sig_genes,
    aes(label = human_gene_symbol),
    size = 3,
    max.overlaps = 20
  ) +
  
  scale_color_manual(
    values = c(
      "DOWN" = "steelblue",
      "UP" = "firebrick",
      "Insignificant" = "grey65"
    )
  ) +
  
  labs(
    x = "meta log2 Fold Change",
    y = "-log10 FDR",
    color = "DE genes",
    title = paste0("Meta-analysis of SIRT6 KO vs WT (", n_all, " DEGs: ", n_up, " UP, ", n_down, " DOWN)")
  )

# Count the number of significant genes (FDR < 0.05)
df <- meta_results %>%
  filter(meta_results$FDR < 0.05) # the number of rows is 605 => 605 significant genes

# Print top 5 significant genes in up and down groups 
top_5_lfc <- sig_genes %>%
  filter(diffexpressed %in% c("UP", "DOWN")) %>%
  group_by(diffexpressed) %>%
  # Sort by descending absolute LFC
  arrange(desc(abs(meta_LFC))) %>%
  slice_head(n = 5) %>%
  ungroup()


######################################################
###### 2. meta_LFC vs heterogeneity (I²) plot ########
######################################################

# Remove failed models
plot_data <- meta_results %>%
  filter(!is.na(meta_LFC))

# Label the most conserved genes
library(ggrepel)
conserved_genes <- plot_data %>%
  filter(I2 < 50 & diffexpressed != "Insignificant") %>%
  arrange(FDR) 

# Scatter plot
ggplot(plot_data, aes(x = meta_LFC, y = I2)) +
  
  geom_point(aes(color = diffexpressed), alpha = 0.7) +
  
  geom_hline(yintercept = 50, linetype = "dashed") +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed") +
  
  # label significant genes with I² < 50
  geom_text_repel(
    data = conserved_genes,
    aes(label = human_gene_symbol),
    size = 3
  ) +
  
  scale_color_manual(
    values = c(
      "DOWN" = "steelblue",
      "UP" = "firebrick",
      "Insignificant" = "grey65"
    )
  ) +
  
  labs(
    x = "Meta log2 Fold Change",
    y = "Heterogeneity (I²)",
    color = "DE genes",
    title = "Conserved SIRT6 transcriptional targets"
  ) +
  
  theme_minimal()


############################################
###### 3. Forest plot for key genes ########
############################################

library(metafor)
library(dplyr)

# Add human_gene_symbol column to meta_ko table
meta_ko <- meta_ko %>%
  left_join(gene_map, by = "human_gene_id")

# Add biological system column to meta_ko table
df <- read.csv("~/SIRT6_db/DE_results/Summary_table_and_circular_plot/summary_DE_table.csv")
# mapping table (exclude GSE102830 experiment id to prevent many-to-many mapping)
system_map <- df %>%
  select(experiment_id, biological_system) %>%
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

# Define genes for the forest plot
genes <- c(
  "H2AC7",
  "GAL", 
  "ADGRE2",
  "PKP1",
  "HSPB1",
  "CCL18"
)

# Create a function for constructing a forest plot
plot_forest_gene <- function(gene_name, data) {
  
  # filter data for a gene
  gene_data <- data %>%
    filter(human_gene_symbol == gene_name) %>%
    arrange(organism, log2FoldChange) # order by organism, then by LFC
  
  # run meta-analysis
  res <- rma(
    yi = gene_data$log2FoldChange,
    sei = gene_data$lfcSE,
    method = "REML"
  )
  
  # define labels for the plot
  labels <- paste0(
    gene_data$organism_short,
    " | ",
    gene_data$experiment_id,
    " | ",
    gene_data$biological_system
  )
  
  # compute the weight for each expreriment that shows how much each experiment contributes to the meta estimate
  # weights <- 1 / (gene_data$lfcSE^2 + res$tau2)
  # weights <- round(100 * weights / sum(weights), 1)
  
  # Add an extra space (bottom, left, top, right)
  par(mar = c(4, 4, 4, 2))
  
  # construct a forest plot with labels
  forest(
    res,
    slab = labels,
    mlab = "", # remove default diamond label
    xlab = "log2 Fold Change (SIRT6 KO vs WT)",
    main = gene_name,
    cex = 0.8,
    refline = 0, # reference line 
  )
  
  # add colored meta-analysis diamond
  addpoly(
    res,
    row = -1,
    mlab = "Meta-analysis", # meta-effect diamond label
    col = "firebrick",
    border = "firebrick"
  )
  }

# run the function for each of the genes from the list
for (g in genes) {
  plot_forest_gene(g, meta_ko)
}


