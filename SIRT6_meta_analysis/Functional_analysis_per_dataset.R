##################################################################################################
######################################## SIRT6 per-dataset GSEA ##################################
##################################################################################################

library(tidyverse)
library(fgsea)
library(msigdbr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(patchwork)

# Explicitly use dplyr::select to avoid conflicts
select <- dplyr::select

# Load data
de <- readRDS("~/SIRT6_db/DE_results/all_53_DE_results_in_one_table.rds")

metadata <- read.csv("~/SIRT6_db/DE_results/Summary_table_and_circular_plot/summary_DE_table.csv")

ortholog <- read.csv("~/SIRT6_db/meta_analysis/ortholog_map_1to1.csv")

metadata <- metadata %>%
  mutate(experiment_id = recode(experiment_id,
                                "GSE130690,2" = "GSE130690,GSE130692",
                                "GSE216185,6" = "GSE216185,GSE216186"
  ))

# Organism names
organism_name_map <- c(
  "drosophila_melanogaster" = "Drosophila melanogaster",
  "homo_sapiens"            = "Homo sapiens",
  "macaca_fascicularis"     = "Macaca fascicularis",
  "mus_musculus"            = "Mus musculus",
  "rattus_norvegicus"       = "Rattus norvegicus",
  "sus_scrofa"              = "Sus scrofa"
)

de <- de %>%
  mutate(organism_std = organism_name_map[organism])

# Contrast names
contrast_map <- c(
  "SIRT6.KO_vs_WT"     = "SIRT6 KO vs WT",
  "SIRT6.OE_vs_WT"     = "SIRT6 OE vs WT",
  "SIRT6.OE.K3R_vs_WT" = "SIRT6 OE K3R vs WT",
  "SIRT6.OE.K3Q_vs_WT" = "SIRT6 OE K3Q vs WT",
  "SIRT6.Het_vs_WT"    = "SIRT6 Het vs WT"
)

de <- de %>%
  mutate(contrast_std = contrast_map[contrast])

# Match the stratum
de <- de %>%
  mutate(stratum_std = stratum %>%
           str_replace("Adult \\(10 days\\)", "Adult 10 days") %>%
           str_replace("Aged \\(40 days\\)",  "Aged 40 days"))

#################################
# 1. Filter to relevant contrasts
#################################

ko_contrasts <- "SIRT6.KO_vs_WT"
oe_contrasts <- c("SIRT6.OE_vs_WT", "SIRT6.OE.K3R_vs_WT")
keep_contrasts <- c(ko_contrasts, oe_contrasts)

de_filtered <- de %>%
  filter(contrast %in% keep_contrasts) %>%
  mutate(contrast_group = if_else(contrast == "SIRT6.KO_vs_WT", "KO", "OE"))

###############################################
# 2. Build organism-specific hallmark gene sets 
###############################################

hallmarks_human <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  select(gs_name, gene_symbol)

# organism column in ortholog table uses "Genus species" format - matches organism_std
build_organism_genesets <- function(org_std, orthologs_df, hallmarks_df) {
  
  if (org_std == "Homo sapiens") {
    # Human: not in the `organism` column. Use human_gene_symbol -> human_gene_id
    # mapping which is present on every row of the ortholog table.
    org_mapping <- orthologs_df %>%
      select(target_gene_id    = human_gene_id,
                    target_gene_symbol = human_gene_symbol) %>%
      filter(!is.na(target_gene_id), !is.na(target_gene_symbol)) %>%
      distinct()
  } else {
    org_mapping <- orthologs_df %>%
      filter(organism == org_std) %>%
      select(target_gene_id    = organism_gene_id,
                    target_gene_symbol = human_gene_symbol) %>%
      filter(!is.na(target_gene_id), !is.na(target_gene_symbol)) %>%
      distinct()
  }
  
  mapped <- hallmarks_df %>%
    inner_join(org_mapping, by = c("gene_symbol" = "target_gene_symbol"),
               relationship = "many-to-many") %>%
    select(gs_name, target_gene_id) %>%
    distinct()
  
  gs_list <- split(mapped$target_gene_id, mapped$gs_name)
  gs_list <- gs_list[sapply(gs_list, length) >= 10]
  
  message(sprintf("  [%s] %d / 50 Hallmark gene sets retained after ortholog mapping",
                  org_std, length(gs_list)))
  return(gs_list)
}

message("Building organism-specific Hallmark gene sets...")

# organisms_in_data now uses the standardized "Genus species" names
organisms_in_data <- unique(de_filtered$organism_std)

organism_genesets <- map(
  set_names(organisms_in_data),
  ~ build_organism_genesets(.x, ortholog, hallmarks_human)
)

###############################################
# 3. Run GSEA per dataset, stratum and contrast
###############################################

run_gsea_unit <- function(df_unit, gs_list) {
  
  ranked <- df_unit %>%
    filter(!is.na(log2FoldChange), !is.na(pvalue), pvalue > 0) %>%
    mutate(pi_score = sign(log2FoldChange) * -log10(pvalue)) %>%
    group_by(gene_id) %>%
    summarise(pi_score = mean(pi_score), .groups = "drop") %>%
    arrange(desc(pi_score))
  
  ranks <- setNames(ranked$pi_score, ranked$gene_id)
  if (length(ranks) < 100) return(NULL)
  
  set.seed(42)
  result <- tryCatch(
    fgsea(pathways = gs_list, stats = ranks,
          minSize = 10, maxSize = 500, nPermSimple = 1000),
    error = function(e) { message("    fgsea error: ", e$message); NULL }
  )
  return(result)
}

message("\nRunning GSEA for all dataset units...")

# Use harmonized columns: organism_std and stratum_std
units <- de_filtered %>%
  distinct(experiment_id, stratum_std, contrast, contrast_group, organism_std)

gsea_results_list <- vector("list", nrow(units))

for (i in seq_len(nrow(units))) {
  
  unit <- units[i, ]
  org  <- unit$organism_std   # "Genus species" format — matches organism_genesets names
  
  if (!org %in% names(organism_genesets)) {
    message(sprintf("  [%d/%d] Skipping %s | %s | %s — no gene sets for %s",
                    i, nrow(units), unit$experiment_id,
                    unit$stratum_std, unit$contrast, org))
    next
  }
  
  df_unit <- de_filtered %>%
    filter(experiment_id == unit$experiment_id,
           stratum_std   == unit$stratum_std,
           contrast      == unit$contrast)
  
  message(sprintf("  [%d/%d] %s | %s | %s | %s (%d genes)",
                  i, nrow(units), org, unit$experiment_id,
                  unit$stratum_std, unit$contrast, nrow(df_unit)))
  
  res <- run_gsea_unit(df_unit, organism_genesets[[org]])
  
  if (!is.null(res)) {
    gsea_results_list[[i]] <- res %>%
      mutate(
        experiment_id  = unit$experiment_id,
        stratum        = unit$stratum_std,   # store harmonized stratum
        contrast       = unit$contrast,
        contrast_group = unit$contrast_group,
        organism       = org
      )
  }
}

gsea_all <- bind_rows(compact(gsea_results_list))
message(sprintf("\nGSEA complete: %d pathway × unit results", nrow(gsea_all)))

library(readr)
write_csv(gsea_all,
          "~/SIRT6_db/meta_analysis/gsea_per_dataset.csv"
          )

################################################
# 4. Join with metadata and prepare for plotting
################################################

# Reverse contrast map so we can translate metadata contrasts back to de format
# for building a shared unit_id key
contrast_map_rev <- setNames(names(contrast_map), unname(contrast_map))

meta_clean <- metadata %>%
  filter(contrast %in% c("SIRT6 KO vs WT", "SIRT6 OE vs WT", "SIRT6 OE K3R vs WT")) %>%
  mutate(
    contrast_group = if_else(contrast == "SIRT6 KO vs WT", "KO", "OE"),
    # Translate metadata contrast to de-style contrast for unit_id alignment
    contrast_de    = contrast_map_rev[contrast],
    unit_id        = paste(experiment_id, stratum, contrast_de, sep = " | ")
  )

# unit_id in gsea_all already uses de-style contrast and harmonized stratum
gsea_all <- gsea_all %>%
  mutate(unit_id = paste(experiment_id, stratum, contrast, sep = " | "))

# Redo the metadata join (drops organism & contrast_group from meta_clean to avoid .x/.y)
gsea_plot <- gsea_all %>%
  left_join(
    meta_clean %>% select(unit_id, biological_system, biological_context,
                                 global_order, n_deg),
    by = "unit_id"
  ) %>%
  mutate(pathway_clean = str_replace(pathway, "HALLMARK_", "") %>%
           str_replace_all("_", " ") %>%
           str_to_title())

colnames(gsea_plot)                        # should show "organism" and "contrast_group" (no .x/.y)

# Sanity check: how many GSEA results failed to join?
n_missing <- sum(is.na(gsea_plot$biological_system))
if (n_missing > 0) {
  message(sprintf("Warning: %d GSEA rows did not match any metadata entry", n_missing))
  # Inspect mismatches:
  # gsea_plot %>% filter(is.na(biological_system)) %>%
  #   distinct(experiment_id, stratum, contrast) %>% print(n = 30)
}

# Should be 0
sum(is.na(gsea_plot$biological_system))





#####################################
# 5. Define display order and filters
#####################################

# 5a. Order datasets by: contrast_group (KO first, OE second),then organism, then global_order from metadata
unit_order <- gsea_plot %>%
  distinct(unit_id, contrast_group, organism, global_order, biological_system) %>%
  arrange(
    contrast_group,
    factor(organism, levels = c("Homo sapiens", "Mus musculus",
                                "Rattus norvegicus", "Macaca fascicularis",
                                "Sus scrofa", "Drosophila melanogaster")),
    global_order
  ) %>%
  mutate(unit_rank = row_number())

gsea_plot <- gsea_plot %>%
  left_join(unit_order %>% select(unit_id, unit_rank), by = "unit_id")

# 5b. Select pathways to display:
# Keep pathways significant (padj < 0.05) in at least 10 units
gsea_plot <- as_tibble(gsea_plot)

sig_pathways <- gsea_plot %>%
  filter(padj < 0.05) %>%
  group_by(pathway) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n >= 10) %>%
  pull(pathway)

length(sig_pathways)   

message(sprintf("Pathways significant in >= 10 units: %d", length(sig_pathways)))

# Order pathways: cluster by NES pattern (use hierarchical clustering)
nes_matrix_for_clustering <- gsea_plot %>%
  filter(pathway %in% sig_pathways) %>%
  select(pathway, unit_id, NES) %>%
  pivot_wider(names_from = unit_id, values_from = NES, values_fill = 0) %>%
  column_to_rownames("pathway") %>%
  as.matrix()

pathway_clust  <- hclust(dist(nes_matrix_for_clustering), method = "ward.D2")
pathway_order  <- rownames(nes_matrix_for_clustering)[pathway_clust$order]



###########################################
# 6. Build NES matrix for heatmap (panel A)
###########################################

# Show all NES values regardless of significance 
# Overlay asterisks for significant cells (* padj<0.05, ** padj<0.01, *** padj<0.001).

# Pathway label cleanup — proper biology formatting
pathway_label_map <- c(
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"            = "TNFA signaling via NFKB",
  "HALLMARK_HYPOXIA"                            = "Hypoxia",
  "HALLMARK_CHOLESTEROL_HOMEOSTASIS"            = "Cholesterol homeostasis",
  "HALLMARK_MITOTIC_SPINDLE"                    = "Mitotic spindle",
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING"         = "Wnt/β-catenin signaling",
  "HALLMARK_TGF_BETA_SIGNALING"                 = "TGF-β signaling",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING"            = "IL6/JAK/STAT3 signaling",
  "HALLMARK_DNA_REPAIR"                         = "DNA repair",
  "HALLMARK_G2M_CHECKPOINT"                     = "G2/M checkpoint",
  "HALLMARK_APOPTOSIS"                          = "Apoptosis",
  "HALLMARK_NOTCH_SIGNALING"                    = "Notch signaling",
  "HALLMARK_ADIPOGENESIS"                       = "Adipogenesis",
  "HALLMARK_ESTROGEN_RESPONSE_EARLY"            = "Estrogen response (early)",
  "HALLMARK_ESTROGEN_RESPONSE_LATE"             = "Estrogen response (late)",
  "HALLMARK_ANDROGEN_RESPONSE"                  = "Androgen response",
  "HALLMARK_MYOGENESIS"                         = "Myogenesis",
  "HALLMARK_SPERMATOGENESIS"                    = "Spermatogenesis",
  "HALLMARK_PANCREAS_BETA_CELLS"                = "Pancreas β-cells",
  "HALLMARK_INFLAMMATORY_RESPONSE"              = "Inflammatory response",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE"          = "IFN-alpha response",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE"          = "IFN-gamma response",
  "HALLMARK_APICAL_JUNCTION"                    = "Apical junction",
  "HALLMARK_APICAL_SURFACE"                     = "Apical surface",
  "HALLMARK_HEDGEHOG_SIGNALING"                 = "Hedgehog signaling",
  "HALLMARK_COMPLEMENT"                         = "Complement",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"          = "Unfolded protein response",
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING"            = "PI3K/AKT/mTOR signaling",
  "HALLMARK_MTORC1_SIGNALING"                   = "mTORC1 signaling",
  "HALLMARK_E2F_TARGETS"                        = "E2F targets",
  "HALLMARK_MYC_TARGETS_V1"                     = "MYC targets V1",
  "HALLMARK_MYC_TARGETS_V2"                     = "MYC targets V2",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"  = "Epithelial-mesenchymal transition",
  "HALLMARK_INFLAMMATORY_RESPONSE"              = "Inflammatory response",
  "HALLMARK_XENOBIOTIC_METABOLISM"              = "Xenobiotic metabolism",
  "HALLMARK_FATTY_ACID_METABOLISM"              = "Fatty acid metabolism",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION"          = "Oxidative phosphorylation",
  "HALLMARK_GLYCOLYSIS"                         = "Glycolysis",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"    = "Reactive oxygen species",
  "HALLMARK_P53_PATHWAY"                        = "p53 pathway",
  "HALLMARK_UV_RESPONSE_UP"                     = "UV response (up)",
  "HALLMARK_UV_RESPONSE_DN"                     = "UV response (down)",
  "HALLMARK_ANGIOGENESIS"                       = "Angiogenesis",
  "HALLMARK_HEME_METABOLISM"                    = "Heme metabolism",
  "HALLMARK_COAGULATION"                        = "Coagulation",
  "HALLMARK_IL2_STAT5_SIGNALING"                = "IL2/STAT5 signaling",
  "HALLMARK_BILE_ACID_METABOLISM"               = "Bile acid metabolism",
  "HALLMARK_PEROXISOME"                         = "Peroxisome",
  "HALLMARK_ALLOGRAFT_REJECTION"                = "Allograft rejection",
  "HALLMARK_KRAS_SIGNALING_UP"                  = "KRAS signaling (up)",
  "HALLMARK_KRAS_SIGNALING_DN"                  = "KRAS signaling (down)",
  "HALLMARK_PROTEIN_SECRETION"                  = "Protein secretion"
)

# Helper: clean any HALLMARK_* name; fall back to title case if not in map
clean_pathway <- function(x) {
  ifelse(x %in% names(pathway_label_map),
         pathway_label_map[x],
         str_replace(x, "HALLMARK_", "") %>% str_replace_all("_", " ") %>% str_to_title())
}

plot_data_A <- gsea_plot %>%
  filter(pathway %in% sig_pathways) %>%
  mutate(pathway_clean = clean_pathway(pathway))

# Build wide NES matrix - all values (no NA filtering by significance)
nes_mat <- plot_data_A %>%
  dplyr::select(pathway_clean, unit_id, NES) %>%
  pivot_wider(names_from = unit_id, values_from = NES) %>%
  column_to_rownames("pathway_clean") %>%
  as.matrix()

# Build parallel padj matrix for asterisk overlay
padj_mat <- plot_data_A %>%
  dplyr::select(pathway_clean, unit_id, padj) %>%
  pivot_wider(names_from = unit_id, values_from = padj) %>%
  column_to_rownames("pathway_clean") %>%
  as.matrix()

# Re-order rows by clustering, columns by unit_rank
pathway_order_clean <- clean_pathway(pathway_order)

col_order <- unit_order %>% arrange(unit_rank) %>% pull(unit_id)
col_order <- col_order[col_order %in% colnames(nes_mat)]

nes_mat  <- nes_mat[pathway_order_clean[pathway_order_clean %in% rownames(nes_mat)],
                    col_order, drop = FALSE]
padj_mat <- padj_mat[rownames(nes_mat), colnames(nes_mat), drop = FALSE]


################################
# 7. Annotation data for panel A
################################

col_anno_df <- unit_order %>%
  filter(unit_id %in% col_order) %>%
  arrange(unit_rank) %>%
  select(unit_id, organism, biological_system, contrast_group) %>%
  column_to_rownames("unit_id")

# Color palettes
organism_colors <- c(
  "Homo sapiens"          = "#E41A1C",
  "Mus musculus"          = "#377EB8",
  "Rattus norvegicus"     = "#4DAF4A",
  "Macaca fascicularis"   = "#FF7F00",
  "Sus scrofa"            = "#984EA3",
  "Drosophila melanogaster" = "#A65628"
)

# 10 biological systems - assign a color palette
bio_sys_vals <- unique(col_anno_df$biological_system)
bio_sys_colors <- setNames(
  colorRampPalette(RColorBrewer::brewer.pal(10, "Set3"))(length(bio_sys_vals)),
  bio_sys_vals
)

contrast_colors <- c("KO" = "palevioletred", "OE" = "chartreuse4") # KO - "#D95F02", OE - "#1B9E77"

col_anno <- HeatmapAnnotation(
  `Contrast`          = col_anno_df[col_order, "contrast_group"],
  `Organism`          = col_anno_df[col_order, "organism"],
  `Biological system` = col_anno_df[col_order, "biological_system"],
  col = list(
    `Contrast`          = contrast_colors,
    `Organism`          = organism_colors,
    `Biological system` = bio_sys_colors
  ),
  annotation_height    = unit(c(4, 4, 4), "mm"),
  show_legend          = TRUE,
  annotation_name_side = "left",
  annotation_name_gp   = gpar(fontsize = 14, fontface = "bold"),
  annotation_legend_param = list(
    `Contrast`          = list(title_gp = gpar(fontsize = 16, fontface = "bold"),
                               labels_gp = gpar(fontsize = 14)),
    `Organism`          = list(title_gp = gpar(fontsize = 16, fontface = "bold"),
                               labels_gp = gpar(fontsize = 14)),
    `Biological system` = list(title_gp = gpar(fontsize = 16, fontface = "bold"),
                               labels_gp = gpar(fontsize = 14))
  )
)

# NES color scale: diverging blue-white-red
col_fun <- colorRamp2(
  c(-3, -1.5, 0, 1.5, 3),
  c("#2166AC", "#92C5DE", "white", "#F4A582", "#B2182B")
)

##################################
# 8. Panel A — Per-dataset heatmap
##################################

# Asterisk overlay function: marks significance level per cell
sig_to_stars <- function(p) {
  if (is.na(p))         return("")
  if (p < 0.0005)        return("***")
  if (p < 0.005)         return("**")
  if (p < 0.05)         return("*")
  return("")
}

cell_fun_stars <- function(j, i, x, y, width, height, fill) {
  star <- sig_to_stars(padj_mat[i, j])
  if (nzchar(star)) {
    grid.text(star, x, y,
              gp = gpar(fontsize = 12, col = "black"),
              vjust = 0.75)
  }
}

hm_A <- Heatmap(
  nes_mat,
  name             = "NES",
  col              = col_fun,
  na_col           = "grey92",
  cluster_rows     = FALSE,
  cluster_columns  = FALSE,
  show_column_names = FALSE,
  row_names_side   = "left",
  row_names_gp     = gpar(fontsize = 14),
  column_split     = col_anno_df[col_order, "contrast_group"],
  column_gap       = unit(3, "mm"),
  top_annotation   = col_anno,
  cell_fun         = cell_fun_stars,
  border           = TRUE,
  width            = unit(ncol(nes_mat) * 6, "mm"),
  height           = unit(nrow(nes_mat) * 9, "mm"),
  heatmap_legend_param = list(
    title          = "NES",
    at             = c(-3, -1.5, 0, 1.5, 3),
    labels         = c("-3", "-1.5", "0", "1.5", "3"),
    legend_height  = unit(45, "mm"),
    title_gp       = gpar(fontsize = 16, fontface = "bold"),
    labels_gp      = gpar(fontsize = 14)
  ),
  column_title_gp  = gpar(fontsize = 14, fontface = "bold"),
  row_title_gp     = gpar(fontsize = 8)
)

###############################################################
# 9. Panel B - per-dataset heatmap grouped by biological system
###############################################################

# The same data as Panel A, but columns ordered/split by biological_system instead of contrast group

# 9a. New unit ordering: by biological_system, then contrast_group, then organism
unit_order_B <- gsea_plot %>%
  distinct(unit_id, contrast_group, organism, biological_system, global_order) %>%
  arrange(
    biological_system,
    contrast_group,
    factor(organism, levels = c("Homo sapiens", "Mus musculus",
                                "Rattus norvegicus", "Macaca fascicularis",
                                "Sus scrofa", "Drosophila melanogaster")),
    global_order
  ) %>%
  mutate(unit_rank_B = row_number())

col_order_B <- unit_order_B %>% arrange(unit_rank_B) %>% pull(unit_id)
col_order_B <- col_order_B[col_order_B %in% colnames(nes_mat)]

# Reorder both matrices to new column order
nes_mat_B  <- nes_mat[,  col_order_B, drop = FALSE]
padj_mat_B <- padj_mat[, col_order_B, drop = FALSE]

# 9b. Annotation data frame in new order
col_anno_df_B <- unit_order_B %>%
  filter(unit_id %in% col_order_B) %>%
  arrange(unit_rank_B) %>%
  select(unit_id, organism, biological_system, contrast_group) %>%
  column_to_rownames("unit_id")

# 9c. Annotation: now Contrast is informative within each biological system block
col_anno_B <- HeatmapAnnotation(
  `Contrast` = col_anno_df_B[col_order_B, "contrast_group"],
  `Organism` = col_anno_df_B[col_order_B, "organism"],
  col = list(
    `Contrast` = contrast_colors,
    `Organism` = organism_colors
  ),
  annotation_height    = unit(c(4, 4), "mm"),
  show_legend          = c(TRUE, TRUE),
  annotation_name_side = "left",
  annotation_name_gp   = gpar(fontsize = 14, fontface = "bold"),
  annotation_legend_param = list(
    `Contrast` = list(title_gp  = gpar(fontsize = 16, fontface = "bold"),
                      labels_gp = gpar(fontsize = 14)),
    `Organism` = list(title_gp  = gpar(fontsize = 16, fontface = "bold"),
                      labels_gp = gpar(fontsize = 14))
  )
)

# 9d. Asterisk overlay using padj_mat_B (different column order from Panel A)
cell_fun_stars_B <- function(j, i, x, y, width, height, fill) {
  star <- sig_to_stars(padj_mat_B[i, j])
  if (nzchar(star)) {
    grid.text(star, x, y,
              gp = gpar(fontsize = 12, col = "black"),
              vjust = 0.75)
  }
}

# 9e. Build Panel B heatmap, split columns by biological_system
hm_B <- Heatmap(
  nes_mat_B,
  name              = "NES",
  col               = col_fun,
  na_col            = "grey92",
  cluster_rows      = FALSE,
  cluster_columns   = FALSE,
  show_column_names = FALSE,
  row_names_side    = "left",
  row_names_gp      = gpar(fontsize = 14),
  column_split      = col_anno_df_B[col_order_B, "biological_system"],
  column_gap        = unit(1.5, "mm"),
  column_title_rot  = 45,
  column_title_gp   = gpar(fontsize = 14, fontface = "bold"),
  top_annotation    = col_anno_B,
  cell_fun          = cell_fun_stars_B,
  border            = TRUE,
  width             = unit(ncol(nes_mat_B) * 6, "mm"),
  height            = unit(nrow(nes_mat_B) * 9, "mm"),
  heatmap_legend_param = list(
    title         = "NES",
    at            = c(-3, -1.5, 0, 1.5, 3),
    labels        = c("-3", "-1.5", "0", "1.5", "3"),
    legend_height = unit(45, "mm"),
    title_gp      = gpar(fontsize = 16, fontface = "bold"),
    labels_gp     = gpar(fontsize = 14)
  )
)

##################################
# 10. Save panels as separate PDFs
##################################

out_dir <- "~/SIRT6_db/meta_analysis"

# --- Panel A ---
pdf(file.path(out_dir, "sirt6_functional_properties_panelA.pdf"),
    width = 20, height = 14)

draw(
  hm_A,
  padding = unit(c(5, 5, 5, 30), "mm"),
  heatmap_legend_side    = "right",
  annotation_legend_side = "right",
  column_title    = "A  Functional properties of SIRT6 KO and OE (by contrast & organism)",
  column_title_gp = gpar(fontsize = 18, fontface = "bold")
)

dev.off()

# --- Panel B ---
pdf(file.path(out_dir, "sirt6_functional_properties_panelB.pdf"),
    width = 20, height = 14)

draw(
  hm_B,
  padding = unit(c(5, 5, 5, 30), "mm"),
  heatmap_legend_side    = "right",
  annotation_legend_side = "right",
  column_title    = "B  Functional properties of SIRT6 KO and OE (by biological system)",
  column_title_gp = gpar(fontsize = 18, fontface = "bold")
)

dev.off()

message("Saved Panel A and Panel B PDFs to: ", out_dir)

#################################################
# 11. Save GSEA results table (for Supplementary)
#################################################

gsea_plot %>%
  select(contrast_group, organism, biological_system, biological_context,
                experiment_id, stratum, contrast, pathway, pathway_clean,
                NES, pval, padj, size, leadingEdge) %>%
  mutate(leadingEdge = map_chr(leadingEdge, paste, collapse = ";")) %>%
  write_csv("~/SIRT6_db/meta_analysis/sirt6_gsea_all_results.csv")

