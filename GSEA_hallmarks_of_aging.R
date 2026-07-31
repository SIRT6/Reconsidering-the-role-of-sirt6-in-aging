graphics.off()
rm(list = ls())

library(fs)
library(nanoparquet)
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
library(purrr)
library(fgsea)
library(org.Hs.eg.db)    
library(AnnotationDbi)
library(dplyr)

select <- dplyr::select

# ---- Ortholog table (1:1, all species -> human) ----
ortho <- read.csv('/tank/projects/public_data/Sirt6_datasets/Expression/DE_results/ortholog_map_1to1.csv')
cat("Organisms in table:\n")
print(unique(ortho$organism))
head(ortho)
tail(ortho)

# ---- Generic OpenGenes tsv reader ----
read_opengenes_tsv <- function(path) {
  df <- read.delim(path, sep = "\t", header = TRUE, quote = "",
                   stringsAsFactors = FALSE, check.names = FALSE)
  colnames(df) <- gsub('^"|"$|^X\\.|\\.$', '', colnames(df))
  df %>% mutate(across(everything(), ~ gsub('^"|"$', '', .x)))
}

# ---- Gene confidence level (human DB) ----
genes_conf <- read_opengenes_tsv("gene-confidence-level.tsv") %>%
  dplyr::select(hgnc, confidence_level)

score_map <- c("highest" = 5, "high" = 4, "moderate" = 3, "low" = 2)

# ---- Hallmarks of aging (human DB, hgnc) ----
hallmarks_raw <- read_opengenes_tsv("gene-aging-mechanisms.tsv") %>%
  dplyr::select(hgnc, hallmarks_of_aging) %>%
  mutate(hallmarks_of_aging = str_trim(hallmarks_of_aging)) %>%
  separate_rows(hallmarks_of_aging, sep = ",") %>%
  mutate(hallmarks_of_aging = str_trim(hallmarks_of_aging)) %>%
  filter(hallmarks_of_aging != "" & !is.na(hallmarks_of_aging))

# ---- Aging score by confidence level (direct join on human_gene_symbol) ----
compute_aging_score_by_confidence <- function(df_sig, genes_conf) {
  if (nrow(df_sig) == 0) {
    return(tibble(total_score = 0, n_genes = 0, n_orth_found = 0, n_scored = 0))
  }
  
  result <- df_sig %>%
    left_join(
      genes_conf %>% mutate(confidence_level = tolower(confidence_level)),
      by = c("human_gene_symbol" = "hgnc")
    ) %>%
    mutate(
      gene_score = ifelse(
        !is.na(confidence_level) & confidence_level %in% names(score_map),
        score_map[confidence_level],
        0
      )
    )
  
  tibble(
    total_score  = sum(result$gene_score),
    n_genes      = nrow(result),
    n_orth_found = sum(!is.na(result$human_gene_symbol)),
    n_scored     = sum(result$gene_score > 0)
  )
}

# ---- Build geneSets for a given organism ----
# hallmarks_raw$hgnc (human) -> organism_gene_symbol via ortho
build_geneSets_for_organism <- function(hallmarks_raw, ortho, organism_name) {
  ortho_sp <- ortho %>% filter(organism == organism_name)
  
  hallmarks_org <- hallmarks_raw %>%
    inner_join(ortho_sp %>% select(human_gene_symbol, organism_gene_symbol),
               by = c("hgnc" = "human_gene_symbol"))
  
  hallmarks_org %>%
    distinct(organism_gene_symbol, hallmarks_of_aging) %>%
    group_by(hallmarks_of_aging) %>%
    summarize(genes = list(unique(organism_gene_symbol))) %>%
    deframe()
}

# ---- Process dataset (shared across organisms) ----
process_dataset <- function(parquet_file, ortho_sp, genes_conf, geneSets, out_dir) {
  
  dataset_id      <- basename(dirname(parquet_file))
  dataset_id_safe <- gsub("[^A-Za-z0-9_-]", "_", dataset_id)
  
  cat("\n============================\n")
  cat("Processing:", dataset_id, "\n")
  
  df <- read_parquet(parquet_file)
  if (!"gene_id" %in% colnames(df)) {
    warning(dataset_id, ": no 'gene_id' column, skipping")
    return(NULL)
  }
  
  # ---- Annotate + map to human nomenclature in one join ----
  df2 <- df %>%
    left_join(
      ortho_sp %>% select(organism_gene_id, organism_gene_symbol, human_gene_id, human_gene_symbol),
      by = c("gene_id" = "organism_gene_id")
    )
  
  n_dup_check <- nrow(df2) - length(unique(df2$gene_id))
  if (n_dup_check > 0) {
    warning(dataset_id, ": ", n_dup_check, " extra rows after annotation (check ortho duplicates)")
  }
  
  sig <- df2 %>% filter(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 0.58)
  cat("  Significant genes:", nrow(sig), "| with ortholog found:",
      sum(!is.na(sig$human_gene_symbol)), "\n")
  
  if (nrow(sig) == 0) {
    warning(dataset_id, ": no significant genes, skipping")
    return(NULL)
  }
  
  # ---- Aging score ----
  aging_score_conf <- compute_aging_score_by_confidence(sig, genes_conf)
  aging_score_conf$dataset <- dataset_id
  cat("  Aging score:", aging_score_conf$total_score,
      "| orthologs found:", aging_score_conf$n_orth_found, "/", aging_score_conf$n_genes,
      "| scored:", aging_score_conf$n_scored, "\n")
  
  # ---- fgsea: stats keyed by organism_gene_symbol (matches geneSets) ----
  if ("stat" %in% colnames(sig)) {
    stats_vec <- sig$stat
    names(stats_vec) <- sig$organism_gene_symbol
  } else if ("log2FoldChange" %in% colnames(sig)) {
    stats_vec <- sig$log2FoldChange
    names(stats_vec) <- sig$organism_gene_symbol
  } else {
    warning(dataset_id, ": no 'stat'/'log2FoldChange', skipping fgsea")
    return(list(aging_score = aging_score_conf, fgsea = NULL))
  }
  
  stats_vec <- stats_vec[!is.na(names(stats_vec)) & !is.na(stats_vec) & names(stats_vec) != ""]
  
  if (any(duplicated(names(stats_vec)))) {
    stats_df <- tibble(gene = names(stats_vec), stat = stats_vec) %>%
      group_by(gene) %>%
      summarize(stat = stat[which.max(abs(stat))], .groups = "drop")
    stats_vec <- stats_df$stat
    names(stats_vec) <- stats_df$gene
  }
  stats_vec <- sort(stats_vec, decreasing = TRUE)
  
  if (length(stats_vec) < 5) {
    warning(dataset_id, ": too few genes for fgsea (<5), skipping")
    return(list(aging_score = aging_score_conf, fgsea = NULL))
  }
  
  set.seed(42)
  fgsea_res <- fgseaMultilevel(pathways = geneSets, stats = stats_vec, minSize = 5, maxSize = 5000)
  fgsea_res_sig <- subset(fgsea_res, padj < 0.05)
  
  if (nrow(fgsea_res_sig) == 0) {
    cat("  ", dataset_id, ": no significant pathways\n")
    return(list(aging_score = aging_score_conf, fgsea = NULL))
  }
  
  fgsea_res_sig_fixed <- fgsea_res_sig %>%
    mutate(
      leadingEdge = map_chr(leadingEdge, ~ paste(.x, collapse = ";")),
      dataset     = dataset_id
    )
  
  out_file <- file.path(out_dir, paste0(dataset_id_safe, "_all_SIRT6-KO_vs_WT.csv"))
  write.csv(fgsea_res_sig_fixed, out_file, row.names = FALSE)
  cat("  Saved:", out_file, "\n")
  
  list(aging_score = aging_score_conf, fgsea = fgsea_res_sig_fixed)
}

# ---- Batch runner for one organism ----
run_organism_batch <- function(species_dir_name, organism_name, ortho, genes_conf, hallmarks_raw) {
  
  de_dir  <- paste0("/tank/projects/public_data/Sirt6_datasets/Expression/DE_results/", species_dir_name, "/results/")
  out_dir <- file.path("/tank/projects/ekashuk/sirt6_biomarker/expression", species_dir_name)
  dir_create(out_dir)
  
  all_parquet <- dir_ls(de_dir, recurse = TRUE, glob = "*.parquet")
  cat("=== ", organism_name, " ===\n")
  cat("Total files:", length(all_parquet), "\n")
  print(basename(unique(dirname(all_parquet))))
  
  # ---- Special case: Homo sapiens — ortho not needed 
  if (organism_name == "Homo sapiens") {
    
    all_ids <- map(all_parquet, ~ read_parquet(.x)$gene_id) %>% unlist() %>% unique() %>% na.omit()
    
    ann <- AnnotationDbi::select(
      org.Hs.eg.db, keys = all_ids, keytype = "ENSEMBL",
      columns = c("SYMBOL", "ENTREZID")
    ) %>%
      distinct(ENSEMBL, .keep_all = TRUE)
    
    cat("Genes annotated:", nrow(ann), "of", length(all_ids), "\n")
    
    ortho_sp <- ann %>%
      transmute(
        organism_gene_id     = ENSEMBL,
        organism_gene_symbol = SYMBOL,
        human_gene_id        = ENSEMBL,
        human_gene_symbol    = SYMBOL
      )
  } else {
    ortho_sp <- ortho %>% filter(organism == organism_name)
  }
  
  cat("Orthologs/annotations for", organism_name, ":", nrow(ortho_sp), "\n")
  
  geneSets <- if (organism_name == "Homo sapiens") {
    hallmarks_raw %>%
      distinct(hgnc, hallmarks_of_aging) %>%
      group_by(hallmarks_of_aging) %>%
      summarize(genes = list(unique(hgnc))) %>%
      deframe()
  } else {
    build_geneSets_for_organism(hallmarks_raw, ortho, organism_name)
  }
  cat("Geneset sizes:\n")
  print(lengths(geneSets))
  
  results <- map(all_parquet, function(f) {
    tryCatch(
      process_dataset(f, ortho_sp = ortho_sp, genes_conf = genes_conf, geneSets = geneSets, out_dir = out_dir),
      error = function(e) { warning("Failed on ", f, ": ", conditionMessage(e)); NULL }
    )
  })
  names(results) <- basename(dirname(all_parquet))
  
  aging_scores_all <- map_dfr(results, ~ .x$aging_score) %>%
    mutate(orth_match_rate = round(n_orth_found / n_genes, 3)) %>%
    relocate(dataset, .before = everything())
  print(aging_scores_all)
  write.csv(aging_scores_all, file.path(out_dir, "aging_scores_all_datasets.csv"), row.names = FALSE)
  
  fgsea_all <- map_dfr(results, ~ .x$fgsea)
  write.csv(fgsea_all, file.path(out_dir, "fgsea_results_all_datasets.csv"), row.names = FALSE)
  
  cat("Done. Processed", length(all_parquet), "files,",
      sum(!map_lgl(results, is.null)), "successful.\n\n")
  
  list(aging_scores = aging_scores_all, fgsea = fgsea_all)
}

## Homo sapiens

res_human <- run_organism_batch(
  species_dir_name = "homo_sapiens",
  organism_name    = "Homo sapiens",
  ortho            = ortho,
  genes_conf       = genes_conf,
  hallmarks_raw    = hallmarks_raw
)

## Mus musculus ----
res_mouse <- run_organism_batch(
  species_dir_name = "mus_musculus",
  organism_name    = "Mus musculus",   
  ortho            = ortho,
  genes_conf       = genes_conf,
  hallmarks_raw    = hallmarks_raw
)

## Macaca fascicularis ----
res_macaque <- run_organism_batch(
  species_dir_name = "macaca_fascicularis",
  organism_name    = "Macaca fascicularis",   
  ortho            = ortho,
  genes_conf       = genes_conf,
  hallmarks_raw    = hallmarks_raw
)

## Rattus norvegicus ----

res_rat <- run_organism_batch(
  species_dir_name = "rattus_norvegicus",
  organism_name    = "Rattus norvegicus",   
  ortho            = ortho,
  genes_conf       = genes_conf,
  hallmarks_raw    = hallmarks_raw
)

## Drosophila melanogaster ----

res_fly <- run_organism_batch(
  species_dir_name = "drosophila_melanogaster",
  organism_name    = "Drosophila melanogaster",   
  ortho            = ortho,
  genes_conf       = genes_conf,
  hallmarks_raw    = hallmarks_raw
)

## Sus scrofa ---- 
res_pig <- run_organism_batch(
  species_dir_name = "sus_scrofa",
  organism_name    = "Sus scrofa",   
  ortho            = ortho,
  genes_conf       = genes_conf,
  hallmarks_raw    = hallmarks_raw
)

## Bubble plot ----

summary_df <- read.csv('/tank/projects/public_data/Sirt6_datasets/Expression/DE_results/summary_DE_table.csv', row.names = 1)

# ============================================================
# 0. Color palettes
# ============================================================

contrast_colors <- c(
  "SIRT6 KO vs WT"      = "palevioletred",
  "SIRT6 Het vs WT"     = "lightcoral",
  "SIRT6 OE vs WT"      = "chartreuse4",
  "SIRT6 OE K3Q vs WT"  = "chartreuse2",
  "SIRT6 OE K3R vs WT"  = "darkolivegreen1"
)

system_colors <- c(
  "Nervous"          = "#66C5CCFF",
  "Cardiovascular"   = "#F6CF71FF",
  "Respiratory"      = "#F89C74FF",
  "Metabolic"        = "#B497E7FF",
  "Gastrointestinal" = "#87C55FFF",
  "Musculoskeletal"  = "#9EB9F3FF",
  "Immune"           = "#FE88B1FF",
  "Developmental"    = "#B3B3B3FF",
  "Cancer"           = "#8BE0A4FF",
  "Unknown"          = "grey70"
)

# ============================================================
# 2. Load fgsea results across all organisms
# ============================================================

base_path <- "/tank/projects/ekashuk/sirt6_biomarker/expression"
organism_dirs <- c(
  "homo_sapiens", "macaca_fascicularis", "mus_musculus",
  "rattus_norvegicus", "sus_scrofa", "drosophila_melanogaster"
)

all_data <- map_dfr(organism_dirs, function(org){
  
  dir_path <- file.path(base_path, org)
  if (!dir.exists(dir_path)) {
    warning("Directory not found: ", dir_path)
    return(NULL)
  }
  
  files <- list.files(dir_path, pattern = "*.csv", full.names = TRUE)
  files <- files[!grepl("all_datasets", basename(files))]  # exclude summary files
  
  map_dfr(files, function(f){
    df <- read_csv(f, show_col_types = FALSE)
    if (!all(c("padj", "pathway", "NES") %in% colnames(df))) {
      warning("Skipping file (missing required columns): ", basename(f))
      return(NULL)
    }
    df %>%
      mutate(
        organism   = org,
        experiment = tools::file_path_sans_ext(basename(f)),
        logFDR     = -log10(padj),
        pathway    = str_remove_all(pathway, "^'|'$")
      )
  })
})


all_data <- all_data %>%
  mutate(
    organism = recode(organism,
                      homo_sapiens = "Homo sapiens",
                      macaca_fascicularis = "Macaca fascicularis",
                      mus_musculus = "Mus musculus",
                      rattus_norvegicus = "Rattus norvegicus",
                      sus_scrofa = "Sus scrofa",
                      drosophila_melanogaster = "Drosophila melanogaster"
    )
  )

# ============================================================
# 3. Recode pathway -> hallmark
# ============================================================

pathway_map <- c(
  "degradation of proteolytic systems"                       = "Loss of proteostasis",
  "impairment of proteins folding and stability"              = "Loss of proteostasis",
  "impairment of the mitochondrial integrity and biogenesis"  = "Mitochondrial dysfunction",
  "accumulation of reactive oxygen species"                   = "Mitochondrial dysfunction",
  "sterile inflammation"                                      = "Chronic inflammation",
  "nuclear DNA instability"                                   = "Genome instability",
  "senescent cells accumulation"                              = "Cellular senescence",
  "changes in the extracellular matrix structure"             = "Extracellular matrix changes",
  "chromatin remodeling"                                      = "Epigenetic alterations",
  "alterations in histone modifications"                      = "Epigenetic alterations",
  "transcriptional alterations"                               = "Epigenetic alterations",
  "intercellular communication impairment"                    = "Intercellular communication",
  "disabled macroautophagy"                                   = "Disabled macroautophagy",
  "stem cell exhaustion"                                      = "Stem cell exhaustion"
)

aging_hallmarks <- c(
  "Genome instability", "Telomere attrition", "Epigenetic alterations",
  "Loss of proteostasis", "Disabled macroautophagy", "Deregulated nutrient sensing",
  "Mitochondrial dysfunction", "Cellular senescence", "Stem cell exhaustion",
  "Intercellular communication", "Chronic inflammation", "Dysbiosis",
  "Extracellular matrix changes", "Psychosocial isolation"
)

all_data <- all_data %>%
  mutate(
    pathway = recode(pathway, !!!pathway_map),
    pathway = factor(pathway, levels = aging_hallmarks)
  )

# ============================================================
# 4. GEO id + contrast from filename (for labels and contrast color)
# ============================================================

extract_geo_id <- function(name) str_extract(name, "^[A-Za-z]+[0-9]+")

extract_contrast <- function(name) {
  m <- str_extract(name, "SIRT6-[A-Za-z0-9]+(_[A-Za-z0-9]+)?_vs_WT")
  m <- gsub("SIRT6-", "SIRT6 ", m)
  m <- gsub("_vs_WT", " vs WT", m)
  m <- gsub("_", " ", m)
  m
}

all_data <- all_data %>%
  mutate(
    geo_id   = extract_geo_id(experiment),
    contrast = extract_contrast(experiment)
  )

# ============================================================
# 5. biological_system from summary_df — simple fuzzy match
#    by experiment_id occurring within experiment
# ============================================================

match_system <- function(exp_name) {
  candidates <- summary_df %>% dplyr::filter(str_detect(exp_name, fixed(experiment_id)))
  
  if (nrow(candidates) == 0) return(NA_character_)
  if (nrow(candidates) == 1) return(candidates$biological_system[1])
  
  # multiple candidates (different strata of same experiment_id) —
  # refine by matching stratum in the experiment name
  refined <- candidates %>% dplyr::filter(str_detect(exp_name, fixed(stratum)))
  if (nrow(refined) >= 1) return(refined$biological_system[1])
  
  candidates$biological_system[1]  # fallback
}

all_data <- all_data %>%
  rowwise() %>%
  mutate(biological_system = match_system(experiment)) %>%
  ungroup()

all_data <- all_data %>%
  mutate(
    biological_system = case_when(
      experiment == "GSE130690_GSE130692_all_SIRT6-KO_vs_WT" ~ "Developmental",
      TRUE ~ biological_system
    )
  )

# ---- Diagnostics ----
unmatched <- all_data %>% dplyr::filter(is.na(biological_system)) %>% distinct(experiment)
if (nrow(unmatched) > 0) {
  cat("Warning: no biological_system found for:\n")
  print(unmatched, n = 50)
} else {
  cat("All experiments matched to a biological_system\n")
}

cat("\nDistribution by system:\n")
print(table(all_data$biological_system, useNA = "always"))

# ============================================================
# 6. Row labels: GEO id + two colored squares (contrast, system)
# ============================================================

label_df <- all_data %>%
  distinct(experiment, geo_id, contrast, biological_system, organism) %>%
  mutate(
    contrast_col = coalesce(contrast_colors[contrast], "grey50"),
    system_col   = coalesce(system_colors[biological_system], "grey70"),
    label = paste0(
      geo_id, " ",
      "<span style='color:", contrast_col, "; padding:6px;'>&#9632;</span>",
      "<span style='color:", system_col, "; padding:6px;'>&#9632;</span>"
    )
  ) %>%
  arrange(organism, geo_id, contrast, experiment)

label_vec <- setNames(label_df$label, as.character(label_df$experiment))

all_data <- all_data %>%
  mutate(experiment = factor(experiment, levels = rev(label_df$experiment)))

ggplot(all_data, aes(x = pathway, y = experiment, size = logFDR)) +
  geom_point(aes(fill = NES), shape = 21, color = "black", stroke = 0.4, alpha = 0.9) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits=c(-3.1,3.1)) +
  scale_y_discrete(labels = label_vec) +
  facet_grid(organism ~ ., scales = "free_y", space = "free_y") +
  theme_bw() +
  theme(
    axis.text.x   = element_text(color = 'black', angle = 60, hjust = 1, size=11),
    axis.text.y   = ggtext::element_markdown(size=10.5, color='black'),
    strip.text.y  = element_text(color = 'black', angle = 0),
    panel.spacing = unit(0.8, "lines")
  ) +
  labs(x = "", y = "", fill = "NES", size = expression(paste(-log[10], "FDR")))

ggsave(
  filename = "GSEA_bubble_plot.pdf",
  plot     = last_plot(),  
  device   = cairo_pdf
)


# Fraction of experiments with negative vs positive NES per hallmark
counts <- all_data %>%
  group_by(pathway) %>%
  summarize(
    n_total    = n(),
    n_negative = sum(NES < 0),
    n_positive = sum(NES > 0),
    pct_negative = round(100 * n_negative / n_total, 1)
  ) %>%
  arrange(desc(n_total)) %>%
  rowwise() %>%
  mutate(p_binom = binom.test(n_negative, n_total, 0.5)$p.value) %>%
  ungroup()
print(counts)

# Check pseudoreplication for the strongest hallmark (Mitochondrial dysfunction)
mito <- all_data %>%
  filter(pathway == "Mitochondrial dysfunction") %>%
  distinct(experiment, geo_id, organism, NES) %>%
  arrange(NES)
print(mito, n = Inf)

mito_by_gse <- mito %>%
  group_by(geo_id) %>%
  summarize(mean_NES = mean(NES), n_tissues = n(), .groups = "drop")
print(mito_by_gse)

bt <- binom.test(sum(mito_by_gse$mean_NES < 0), nrow(mito_by_gse), 0.5)
cat(sprintf("By independent experiments: %d/%d, p=%.4f\n",
            sum(mito_by_gse$mean_NES < 0), nrow(mito_by_gse), bt$p.value))